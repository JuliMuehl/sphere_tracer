#include "kdtree.h"
int cmp_x(const void* v_,const void* w_){
  const struct Photon* v = (const struct Photon*)v_,*w = (const struct Photon*)w_;
  return v->point.x > w->point.x;
}

int cmp_y(const void* v_,const void* w_){
  const struct Photon* v = (const struct Photon*)v_,*w = (const struct Photon*)w_;
  return v->point.y > w->point.y;
}

int cmp_z(const void* v_,const void* w_){
  const struct Photon* v = (const struct Photon*)v_,*w = (const struct Photon*)w_;
  return v->point.z > w->point.z;
}

struct KDTreeNode* build_kdtree_recursive(struct Photon* photons,uint32_t num_photons,uint32_t depth){
  if(num_photons == 0) return NULL;
  struct KDTreeNode* node = (struct KDTreeNode*)malloc(sizeof(*node));
  uint32_t axis = depth % 3;
  __compar_fn_t cmp_funcs[3] = {cmp_x,cmp_y,cmp_z};
  qsort(photons,num_photons,sizeof(*photons),cmp_funcs[axis]);
  node->axis = axis;
  node->photon = photons[num_photons/2];
  if(num_photons > 1){ // => num_photons/2 >= 1
    node->left = build_kdtree_recursive(photons,num_photons/2,depth+1);
    node->right = build_kdtree_recursive(photons+num_photons/2+1,num_photons - num_photons/2 - 1,depth+1);
  }else{
    node->left = NULL;
    node->right = NULL;
  }
  return node;
}

bool element_of_aabb(struct Vector3f point,const struct AABB* b){
  return b->xmin <= point.x && point.x <= b->xmax &&
         b->ymin <= point.y && point.y <= b->ymax && 
         b->zmin <= point.z && point.z <= b->zmax;
}

struct AABB left_node_aabb(const struct AABB* parent,float split_value,uint32_t splitting_axis){
  float bounds[] = {parent->xmin,parent->xmax,
                    parent->ymin,parent->ymax,
                    parent->zmin,parent->zmax};
  bounds[2*splitting_axis + 1] = split_value;
  struct AABB res = {
    .xmin = bounds[0],.xmax = bounds[1],
    .ymin = bounds[2],.ymax = bounds[3],
    .zmin = bounds[4],.zmax = bounds[5],
  };
  return res;
}

struct AABB right_node_aabb(const struct AABB *parent,float splitting_value,uint32_t splitting_axis){
  float bounds[] = {parent->xmin,parent->xmax,
                    parent->ymin,parent->ymax,
                    parent->zmin,parent->zmax};
  bounds[2*splitting_axis] = splitting_value;
  struct AABB res = {
    .xmin = bounds[0],.xmax = bounds[1],
    .ymin = bounds[2],.ymax = bounds[3],
    .zmin = bounds[4],.zmax = bounds[5],
  };
  return res;
}

float dist_point_aabb(struct Vector3f point,const struct AABB* aabb){
  if(element_of_aabb(point,aabb)) return 0.0f;
  return fmin(fmin(fmin(fabs(point.x-aabb->xmin),fabs(point.x-aabb->xmax)),
                   fmin(fabs(point.y-aabb->ymin),fabs(point.y-aabb->ymax))),
                   fmin(fabs(point.z-aabb->zmin),fabs(point.z-aabb->zmax)));
}

void nn_search_kdtree(struct KDTreeNode* node,
                      struct Vector3f point,
                      struct Photon* best_photon,
                      float* min_dist){
  if(node == NULL) return;
  static uint32_t count = 0;
  count++;
  //printf("count = %u\n",count);
  struct Vector3f root = node->photon.point;
  float root_dist = norm_v3f(sub_v3f(root,point));
  if(root_dist < *min_dist){
    *best_photon = node->photon;
    *min_dist = root_dist; 
  }
  float point_coords[3] = {point.x,point.y,point.z};
  float root_coords[3] = {root.x,root.y,root.z};
  if(node->left && !(point_coords[node->axis] > root_coords[node->axis] 
                      && point_coords[node->axis] - root_coords[node->axis] > *min_dist)){
    nn_search_kdtree(node->left,point,best_photon,min_dist);
  }
  if(node->right &&!(point_coords[node->axis] < root_coords[node->axis] 
                      && point_coords[node->axis] - root_coords[node->axis] < -*min_dist)){
    nn_search_kdtree(node->right,point,best_photon,min_dist);
  }
}

void knn_search_kdtree(struct KDTreeNode* node,
                      struct Vector3f point,
                      uint32_t k,
                      struct Photon* best_photons,
                      float* min_dists){
  struct Vector3f root = node->photon.point;
  float root_dist = norm_v3f(sub_v3f(root,point));
  if(root_dist < min_dists[k-1]){
    min_dists[k-1] = root_dist;
    for(uint32_t i = 0;i<k;i++){
      if(root_dist < min_dists[i]){
        for(uint32_t j = k-1;j>i;j--){
          best_photons[j] = best_photons[j-1];
          min_dists[j] = min_dists[j-1];
        }
        best_photons[i] = node->photon;
        min_dists[i] = root_dist; 
        break;
      }
    }
  }
  float point_coords[3] = {point.x,point.y,point.z};
  float root_coords[3] = {root.x,root.y,root.z};
  float min_dist = min_dists[k-1];
  if(node->left && !(point_coords[node->axis] > root_coords[node->axis] 
                      && point_coords[node->axis] - root_coords[node->axis] > min_dist)){
    knn_search_kdtree(node->left,point,k,best_photons,min_dists);
  }
  if(node->right &&!(point_coords[node->axis] < root_coords[node->axis] 
                      && point_coords[node->axis] - root_coords[node->axis] < -min_dist)){
    knn_search_kdtree(node->right,point,k,best_photons,min_dists);
  }
}
