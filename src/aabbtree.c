#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <float.h>
#include <string.h>
#include <pthread.h>
#include <time.h>
#include <stdalign.h>
#include <assert.h>

#include "aabbtree.h"

#include "geom.h"

void print_v3f(struct Vector3f v){
  printf("(%f,%f,%f)\n",v.x,v.y,v.z);
}

struct AABB{
  float xmin,ymin,zmin,xmax,ymax,zmax;
};

float min(float a,float b){
  return a <= b ? a:b;
}

float max(float a,float b){
  return a <= b ? b:a;
}

struct AABB aabb_from_triangle(struct Triangle* triangle){
  float xmin = min(min(triangle->p1.x,triangle->p2.x),triangle->p3.x);
  float ymin = min(min(triangle->p1.y,triangle->p2.y),triangle->p3.y);
  float zmin = min(min(triangle->p1.z,triangle->p2.z),triangle->p3.z);
  float xmax = max(max(triangle->p1.x,triangle->p2.x),triangle->p3.x);
  float ymax = max(max(triangle->p1.y,triangle->p2.y),triangle->p3.y);
  float zmax = max(max(triangle->p1.z,triangle->p2.z),triangle->p3.z);
  struct AABB res = {xmin,ymin,zmin,xmax,ymax,zmax};
  return res;
}

struct AABB union_aabb(struct AABB* a,struct AABB* b){
  return (struct AABB){
    min(a->xmin,b->xmin),
    min(a->ymin,b->ymin),
    min(a->zmin,b->zmin),
    max(a->xmax,b->xmax),
    max(a->ymax,b->ymax),
    max(a->zmax,b->zmax),
  };
}

struct Vector3f centroid_aabb(struct AABB* a){
  struct Vector3f res = {
    0.5 * (a->xmin + a->xmax),
    0.5 * (a->ymin + a->ymax),
    0.5 * (a->zmin + a->zmax),
  };
  return res;
}

float volume_aabb(const struct AABB* a){
  return (a->xmax - a->xmin) * (a->ymax - a->ymin) * (a->zmax - a->zmin);
}

float surface_area_aabb(const struct AABB* a){
  struct Vector3f d = {a->xmax - a->xmin,a->ymax-a->ymin,a->zmax-a->zmin};
  return 2.0f * (d.x*d.y + d.x*d.z + d.y*d.z);
}

struct AABBTreeNode{
  struct AABBTreeNode *left,*right;
  struct AABB aabb;
  struct Triangle* triangle;
  uint8_t axis;
};

int _aabb_node_cmp_x(const void* _left,const void* _right){
  const struct AABBTreeNode* left = *(const struct AABBTreeNode**) _left;
  const struct AABBTreeNode* right = *(const struct AABBTreeNode**) _right;
  return (left->aabb.xmin + left->aabb.xmax) - (right->aabb.xmin + right->aabb.xmax);
}

int _aabb_node_cmp_y(const void* _left,const void* _right){
  const struct AABBTreeNode* left = *(const struct AABBTreeNode**) _left;
  const struct AABBTreeNode* right = *(const struct AABBTreeNode**) _right;
  return (left->aabb.ymin + left->aabb.ymax) - (right->aabb.ymin + right->aabb.ymax);
}

int _aabb_node_cmp_z(const void* _left,const void* _right){
  const struct AABBTreeNode* left = *(const struct AABBTreeNode**) _left;
  const struct AABBTreeNode* right = *(const struct AABBTreeNode**) _right;
  return ((left->aabb.zmin + left->aabb.zmax) - (right->aabb.zmin + right->aabb.zmax));
}

int (*_aabb_node_cmp_funcs[3])(const void*,const void*) = {_aabb_node_cmp_x,_aabb_node_cmp_y,_aabb_node_cmp_z};

struct AABBTreeNode* build_aabbtree_from_leaf_nodes(struct AABBTreeNode** leaf_nodes,size_t slice_start,size_t slice_end){
  if(slice_start >= slice_end){
    return NULL;
  }
  if(slice_start == slice_end - 1){
    return leaf_nodes[slice_start];
  }
  int axis = rand()%3;
  qsort(&leaf_nodes[slice_start],
        slice_end-slice_start,
        sizeof(struct AABBTreeNode*),
        _aabb_node_cmp_funcs[axis]);
  struct AABBTreeNode* node = malloc(sizeof(struct AABBTreeNode));
  node->axis = axis;
  node->left = build_aabbtree_from_leaf_nodes(leaf_nodes,slice_start,slice_start + (slice_end-slice_start) / 2);
  node->right = build_aabbtree_from_leaf_nodes(leaf_nodes,slice_start + (slice_end-slice_start) / 2,slice_end);
  node->aabb = leaf_nodes[slice_start]->aabb;
  for(int i = slice_start+1;i<slice_end;i++){
    node->aabb = union_aabb(&node->aabb,&leaf_nodes[i]->aabb);
  }
  node->triangle = NULL;
  return node;
}

struct AABBTreeNode* build_aabbtree_from_leaf_nodes_sah(struct AABBTreeNode** leaf_nodes,size_t slice_start,size_t slice_end){
  //printf("build_aabbtree_from_leaf_nodes_sah(%p,%lu,%lu) called\n",leaf_nodes,slice_start,slice_end);
  if(slice_start >= slice_end){
    return NULL;
  }
  if(slice_start == slice_end - 1){
    return leaf_nodes[slice_start];
  }
  struct AABB bounds = {0.0f};
  for(size_t i = slice_start;i<slice_end;i++){
    bounds = union_aabb(&bounds,&leaf_nodes[i]->aabb);
  }
  float sizes[3] = {bounds.xmax - bounds.xmin,bounds.ymax - bounds.ymin,bounds.zmax - bounds.zmin};
  int axis = 0;
  for(int j = 1;j<3;j++){
    if(sizes[axis] <= sizes[j])
      axis = j;
  }
  int nBuckets = 12;
  struct BucketInfo {
     int count;
     struct AABB bounds;
  };
  struct BucketInfo buckets[nBuckets];
  memset(&buckets[0],0,sizeof(buckets));
  float mins[3] = {bounds.xmin,bounds.ymin,bounds.zmin};
  for (size_t i = slice_start; i < slice_end; i++) {
    struct Vector3f c = centroid_aabb(&leaf_nodes[i]->aabb);
    float cvec[3] = {c.x,c.y,c.z};
    int b = nBuckets * (cvec[axis] - mins[axis]) / sizes[axis];
    if (b == nBuckets) b = nBuckets - 1;
    buckets[b].count++;
    buckets[b].bounds = union_aabb(&buckets[b].bounds, &leaf_nodes[i]->aabb);
  }
  float cost[nBuckets - 1];
  //printf("Bounds surface area = %f\n",surface_area_aabb(&bounds));
  for (int i = 0; i < nBuckets - 1; ++i) {
    struct AABB b0 = {0.0f}, b1 = {0.0f};
    int count0 = 0, count1 = 0;
    for (int j = 0; j <= i; ++j) {
       b0 = union_aabb(&b0, &buckets[j].bounds);
       count0 += buckets[j].count;
    }
    for (int j = i+1; j < nBuckets; ++j) {
       b1 = union_aabb(&b1, &buckets[j].bounds);
       count1 += buckets[j].count;
    }
    //printf("b0={(%f,%f,%f),(%f,%f,%f)}\n",b0.xmin,b0.ymin,b0.zmin,b0.xmax,b0.ymax,b0.zmax);
    //printf("b1={(%f,%f,%f),(%f,%f,%f)}\n",b0.xmin,b0.ymin,b0.zmin,b0.xmax,b0.ymax,b0.zmax);
    cost[i] = .125f + (count0 * surface_area_aabb(&b0) +
                      count1 * surface_area_aabb(&b1)) / surface_area_aabb(&bounds);
  }

  float minCost = cost[0];
  int minCostSplitBucket = 0;
  for (int i = 1; i < nBuckets - 1; ++i) {
    if (cost[i] < minCost) {
      minCost = cost[i];
      minCostSplitBucket = i;
    }
  }

  size_t i = slice_start,j=slice_end-1;
  while(i < j){
    struct Vector3f c = centroid_aabb(&leaf_nodes[i]->aabb);
    float cvec[3] = {c.x,c.y,c.z};
    int b = nBuckets * (cvec[axis] - mins[axis]) / sizes[axis];
    if (b == nBuckets) b = nBuckets - 1;
    if(b <= minCostSplitBucket){
      i++;
    }else{
      struct AABBTreeNode* tmp = leaf_nodes[i];
      leaf_nodes[i] = leaf_nodes[j];
      leaf_nodes[j] = tmp;
      j--;
    }
  }
  //Fall back to splitting in the middle in case we don't get any nodes on one side
  //TODO we might replace this by just creating a leaf node with multiple primitives
  if(i == slice_start || i == slice_end-1){
    return build_aabbtree_from_leaf_nodes(leaf_nodes,slice_start,slice_end);
  }
  struct AABBTreeNode* node = malloc(sizeof(struct AABBTreeNode));
  //printf("i=%lu\n",i);
  node->axis = axis;
  node->left = build_aabbtree_from_leaf_nodes_sah(leaf_nodes,slice_start,i);
  node->right = build_aabbtree_from_leaf_nodes_sah(leaf_nodes,i,slice_end);
  
  node->aabb = bounds;
  node->triangle = NULL;
  return node;
}

void destroy_aabbtree_recursive(struct AABBTreeNode* node){
  if(node == NULL) return;
  destroy_aabbtree_recursive(node->left);
  destroy_aabbtree_recursive(node->right);
  free(node);
}

struct Intersection intersect_ray_aabb_fast(struct AABB* aabb,struct Ray* ray,float tmax){
  struct Intersection res;
  float max_values[3] = {aabb->xmax,aabb->ymax,aabb->zmax};
  float min_values[3] = {aabb->xmin,aabb->ymin,aabb->zmin};
  float d[3] = {ray->direction.x,ray->direction.y,ray->direction.z};
  float o[3] = {ray->origin.x,ray->origin.y,ray->origin.z};
  res.t = tmax;
  int res_axis = 0;
  float t0=ray->tmin,t1=tmax;
  for(int i = 0;i<3;i++){
    float t_near = (min_values[i] - o[i]) / d[i];
    float t_far = (max_values[i] - o[i]) / d[i];
    if(t_near > t_far){
      float tmp = t_near;t_near = t_far;t_far=tmp;
    }
    t0 = t_near > t0 ? t_near : t0;
    t1 = t_far < t1 ? t_far  : t1;
    if(t0 > t1){
      res.hit = false;
      return res;
    }
  }
  res.hit = true;
  res.t = t0;
  return res;
}

float signf(float x){
  if(x >= 0)
    return 1;
  return -1;
}

struct Intersection intersect_ray_triangle(struct Triangle* triangle,struct Ray* ray){
  struct Intersection res;
  float t = -(dot_v3f(ray->origin,triangle->normal) - dot_v3f(triangle->p1,triangle->normal))/dot_v3f(ray->direction,triangle->normal);
  res.t = t;
  res.point = add_v3f(ray->origin,smul_v3f(t,ray->direction));
  struct Vector3f n1 = normalize_v3f(cross_v3f(sub_v3f(triangle->p1,triangle->p2),triangle->normal));
  struct Vector3f n2 = normalize_v3f(cross_v3f(sub_v3f(triangle->p1,triangle->p3),triangle->normal));
  struct Vector3f n3 = normalize_v3f(cross_v3f(sub_v3f(triangle->p2,triangle->p3),triangle->normal));
  float s1 = signf(dot_v3f(sub_v3f(triangle->p3,triangle->p1),n1));
  float s2 = signf(dot_v3f(sub_v3f(triangle->p2,triangle->p1),n2));
  float s3 = signf(dot_v3f(sub_v3f(triangle->p1,triangle->p2),n3));
  float c1 = s1 * dot_v3f(sub_v3f(res.point,triangle->p1),n1);
  float c2 = s2 * dot_v3f(sub_v3f(res.point,triangle->p1),n2);
  float c3 = s3 * dot_v3f(sub_v3f(res.point,triangle->p2),n3);
  res.hit = 0.0f <= c1  && 
            0.0f <= c2  && 
            0.0f <= c3  && 
            ray->tmin <= t;
  res.normal = triangle->normal;
  res.point = add_v3f(ray->origin,smul_v3f(res.t,ray->direction));
  return res;
}

struct AABBTreeNode* build_aabbtree_from_triangle_array(struct TriangleArray triangles){
  struct AABBTreeNode** leaf_nodes = malloc(sizeof(struct AABBTreeNode*) * triangles.len);
  for(size_t i = 0;i<triangles.len;i++){
    leaf_nodes[i] = malloc(sizeof(struct AABBTreeNode));
    leaf_nodes[i]->aabb = aabb_from_triangle(&triangles.data[i]);
    leaf_nodes[i]->triangle = &triangles.data[i];
    leaf_nodes[i]->left = NULL;
    leaf_nodes[i]->right = NULL;
  }
  srand(0);
  struct AABBTreeNode* root = build_aabbtree_from_leaf_nodes_sah(leaf_nodes,0,triangles.len);
  return root;
}

struct LinearAABBNode{
  alignas(32) struct AABB aabb;
  alignas(32) struct {
    union{
      unsigned int triangle_offset:29;
      unsigned int second_child_offset:29;
    };
    uint8_t axis:2;
    bool is_leaf:1;
  };
};

struct LinearAABBTree{
  struct LinearAABBNode* arr;
  size_t len,cap;
  struct Triangle* triangles;
  struct Frame frame;
  struct Vector3f position;
};

size_t flatten_aabbtree_recursive(struct LinearAABBTree* tree,struct AABBTreeNode* node,struct Triangle* triangles){
  while(tree->len >= tree->cap){
    tree->cap *= 2;
    tree->arr = reallocarray(tree->arr,tree->cap,sizeof(tree->arr[0]));
  }
  tree->arr[tree->len].aabb = node->aabb;
  tree->arr[tree->len].axis = node->axis;
  if(node->triangle != NULL){
    tree->arr[tree->len].is_leaf = true;
    tree->arr[tree->len].triangle_offset = node->triangle - triangles;
    return tree->len++;
  }else{
    assert(node->left != NULL && node->right != NULL);
    size_t node_idx = tree->len;
    tree->arr[node_idx].is_leaf = false;
    tree->len++;
    flatten_aabbtree_recursive(tree,node->left,triangles);
    size_t right_index = flatten_aabbtree_recursive(tree,node->right,triangles);
    tree->arr[node_idx].second_child_offset = right_index;
    return node_idx;
  }
}

struct LinearAABBTree* flatten_aabbtree(struct AABBTreeNode* root,struct Triangle* triangles){
  struct LinearAABBTree* res = malloc(sizeof(*res));
  res->triangles = triangles;
  res->cap = 16;
  res->arr = malloc(sizeof(res->arr[0]) * res->cap);
  res->len = 0;
  flatten_aabbtree_recursive(res,root,triangles);
  res->cap = res->len;
  res->arr = reallocarray(res->arr,res->cap,sizeof(res->arr[0]));
  return res;
}

struct Intersection intersect_ray_linear_aabbnode(struct LinearAABBTree* tree,size_t node_idx,float t,struct Ray* ray){
  struct LinearAABBNode* node = &tree->arr[node_idx];
  if(node->is_leaf){
    struct Intersection res = intersect_ray_triangle(&tree->triangles[node->triangle_offset],ray);
    return res;
  }
  struct Intersection res_left;
  struct Intersection res_right;
  res_left.hit = false;
  res_right.hit = false;
  size_t left_idx = node_idx+1;
  size_t right_idx = node->second_child_offset;
  struct Intersection it_left_aabb = intersect_ray_aabb_fast(&tree->arr[left_idx].aabb,ray,t);
  struct Intersection it_right_aabb = intersect_ray_aabb_fast(&tree->arr[right_idx].aabb,ray,t);
  if(it_left_aabb.hit && it_right_aabb.hit && it_left_aabb.t > it_right_aabb.t){
    size_t tmp = left_idx;left_idx=right_idx;right_idx=tmp;
  }
  if(it_left_aabb.hit && it_left_aabb.t <= t){
    res_left = intersect_ray_linear_aabbnode(tree,left_idx,t,ray);
  }else{
    res_left.hit = false;
  }
  if(res_left.hit)
    t = res_left.t;
  if(it_right_aabb.hit && it_right_aabb.t <= t){
    res_right = intersect_ray_linear_aabbnode(tree,right_idx,t,ray);
  }else{
    res_right.hit = false;
  }
  if(!res_left.hit){
    return res_right;
  }
  if(!res_right.hit){
    return res_left;
  }
  if(res_left.t <= res_right.t){
    return res_left;
  }
  return res_right;
}

struct Intersection intersect_ray_linear_aabbtree(struct LinearAABBTree* tree,struct Ray* ray){
  float t = FLT_MAX;
  struct Intersection aabb_intersection = intersect_ray_aabb_fast(&tree->arr[0].aabb,ray,t);
  if(aabb_intersection.hit)
    return intersect_ray_linear_aabbnode(tree,0,t,ray);
  else
    return aabb_intersection;
}

AABBTree* create_aabbtree(struct TriangleArray triangles){
  struct AABBTreeNode* root = build_aabbtree_from_triangle_array(triangles);
  struct LinearAABBTree* tree = flatten_aabbtree(root,triangles.data);
  tree->frame.tu = (struct Vector3f){1.0f,0.0f,0.0f};
  tree->frame.tv = (struct Vector3f){0.0f,1.0f,0.0f};
  tree->frame.n = (struct Vector3f){0.0f,0.0f,1.0f};
  tree->position = (struct Vector3f){0.0f,0.0f,0.0f};
  destroy_aabbtree_recursive(root);
  return (AABBTree*)tree;
}

void set_aabbtree_frame(AABBTree* aabbtree,struct Frame frame){
  struct LinearAABBTree* tree = (struct LinearAABBTree*) aabbtree;
  tree->frame = frame;
}
void set_aabbtree_position(AABBTree* aabbtree,struct Vector3f position){
  struct LinearAABBTree* tree = (struct LinearAABBTree*) aabbtree;
  tree->position = position;
}

struct Intersection intersect_ray_aabbtree(struct Ray* ray,AABBTree* aabbtree){
  struct LinearAABBTree* tree = (struct LinearAABBTree*) aabbtree;
  struct Ray transformed_ray;
  transformed_ray.tmin = ray->tmin;
  transformed_ray.origin = sub_v3f(ray->origin,tree->position);
  transformed_ray.direction = frame_to_local(&tree->frame, ray->direction);
  struct Intersection it = intersect_ray_linear_aabbtree(tree,&transformed_ray);
  it.normal = it.normal;
  it.point = add_v3f(tree->position,it.point);
  return it;
}
