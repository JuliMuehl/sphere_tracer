#ifndef WORLD_GUARD_H
#define WORLD_GUARD_H
#include "aabbtree.h"
#include "geom.h"

struct World{
  uint32_t num_spheres;
  struct Sphere* spheres;
  float plane_y;
  struct BSDF plane_bsdf,aabbtree_bsdf;
  struct AreaLight light;
  struct Vector3f sky_color;
  AABBTree* aabbtree;
};

static inline struct Intersection intersect_ray_world(struct Ray* ray,struct World* world){
  struct Intersection it = intersect_ray_y_plane(ray,world->plane_y);
  it.emission = (struct Vector3f){0.0f,0.0f,0.0f};
  it.bsdf = world->plane_bsdf;
  struct Intersection light_it = intersect_ray_area_light(ray,&world->light,false);
  if(!it.hit || light_it.hit && light_it.t < it.t ){
    it = light_it;
  }
  for(uint32_t i = 0;i<world->num_spheres;i++){
    struct Intersection sphere_it = intersect_ray_sphere(ray,&world->spheres[i]);
    if(sphere_it.hit){
      if(!it.hit || sphere_it.t < it.t){
        it = sphere_it;
      }
    }
  }
  struct Intersection aabbtree_it = intersect_ray_aabbtree(ray,world->aabbtree);
  if(aabbtree_it.hit){
    if(!it.hit || aabbtree_it.t < it.t){
      it = aabbtree_it;
      //it.normal = smul_v3f(-1.0f,it.normal);
      //printf("it.normal = %f %f %f\n",it.normal.x,it.normal.y,it.normal.z);
      //printf("it.point = %f %f %f\n",it.point.x,it.point.y,it.point.z);
      it.bsdf = world->aabbtree_bsdf;
      it.emission = (struct Vector3f){0.0f,0.0f,0.0f};
    }
  }
  return it;
}
#endif
