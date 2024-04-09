#ifndef GEOM_H_GUARD
#define GEOM_H_GUARD
#include <inttypes.h> 
#include <stdbool.h>
#include "vectors.h"
#include "bsdf.h"

struct Ray{
  float tmin;
  struct Vector3f origin;
  struct Vector3f direction;
};

struct Sphere{
  float radius;
  struct Vector3f center;
  struct BSDF bsdf;
  struct Vector3f emission;
};

struct AreaLight{
  float xmin,zmin,xmax,zmax;
  float y;
  struct Vector3f emission_intensity;
};

static inline struct Vector3f warp_square_to_area_light(struct AreaLight* light,struct Vector2f sample){
  struct Vector3f res;
  res.x = light->xmin + (light->xmax-light->xmin) * sample.x;
  res.z = light->zmin + (light->zmax-light->zmin) * sample.y;
  res.y = light->y;
  return res;
}

struct Intersection{
  bool hit;
  float t;
  struct Vector3f point;
  struct Vector3f normal;
  struct BSDF bsdf;
  struct Vector3f emission;
};

static inline struct Intersection intersect_ray_sphere(struct Ray* ray,struct Sphere* sphere){
  struct Intersection it;
  struct Vector3f o = sub_v3f(ray->origin,sphere->center);
  struct Vector3f d = ray->direction;
  //r*r = <o + t*d,o+t*d> = <o,o> + 2*t*<d,o> + t*t<d,d>;
  float c = dot_v3f(o,o) - sphere->radius*sphere->radius; 
  float b = 2 * dot_v3f(d,o);
  float a = dot_v3f(d,d);
  float r = b*b - 4 * a * c;
  if(r < 0.0f){
    it.hit = false;
    it.t = 0;
    return it;
  }
  float root = sqrtf(r);
  float tp = (-b + root) / (2*a);
  float tm = (-b - root) / (2*a);
  it.t = tm >= ray->tmin ? tm : tp;
  it.hit = it.t >= ray->tmin;
  it.point = add_v3f(ray->origin,smul_v3f(it.t,ray->direction));
  it.normal = normalize_v3f(sub_v3f(it.point,sphere->center));
  it.bsdf = sphere->bsdf;
  it.emission = sphere->emission;
  return it;
}

static inline struct Intersection intersect_ray_y_plane(struct Ray* ray,float y){
  struct Intersection it;
  it.normal = (struct Vector3f){0.0f,1.0f,0.0f};
  if(ray->direction.y == 0){
    it.hit = false;
    it.t = 0;
    it.point = ray->origin;
  }else{
    it.t = (y-ray->origin.y) / ray->direction.y;
    it.point = add_v3f(ray->origin , smul_v3f(it.t , ray->direction));
    it.hit = it.t >= ray->tmin;
  }
  return it;
}

static inline struct Intersection intersect_ray_area_light(struct Ray* ray,struct AreaLight* area_light,bool debug){
  struct Intersection it;
  it.normal = (struct Vector3f){0.0f,ray->origin.y < area_light->y ? -1.0f:1.0f,0.0f};
  //if(ray->direction.y == 0){
    //it.hit = false;
    //it.t = 0;
    //it.point = ray->origin;
  //}else{
    it.t = (area_light->y-ray->origin.y) / ray->direction.y;
    it.point = add_v3f(ray->origin, smul_v3f(it.t, ray->direction));
    
    it.hit =  (area_light->xmin <= it.point.x) && (it.point.x <= area_light->xmax) && 
              (area_light->zmin <= it.point.z) && (it.point.z <= area_light->zmax) &&
              it.t >= ray->tmin;
    it.emission = area_light->emission_intensity;
    it.bsdf.type = BSDF_ZERO;
  //}
  return it;
}

struct Frame{
  struct Vector3f n,tu,tv;
};

static inline struct Frame frame_from_normal(struct Vector3f normal){
  struct Frame res;
  res.n = normalize_v3f(normal);
  if(fabs(res.n.x) >= fabs(res.n.y)){
    float len_2d = sqrtf(res.n.x*res.n.x + res.n.z * res.n.z);
    res.tu = (struct Vector3f){res.n.z/len_2d,0.0f,-res.n.x/len_2d};
  }else{
    float len_2d = sqrtf(res.n.y*res.n.y + res.n.z * res.n.z);
    res.tu = (struct Vector3f){0.0f,res.n.z/len_2d,-res.n.y/len_2d};
  }
  res.tv = cross_v3f(res.n,res.tu);
  return res;
}

static inline struct Vector3f frame_to_local(struct Frame* frame,struct Vector3f v){
  struct Vector3f res;
  res.x = dot_v3f(frame->tu,v);
  res.y = dot_v3f(frame->tv,v);
  res.z = dot_v3f(frame->n,v);
  return res;
}

static inline struct Vector3f frame_to_global(struct Frame* frame,struct Vector3f v){
  struct Vector3f res;
  res = smul_v3f(v.x,frame->tu);
  res = add_v3f(res,smul_v3f(v.y,frame->tv));
  res = add_v3f(res,smul_v3f(v.z,frame->n));
  return res;
}

#endif
