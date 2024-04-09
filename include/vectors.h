#ifndef VECTORS_H_GUARD
#define VECTORS_H_GUARD
#include <math.h>
struct Vector3f{
  float x,y,z;
};

struct Vector2f{
  float x,y;
};

static inline float dot_v3f(struct Vector3f v,struct Vector3f w){
  return v.x*w.x + v.y*w.y + v.z*w.z;
}

static inline struct Vector3f smul_v3f(float s,struct Vector3f v){
  struct Vector3f res = {s*v.x,s*v.y,s*v.z};
  return res;
}

static inline struct Vector3f add_v3f(struct Vector3f v,struct Vector3f w){
  struct Vector3f res = {v.x+w.x,v.y+w.y,v.z+w.z};
  return res;
}

static inline float norm_v3f(struct Vector3f v){
  return sqrtf(dot_v3f(v,v));
}

static inline struct Vector3f normalize_v3f(struct Vector3f v){
  return smul_v3f(1.0f / norm_v3f(v),v);
}

static inline struct Vector3f sub_v3f(struct Vector3f v,struct Vector3f w){
  struct Vector3f res = {v.x-w.x,v.y-w.y,v.z-w.z};
  return res;
}

static inline struct Vector3f cross_v3f(struct Vector3f v,struct Vector3f w){
  struct Vector3f res;
  res.x = v.y*w.z - v.z*w.y;
  res.y = v.z*w.x - v.x*w.z;
  res.z = v.x*w.y - v.y*w.x;
  return res;
}

#endif
