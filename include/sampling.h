#ifndef SAMPLING_H_GUARD
#define SAMPLING_H_GUARD
#include <math.h>
#include "vectors.h"

static inline struct Vector3f warp_square_to_uniform_hemisphere(struct Vector2f sample){
  float cos_theta = sample.x;
  float phi = 2.0f * M_PI * sample.y;
  float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);
  struct Vector3f res = {
    cosf(phi) * sin_theta,
    sinf(phi) * sin_theta,
    cos_theta
  };
  return res;
}

static inline struct Vector3f warp_square_to_cosine_hemisphere(struct Vector2f sample){
  float cos_theta = sqrtf(1.0f-sample.x);
  float phi = 2.0f * M_PI * sample.y;
  float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);
  struct Vector3f res = {
    cosf(phi) * sin_theta,
    sinf(phi) * sin_theta,
    cos_theta
  };
  return res;
}

static inline struct Vector2f warp_square_to_unit_gaussian(struct Vector2f sample){
  struct Vector2f res;
  float r = sqrt(-2.0f*logf(1.0f-sample.x));
  float theta = 2.0f*M_PI*sample.y;
  res.x = r*sinf(theta);
  res.y = r*cosf(theta);
  return res;
}
#endif
