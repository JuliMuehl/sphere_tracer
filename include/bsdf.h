#ifndef BSDF_H_GUARD
#define BSDF_H_GUARD
#include "vectors.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

enum BSDFType{BSDF_ZERO,BSDF_DIFFUSE,BSDF_MIRROR,BSDF_COOK_TORRANCE,BSDF_DIELECTRIC};

struct BSDF{
  enum BSDFType type;
  union{
    struct {
      struct Vector3f albedo;
    } diffuse;
    struct {
      struct Vector3f eta;
      struct Vector3f kappa;
      float alpha;
    }cook_torrance;
    struct {
      float eta_inside;
      float eta_outside;
    }dielectric;
  };
};

struct BSDFQueryRecord{
  struct Vector3f wi;
  struct Vector3f wo;
};

const char* bsdf_type_string(enum BSDFType type);
struct Vector3f eval_bsdf(struct BSDF* bsdf,struct BSDFQueryRecord* rec);
float pdf_bsdf(struct BSDF* bsdf,struct BSDFQueryRecord* rec);
struct Vector3f sample_bsdf(struct BSDF* bsdf,struct BSDFQueryRecord* rec,struct Vector2f sample);

#endif
