#include <stdbool.h>
#include "vectors.h"
#include "sampling.h"
#include "bsdf.h"

const char* bsdf_type_strings[] = {"BSDF_ZERO",
                                   "BSDF_DIFFUSE",
                                   "BSDF_MIRROR",
                                   "BSDF_COOK_TORRANCE",
                                   "BSDF_DIELECTRIC"};

const char* bsdf_type_string(enum BSDFType type){
  return bsdf_type_strings[type];
}

float beckmann_density(struct Vector3f wh,float alpha){
  float cos_theta = wh.z;
  float mx=wh.x/wh.z,my=wh.y/wh.z;
  if(cos_theta <= 0.0f) return 0.0f; // lower hemisphere has probability mass 0
  float cos_theta_squared = cos_theta*cos_theta;
  float cos_theta_fourth_power = cos_theta_squared*cos_theta_squared;
  float D = expf(-(mx*mx+my*my) / (alpha*alpha)) / (M_PI*alpha*alpha*cos_theta_fourth_power);
  return D;
}

float G1(struct Vector3f wo,struct Vector3f wh,float alpha){
  float cos_theta = wo.z;
  if(dot_v3f(wo,wh) / cos_theta <= 0.0f ) return 0.0f;
  float sin_theta = sqrtf(1.0f - cos_theta*cos_theta);
  float tan_theta = sin_theta / cos_theta;
  float a = 1.0f/(tan_theta * alpha);
  //float lambda = (erf(a) - 1.0f) / 2.0f + 1.0f/(2.0f*a*sqrtf(M_PI)) * expf(-a*a);
  //float lambda = 0.0f;
  float lambda = (a < 1.6) * (1.0f - 1.259f*a + 0.396f*a*a) / (3.535f*a + 2.181f*a*a);
  return 1.0f/(1.0f + lambda);
}

float geometry_term(struct Vector3f wi,struct Vector3f wo,struct Vector3f wh,float alpha){
  return G1(wi,wh,alpha) * G1(wo,wh,alpha);
}

//Fresnel term includes an absorption coefficient in the imaginary part of index of refraction N = eta - i kappa
//We assume the light doesn't change medium when hitting a microfacet however the facet absorbs light and the air/vacuum doesn't
//Thus eta1 = eta2 = 1.0f, kappa1 = 0.0f, kappa2 = bsdf.cook_torrance.kappa2.
float fresnel(float eta1,float kappa1,float eta2,float kappa2,float cos_r){
  _Complex float N1 = eta1 + kappa1 * _Complex_I;
  _Complex float N2 = eta2 + kappa2 * _Complex_I;
  _Complex float cos_t = csqrtf(1 - (N1/N2)*(1-cos_r*cos_r));
  float sqrt_Rs = cabsf((N1*cos_r - N2*cos_t) / (N1*cos_r + N2*cos_t));
  float sqrt_Rp = cabsf((N1*cos_t - N2*cos_r) / (N1*cos_t + N2*cos_r));
  float R = 0.5 * (sqrt_Rs*sqrt_Rs + sqrt_Rp*sqrt_Rp);
  return R >= 1.0f ? 1.0f:R;
}

struct Vector3f fresnel_rgb(struct Vector3f eta1,struct Vector3f kappa1,struct Vector3f eta2,struct Vector3f kappa2,float cos_r){
  struct Vector3f res;
  res.x = fresnel(eta1.x,kappa1.x,eta2.x,kappa2.x,cos_r);
  res.y = fresnel(eta1.y,kappa1.y,eta2.y,kappa2.y,cos_r);
  res.z = fresnel(eta1.z,kappa1.z,eta2.z,kappa2.z,cos_r);
  return res;
}

struct Vector3f eval_bsdf(struct BSDF* bsdf,struct BSDFQueryRecord* rec){
  //if(rec->wi.z <= 0 || rec->wo.z <= 0) return (struct Vector3f){0.0f,0.0f,0.0f};
  struct Vector3f color = {0.0f,0.0f,0.0f};
  struct Vector3f wh = normalize_v3f(add_v3f(rec->wo,rec->wi));
  switch(bsdf->type){
    case BSDF_DIFFUSE:{
      color = bsdf->diffuse.albedo;
      break;
    }
    case BSDF_DIELECTRIC:{
      color = (struct Vector3f){1.0f,1.0f,1.0f};
      break;
    }
    case BSDF_MIRROR:{
      color = (struct Vector3f){1.0f,1.0f,1.0f};
      break;
    }
    case BSDF_COOK_TORRANCE:{
      float D = beckmann_density(wh,bsdf->cook_torrance.alpha);
      struct Vector3f eta1 = {1.0f,1.0f,1.0f};
      struct Vector3f kappa1 = {0.0f,0.0f,0.0f}; 
      struct Vector3f eta2 = bsdf->cook_torrance.eta;
      struct Vector3f kappa2 = bsdf->cook_torrance.kappa;
      float G = geometry_term(rec->wi,rec->wo,wh,bsdf->cook_torrance.alpha);
      struct Vector3f F = fresnel_rgb(eta1,eta2,kappa1,kappa2,rec->wo.z);
      float denom = 4.0f*rec->wo.z*rec->wi.z;
      color.x = G*D*F.x/denom;
      color.y = G*D*F.y/denom;
      color.z = G*D*F.z/denom;
      break;
    }
    case BSDF_ZERO:{
      color = (struct Vector3f){0.0f,0.0f,0.0f};
      break;
    }
  }
  return color;
}

float pdf_bsdf(struct BSDF* bsdf,struct BSDFQueryRecord* rec){
  //if(rec->wi.z <= 0 || rec->wo.z <= 0) return 0.0f;
  float pdf = 0;
  switch(bsdf->type){
    case BSDF_DIFFUSE:{
      pdf = rec->wo.z / M_PI;
      break;
    }
    case BSDF_MIRROR:{
      pdf = 0;
      break;
    }
    case BSDF_COOK_TORRANCE:{
      float alpha = bsdf->cook_torrance.alpha;//(1.2 - 0.2 * sqrt(rec->wi.z)) * bsdf->cook_torrance.alpha;
      float alpha_squared = alpha * alpha;
      struct Vector3f wh = normalize_v3f(add_v3f(rec->wo,rec->wi));
      float mx=wh.x/wh.z,my=wh.y/wh.z;
      pdf = 1/(M_PI*alpha_squared)*expf(-(mx*mx + my*my) / (alpha_squared)) / (wh.z*wh.z*wh.z) / (4.0 * dot_v3f(rec->wi,wh));
      break;
    }
    case BSDF_DIELECTRIC:{
      pdf = 0.0f;
      break;
    }
    case BSDF_ZERO:{
      pdf = 0;
      break;
    }
  }
  return pdf;
}

struct Vector3f sample_bsdf(struct BSDF* bsdf,struct BSDFQueryRecord* rec,struct Vector2f sample){
  struct Vector3f color = {0.0f,0.0f,0.0f};
  switch(bsdf->type){
    case BSDF_DIFFUSE:{
      rec->wo = warp_square_to_cosine_hemisphere(sample);
      color = smul_v3f(fabs(rec->wo.z)/pdf_bsdf(bsdf,rec),eval_bsdf(bsdf,rec));
    break;
    }
    case BSDF_MIRROR:{
      rec->wo = rec->wi;
      rec->wo.x *= -1.0f;
      rec->wo.y *= -1.0f;
      color = (struct Vector3f){1.0f,1.0f,1.0f};
      break;
    }
    case BSDF_DIELECTRIC:{
      float eta1 = bsdf->dielectric.eta_outside,eta2 = bsdf->dielectric.eta_inside;
      if(rec->wi.z < 0){
        eta1 = bsdf->dielectric.eta_inside,eta2 = bsdf->dielectric.eta_outside;
      }
      
      float cos_i = fabs(rec->wi.z);
      float cos_t = sqrtf(1.0f - (eta1 / eta2) * (eta1 / eta2) * (1-cos_i*cos_i));
      float rs = (eta1*cos_i - eta2*cos_t) / (eta1*cos_i + eta2*cos_t);
      float rp = (eta1*cos_t - eta2*cos_i) / (eta1*cos_t + eta2*cos_i);
      float fresnel = 0.5f * (rs*rs+rp*rp);
      if (eta1 > eta2)
      {
        const float sin_total_internal_reflection = eta2 / eta1;
        const float sin_i = sqrt(1.f - cos_i * cos_i);
        if (sin_i > sin_total_internal_reflection)
          fresnel = 1.0f;
      }
      //printf("%f\n",fresnel);
      bool reflecting = sample.x <= fresnel;
      //reflecting = true;
      if(reflecting){
        rec->wo = rec->wi;
        rec->wo.x *= -1.0f;
        rec->wo.y *= -1.0f;
      }else{
        rec->wo = smul_v3f(-eta1/eta2,rec->wi);
        if(rec->wi.z > 0){
          rec->wo.z += eta1/eta2*cos_i - cos_t;
        }else{
          rec->wo.z -= eta1/eta2*cos_i - cos_t;
        }
      }
      color = (struct Vector3f){1.0f,1.0f,1.0f};
      break;
    }
    case BSDF_COOK_TORRANCE:{
      struct Vector2f m = warp_square_to_unit_gaussian(sample);
      float alpha = bsdf->cook_torrance.alpha;
      m.x *= (alpha / sqrtf(2.0f));
      m.y *= (alpha / sqrtf(2.0f));
      struct Vector3f wh = normalize_v3f((struct Vector3f){m.x,m.y,1.0f});
      rec->wo = sub_v3f(smul_v3f(2.0f * dot_v3f(rec->wi,wh),wh),rec->wi);
      float pdf = pdf_bsdf(bsdf,rec);
      color = smul_v3f(rec->wo.z/pdf,eval_bsdf(bsdf,rec));
      if(rec->wo.z < 0 || isnan(color.x) || isnan(color.y) || isnan(color.z)){
        color = (struct Vector3f){0.0f,0.0f,0.0f};
      }
    break;
    }
    case BSDF_ZERO:{
      printf("Error: Attempting to sample bsdf of type BSDF_ZERO(constant and equal to zero)!\n");
      exit(1);
    break;
    }
  }
  return color;
}

