#ifndef KDTREE_H_GUARD
#define KDTREE_H_GUARD
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include "vectors.h" 

struct Photon{
  struct Vector3f point;
  struct Vector3f flux;
  struct Vector3f wi_global;
};

struct KDTreeNode{
  struct KDTreeNode *left,*right;
  struct Photon photon;
  uint32_t axis;
};

struct AABB{
  float xmin,ymin,zmin,xmax,ymax,zmax;
};

struct KDTreeNode* build_kdtree_recursive(struct Photon* photons,uint32_t num_photons,uint32_t depth);
void nn_search_kdtree(struct KDTreeNode* node,
                      struct Vector3f point,
                      struct Photon* best_photon,
                      float* min_dist);
void knn_search_kdtree(struct KDTreeNode* node,
                      struct Vector3f point,
                      uint32_t k,
                      struct Photon* best_photons,
                      float* min_dists);
#endif
