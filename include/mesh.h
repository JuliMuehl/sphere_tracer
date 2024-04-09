#ifndef MESH_H_GUARD
#define MESH_H_GUARD
#include <stddef.h>
#include "vectors.h"
struct Triangle{
  struct Vector3f p1,p2,p3;
  struct Vector3f normal;
};
struct TriangleArray{
  struct Triangle* data;
  size_t len;
};

struct TriangleArray load_stl(const char* fname);
#endif
