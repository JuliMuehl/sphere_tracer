#ifndef AABBTREE_GUARD
#define AABBTREE_GUARD
#include <stdalign.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#include "geom.h"
#include "mesh.h"

typedef void AABBTree;
AABBTree* create_aabbtree(struct TriangleArray triangles);
struct Intersection intersect_ray_aabbtree(struct Ray* ray,AABBTree* aabbtree);

void set_aabbtree_frame(AABBTree* aabbtree,struct Frame frame);
void set_aabbtree_position(AABBTree* aabbtree,struct Vector3f position);
#endif
