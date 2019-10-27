#ifndef AABBTREE_H_
#define AABBTREE_H_

#include <memory>
#include "Geometry.h"

class AABBTree
{
public:
	AABBTree();
	AABBTree(std::vector<StructGeom::Triangle>& triangles, unsigned int depthLeft);
	~AABBTree();

	void buildNode(std::vector<StructGeom::Triangle>& triangles, unsigned int depthLeft);
};

struct Node
{
	bool IsLeaf(void) const
	{
		// The right leaf does not use the same memory as the userdata,
		// and will always be Null (no children)
		return right == Null;
	}

	// Fat AABB for leafs, bounding AABB for branches
	AABBTree aabb;

	union
	{
		int32_t parent;
		int32_t next; // free list
	};

	union
	{
		// Child indices
		struct
		{
			int32_t left;
			int32_t right;
		};

		// Since only leaf nodes hold userdata, we can use the
		// same memory used for left/right indices to store
		// the userdata void pointer
		void* userData;
	};

	// leaf = 0, free nodes = -1
	int32_t height;

	static const int32_t Null = -1;
};

#endif