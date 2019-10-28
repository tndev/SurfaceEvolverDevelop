#ifndef AABBTREE_H_
#define AABBTREE_H_

#include <memory>
#include "PrimitiveBox.h"

#define uint unsigned int
#define uint32 uint32_t

using Tri = StructGeom::Triangle;

class AABBTree
{
public:
	Box3 bbox;

	AABBTree* left = nullptr;
	AABBTree* right = nullptr;

	// split axis
	// x = 0, y = 1, z = 2
	uint axis = 2;
	float splitPosition = 0.0f;
	uint depth = 0;

	std::vector<Tri> triangles = {};

	AABBTree();
	AABBTree(std::vector<Tri>& triangles, Box3 bbox, uint depthLeft);
	~AABBTree();

	bool isLeaf();

	void construct(std::vector<Tri>& triangles, uint depthLeft);
	std::vector<Geometry> getAABBGeomsOfDepth(uint depth); // for visualisation
	std::vector<Geometry> getAABBLeafGeoms(); // for visualisation
private:
	const uint MAX_DEPTH = 100;

	float getSplitPosition(std::vector<Tri>& triangles, std::vector<Tri>* out_left, std::vector<Tri>* out_right);
	float getCostEstimate(float splitPos, uint nLeft, uint nRight);
	bool hasEnoughBranching(size_t nLeftTris, size_t nRightTris, size_t nTris);
	void filterTriangles();
};

#endif