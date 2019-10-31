
#ifndef OCTREE_H_
#define OCTREE_H_

#include <iostream>
#include <stack>
#include "AABBTree.h"
#include "../GeometryObject/PrimitiveBox.h"

#define MAX_OCTREE_DEPTH 10

#define uint unsigned int

using Tri = StructGeom::Triangle;

class Octree
{
public:
	struct OctreeNode {
		std::vector<OctreeNode*> children = {};

		OctreeNode* parent = nullptr;
		Octree* tree = nullptr;

		Box3 box;

		OctreeNode();
		OctreeNode(Octree* tree, Box3 box, OctreeNode* parent = nullptr, uint depthLeft = MAX_OCTREE_DEPTH);
		bool intersectsTriangles(Box3* box);
		bool isLargerThanLeaf(Vector3* size);
		bool isALeaf();
		std::vector<Box3> getOctantBoxes(Vector3* size);

		std::vector<OctreeNode> getLeafNodes();
	};

	OctreeNode* root = nullptr;
	AABBTree* aabbTree = nullptr;
	Box3 cubeBox;

	uint depth = 0;
	float leafSize = 1.0f;

	uint nodeCount = 0;

	Octree();
	// expecting a constructed AABBTree for fast lookup
	Octree(AABBTree* aabbTree, Box3 bbox, uint resolution);
	~Octree();

	std::vector<Geometry> getLeafBoxGeoms(); // for visualisation
};

#endif