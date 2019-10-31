
#ifndef OCTREE_H_
#define OCTREE_H_

#include <stack>
#include "AABBTree.h"
#include "PrimitiveBox.h"

#define uint unsigned int

using Tri = StructGeom::Triangle;

class Octree
{
public:
	struct OctreeNode {
		std::vector<OctreeNode> children = {};

		OctreeNode* parent = nullptr;
		Octree* tree = nullptr;

		Box3 box;

		OctreeNode(Octree* tree, Box3 box, OctreeNode* parent = nullptr);
		bool shouldSubdivide();
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

	Octree();
	// expecting a constructed AABBTree for fast lookup
	Octree(AABBTree* aabbTree, Box3 bbox, float leafSize);
	~Octree();
	std::vector<OctreeNode> getLeafNodes();

	std::vector<Geometry> getLeafBoxGeoms(); // for visualisation
};

#endif