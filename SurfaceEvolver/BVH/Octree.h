#ifndef OCTREE_H_
#define OCTREE_H_

#include <iostream>
#include <stack>
#include <chrono>
#include "AABBTree.h"
#include "../ExportImport/VTKExporter.h"
#include "../GeometryObject/PrimitiveBox.h"
#include "../SDF/Grid.h"

#define MAX_OCTREE_DEPTH 100

#define uint unsigned int

using Tri = StructGeom::Triangle;

class Octree
{
public:
	struct OctreeNode {
		OctreeNode** children = nullptr;
		// std::vector<OctreeNode*> children = {};

		OctreeNode* parent = nullptr;
		Octree* tree = nullptr;

		Box3 box;

		float centroidDistance = INFINITY; // infinity unless it's a leaf
		uint depthLeft = MAX_OCTREE_DEPTH;

		OctreeNode();
		OctreeNode(const OctreeNode& other);
		OctreeNode(Octree* tree, Box3 box, OctreeNode* parent = nullptr, uint depthLeft = MAX_OCTREE_DEPTH);
		bool intersectsPrimitives(Box3* box);
		bool isLargerThanLeaf(float* size);
		bool isALeaf();

		void getLeafNodes(std::vector<OctreeNode*>* leafBuffer);
		void getLeafBoxes(std::vector<Box3*>* boxBuffer);
		void getLeafBoxesAndValues(std::vector<Box3*>* boxBuffer, std::vector<float>* valueBuffer);

		void applyMatrix(Matrix4& m);
	};

	OctreeNode* root = nullptr;
	AABBTree* aabbTree = nullptr;
	Box3 cubeBox;
	Box3 bbox;

	uint depth = 0;
	float leafSize = 1.0f;

	float leaf_retrieve_time;

	uint nodeCount = 0;

	Octree();
	// expecting a constructed AABBTree for fast lookup
	Octree(const Octree& other);
	Octree(AABBTree* aabbTree, Box3 bbox, uint resolution);
	~Octree();

	void getAllNodes(std::vector<OctreeNode>* nodeBuffer);

	void getLeafBoxGeoms(std::vector<Geometry>* geoms); // for visualisation
	void GenerateFullOctreeBoxVisualisation(VTKExporter& e);
	void GenerateLeafCellVisualisation(VTKExporter& e, bool visualizeCentroids = true);

	void setLeafValueToScalarGrid(Grid* grid);
	void setConstantValueToScalarGrid(Grid* grid, float value);

	void applyMatrix(Matrix4& m);
};

#endif