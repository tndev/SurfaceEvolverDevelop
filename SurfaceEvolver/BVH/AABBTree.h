#ifndef AABBTREE_H_
#define AABBTREE_H_

#include <stack>
#include "../GeometryObject/PrimitiveBox.h"

#define MAX_DEPTH 20

#define uint unsigned int
#define uint32 uint32_t

using Tri = StructGeom::Triangle;

class AABBTree
{
public:
	struct AABBNode {
		AABBNode* parent = nullptr;
		AABBNode* left = nullptr;
		AABBNode* right = nullptr;

		AABBTree* tree = nullptr;

		Box3 bbox;
		uint axis = 2;
		uint depth = 0;
		float splitPosition = 0.0f;

		std::vector<uint> triangles = {};

		AABBNode();
		AABBNode(std::vector<uint>* triangles, Box3& bbox, AABBTree* tree, uint depthLeft = MAX_DEPTH, AABBNode* parent = nullptr);
		~AABBNode();

		void construct(std::vector<uint>* triangles, uint depthLeft);
		bool isALeaf();
		bool isALeafWithTriangles();

		float getSplitPosition(std::vector<uint>& triangles, std::vector<uint>* out_left, std::vector<uint>* out_right);
		float getCostEstimate(float splitPos, uint nLeft, uint nRight);
		bool hasEnoughBranching(size_t nLeftTris, size_t nRightTris, size_t nTris);
		void filterTriangles();
	};

	struct AABBRay {
		Vector3 start;
		Vector3 direction;
		Vector3 invDirection;

		float minParam = 0.0f;
		float maxParam = 0.0f;

		AABBNode* hitNode = nullptr;
		uint hitFace;

		AABBRay(Vector3& rayStart, Vector3& rayDirection, float minParam, float maxParam);
		~AABBRay();
	};

	Box3 bbox;

	uint depth = 0;

	std::vector<Tri> triangles = {};
	AABBNode* root = nullptr;
	Geometry* geom = nullptr;

	AABBTree();
	AABBTree(Geometry* geom);
	~AABBTree();

	bool hasTriangles();
	bool boxIntersectsATriangle(Box3* box);
	float boxIntersectsATriangleAtDistance(Box3* box);

	std::vector<AABBNode> flatten();
	std::vector<AABBNode> flattenToDepth(uint depth);
	std::vector<Tri> getTrianglesInABox(Box3* box);

	std::pair<bool, std::vector<float>> getRayBoxIntersection(AABBRay* ray, Box3* bbox);
	AABBRay rayIntersect(Vector3& rayStart, Vector3& rayDirection, float minParam, float maxParam);

	std::vector<Geometry> getAABBGeomsOfDepth(uint depth); // for visualisation
	std::vector<Geometry> getAABBLeafGeoms(); // for visualisation
	std::vector<Geometry> getAABBTrianglesOfDepth(uint depth); // for visualisation
};

uint depth(AABBTree* root);

#endif