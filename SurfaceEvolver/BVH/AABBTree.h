#ifndef AABBTREE_H_
#define AABBTREE_H_

#include <stack>
#include <nmmintrin.h>
#include <immintrin.h>
#include "../GeometryObject/PrimitiveBox.h"
#include "../ExportImport/VTKExporter.h"

#define MAX_DEPTH 20

#define uint unsigned int
#define uint32 uint32_t

/*
	An AABB (Axis-Aligned Bounding Box) tree is a binary space partitioning structure which accelerates
	intersection and distance queries of the source geometry primitives (vertices, edges, triangles) to 
	O(log(N)) time
*/

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

		std::vector<uint> primitiveIds = {};

		AABBNode();
		AABBNode(const AABBNode& other);
		AABBNode(std::vector<uint>* primitiveIds, Box3& bbox, AABBTree* tree, uint depthLeft = MAX_DEPTH, AABBNode* parent = nullptr);
		~AABBNode();

		void construct(std::vector<uint>* primitiveIds, uint depthLeft);
		bool isALeaf();
		bool isALeafWithPrimitives();

		float getSplitPosition(std::vector<uint>& primitiveIds, std::vector<uint>* out_left, std::vector<uint>* out_right);
		float getAdaptivelyResampledSplitPosition(std::vector<uint>& primitiveIds, std::vector<uint>* out_left, std::vector<uint>* out_right);
		float getCostEstimate(float splitPos, uint nLeft, uint nRight);
		bool hasEnoughBranching(size_t nLeftPrims, size_t nRightPrims, size_t nPrims);
		void filterPrimitives();  // returns only primitives which actually intersect leaf box

		void applyMatrix(Matrix4& m);
	};

	Box3 bbox;

	uint depth = 0;
	PrimitiveType type = PrimitiveType::tri;

	std::vector<Primitive> primitives = {};
	AABBNode* root = nullptr;
	Geometry* geom = nullptr;

	AABBTree();
	AABBTree(const AABBTree& other);
	AABBTree(Geometry* geom, PrimitiveType type = PrimitiveType::tri);
	~AABBTree();

	uint getDepth();

	bool hasPrimitives();
	bool boxIntersectsAPrimitive(Box3* box);
	float boxIntersectsAPrimitiveAtDistance(Box3* box);

	std::vector<AABBNode> flatten();
	void getAllNodes(std::vector<AABBNode>* buffer);
	std::vector<AABBNode> flattenToDepth(uint depth);
	void getPrimitivesInABox(Box3* box, std::vector<uint>* primIdBuffer);

	AABBNode* getClosestNode(Vector3& point);
	int getClosestPrimitiveIdAndDist(Vector3& point, float* result);

	void applyMatrix(Matrix4& m);

	// ---- Visualisations ----
	std::vector<Geometry> getAABBGeomsOfDepth(uint depth); // for visualisation
	std::vector<Geometry> getAABBLeafGeoms(); // for visualisation
	std::vector<Geometry> getAABBPrimitivesOfDepth(uint depth); // for visualisation

	void GenerateFullTreeBoxVisualisation(VTKExporter& e);
	void GenerateFullLeafBoxVisualisation(VTKExporter& e);
	void GenerateStepwiseLeafBoxVisualisation(VTKExporter& e);
	// ------------------------
};

uint depth(AABBTree* root);

#endif