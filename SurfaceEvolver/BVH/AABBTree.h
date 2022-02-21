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
		std::shared_ptr<AABBNode> left = nullptr;
		std::shared_ptr<AABBNode> right = nullptr;

		AABBTree* tree = nullptr;

		Box3 bbox;
		uint axis = 2;
		uint depth = 0;
		float splitPosition = 0.0f;

		std::vector<uint> primitiveIds = {};

		AABBNode() = default;
		AABBNode(const AABBNode& other);
		AABBNode(const std::vector<uint>& primitiveIds, const Box3& bbox, AABBTree* tree, uint depthLeft);
		~AABBNode() = default;

		void construct(const std::vector<uint>& primitiveIds, uint depthLeft);
		bool isALeaf();
		bool isALeafWithPrimitives();

		double getSplitPosition(const std::vector<uint>& primitiveIds, std::vector<uint>& out_left, std::vector<uint>& out_right);
		float getAdaptivelyResampledSplitPosition(const std::vector<uint>& primitiveIds, std::vector<uint>& out_left, std::vector<uint>& out_right);
		float getCostEstimate(float splitPos, uint nLeft, uint nRight);
		bool hasEnoughBranching(size_t nLeftPrims, size_t nRightPrims, size_t nPrims);
		void filterPrimitives();  // filters only primitives which actually intersect leaf box

		void applyMatrix(Matrix4& m);

		bool useIntrinscs() const
		{
			return tree->useIntrinsics;
		}
	};

	struct AABBRay {
		Vector3 start = Vector3();
		Vector3 direction = Vector3(1.0f, 0.0f, 0.0f);
		Vector3 invDirection = Vector3(1.0f / direction.x, FLT_MAX, FLT_MAX);

		float minParam = 0.0f;
		float maxParam = FLT_MAX;
		uint hitCount = 0;
		// std::vector<uint> hitIds = {};
		// std::vector<AABBNode*> hitNodes = {};

		AABBRay();
		AABBRay(Vector3* start, Vector3* dir = nullptr);
	};

	Box3 bbox;

	uint depth = 0;
	PrimitiveType type = PrimitiveType::tri;

	std::vector<Primitive> primitives = {};
	std::shared_ptr<AABBNode> root = nullptr;
	std::shared_ptr<Geometry> geom = nullptr;

	AABBTree() = default;
	AABBTree(const AABBTree& other);
	AABBTree(const Geometry& geom, PrimitiveType type = PrimitiveType::tri);
	~AABBTree() = default;

	bool useIntrinsics = true;

	uint getDepth();

	bool hasPrimitives();
	bool boxIntersectsAPrimitive(Box3* box);
	float boxIntersectsAPrimitiveAtDistance(Box3* box);
	// returns the number of ray-mesh intersections
	uint rayIntersect(Vector3& rayOrigin, Vector3& rayDirection, float rayMinParam = 0.0f, float rayMaxParam = FLT_MAX);
	bool boolRayIntersect(Vector3& rayOrigin, Vector3& rayDirection, float rayMinParam = 0.0f, float rayMaxParam = FLT_MAX);

	std::vector<AABBNode> flatten();
	void getAllNodes(std::vector<AABBNode>* buffer);
	std::vector<AABBNode> flattenToDepth(uint depth);
	void getPrimitivesInABox(Box3* box, std::vector<uint>* primIdBuffer);

	AABBNode* getClosestNode(Vector3& point);
	int getClosestPrimitiveIdAndDist(Vector3& point, double* result);

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
bool getRayBoxIntersection(const AABBTree::AABBRay& r, Vector3* boxMin, Vector3* boxMax, float* hitParam);

#endif