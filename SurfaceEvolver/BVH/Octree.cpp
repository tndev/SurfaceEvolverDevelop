#include "Octree.h"

Octree::OctreeNode::OctreeNode()
{
}

Octree::OctreeNode::OctreeNode(const OctreeNode& other)
{
	parent = other.parent;
	tree = other.tree;
	box = other.box;
	centroidDistance = other.centroidDistance;
	children = new OctreeNode*[8];
	for (uint i = 0; i < 8; i++) children[i] = other.children[i];
}

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

// ======================== Helper macros for Octree subdivision ======================

#define SET_BOX_MIN_COORD(b, B, i, j, k)									\
		b.min.set(															\
				B.min.x + ((float)i / 2.0f) * size,							\
				B.min.y + ((float)j / 2.0f) * size,							\
				B.min.z + ((float)k / 2.0f) * size							\
		);

#define SET_BOX_MAX_COORD(b, B, i, j, k)									\
		b.max.set(															\
				B.min.x + ((float)(i + 1) / 2.0f) * size,					\
				B.min.y + ((float)(j + 1) / 2.0f) * size,					\
				B.min.z + ((float)(k + 1) / 2.0f) * size					\
		);

#define GET_CUBE_SIZE(b)													\
		b.max.x - b.min.x;

Octree::OctreeNode::OctreeNode(Octree* tree, Box3 box, OctreeNode* parent, uint depthLeft)
{
	this->tree = tree; // so it knows what tree it belongs to
	this->parent = parent;
	this->box = box;
	float size = GET_CUBE_SIZE(this->box);
	this->depthLeft == depthLeft;

	if (this->isLargerThanLeaf(&size) && depthLeft > 0) {
		this->children = new OctreeNode * [8];
		uint chId = 0;
		Box3 childBox = Box3();
		int i, j, k;
		for (i = 0; i < 2; i++) {
			for (j = 0; j < 2; j++) {
				for (k = 0; k < 2; k++) {
					SET_BOX_MIN_COORD(childBox, this->box, i, j, k);
					SET_BOX_MAX_COORD(childBox, this->box, i, j, k);

					if (intersectsPrimitives(&childBox)) {
						this->children[chId] = new OctreeNode(this->tree, childBox, this, depthLeft - 1);
					}
					else {
						this->children[chId] = nullptr;
					}
					chId++;
				}
			}
		}
	}
	else {
		// complete leaf construction by computing the mesh distance
		std::vector<uint> intersectedTriangleIds = {};
		this->tree->aabbTree->getPrimitivesInABox(&this->box, &intersectedTriangleIds);
		float distSq, resultDistSq = FLT_MAX;
		Vector3 center = this->box.getCenter();
		uint id = 0;

		for (auto&& ti : intersectedTriangleIds) {
			distSq = getDistanceToAPrimitiveSq(this->tree->aabbTree->primitives[ti], center);
			resultDistSq = distSq < resultDistSq ? distSq : resultDistSq;
		}
		this->centroidDistance = sqrt(resultDistSq);
	}

	this->tree->nodeCount++;
}

bool Octree::OctreeNode::intersectsPrimitives(Box3* box)
{
	return this->tree->aabbTree->boxIntersectsAPrimitive(box);
}

bool Octree::OctreeNode::isLargerThanLeaf(float* size)
{
	return (*size > this->tree->leafSize); // they're all cubes
}

bool Octree::OctreeNode::isALeaf()
{
	float size = this->box.max.x - this->box.min.x;
	return (this->children == nullptr && !this->isLargerThanLeaf(&size));
}

using Leaf = Octree::OctreeNode;
void Octree::OctreeNode::getLeafNodes(std::vector<Leaf*>* leafBuffer)
{
	std::stack<Leaf*> stack = {};
	stack.push(this);
	uint i;

	while (stack.size()) {
		Leaf* item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
			leafBuffer->push_back(item);
		} else {
			for (i = 0; i < 8; i++) {
				if (item->children[i]) stack.push(item->children[i]);
			}
			/* for (i = 0; i < item->children.size(); i++) {
				stack.push(item->children[i]);
			}*/
		}
	}
}

void Octree::OctreeNode::getLeafBoxes(std::vector<Box3*>* boxBuffer)
{
	std::stack<Leaf*> stack = {};
	stack.push(this);
	uint i;

	while (stack.size()) {
		Leaf* item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
			boxBuffer->push_back(&item->box);
		}
		else {
			for (i = 0; i < 8; i++) {
				if (item->children[i]) stack.push(item->children[i]);
			}
		}
	}
}

void Octree::OctreeNode::getLeafBoxesAndValues(std::vector<Box3*>* boxBuffer, std::vector<float>* valueBuffer)
{
	std::stack<Leaf*> stack = {};
	stack.push(this);
	uint i;

	while (stack.size()) {
		Leaf* item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
			boxBuffer->push_back(&item->box);
			valueBuffer->push_back(item->centroidDistance);
		}
		else {
			for (i = 0; i < 8; i++) {
				if (item->children[i]) stack.push(item->children[i]);
			}
		}
	}
}

void Octree::OctreeNode::applyMatrix(Matrix4& m)
{
	box.min.applyMatrix4(m);
	box.max.applyMatrix4(m);
}

Octree::Octree()
{
}

Octree::Octree(const Octree& other)
{
	root = other.root;
	aabbTree = other.aabbTree;
	cubeBox = other.cubeBox;
	depth = other.depth;
	leafSize = other.leafSize;

	leaf_retrieve_time = other.leaf_retrieve_time;
	nodeCount = other.nodeCount;
}

Octree::Octree(AABBTree* aabbTree, Box3 bbox, uint resolution)
{
	this->bbox = bbox;
	Vector3 size = bbox.getSize();
	float maxDim = std::max({ size.x, size.y, size.z });

	// this cube box will be subdivided
	Box3 cubeBox = Box3(bbox.min, bbox.min + Vector3(maxDim, maxDim, maxDim));
	cubeBox.expandByFactor(1.1f);
	this->bbox.expandByFactor(1.1f);

	this->leafSize = maxDim / resolution;
	this->cubeBox = cubeBox;
	this->aabbTree = aabbTree; // for fast lookup

	this->root = new OctreeNode(this, this->cubeBox);
}

Octree::~Octree()
{
	delete this->root;
}

void Octree::getAllNodes(std::vector<OctreeNode>* nodeBuffer)
{
	std::stack<OctreeNode*> stack = {};
	stack.push(this->root);
	uint i;

	while (stack.size()) {
		OctreeNode* item = stack.top();
		stack.pop();

		nodeBuffer->push_back(*item);

		for (i = 0; i < 8; i++) {
			if (item->children[i]) stack.push(item->children[i]);
		}
	}
}

void Octree::getLeafBoxGeoms(std::vector<Geometry>* geoms)
{
	auto startGetLeaves = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	std::vector<Box3*> boxBuffer = {};
	std::vector<float> valueBuffer = {};
	this->root->getLeafBoxes(&boxBuffer);
	// === Timed code ============
	auto endGetLeaves = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedAABBLeaves = (endGetLeaves - startGetLeaves);
	std::cout << "Octree leaf nodes retrieved after " << elapsedAABBLeaves.count() << " seconds" << std::endl;
	
	for (auto&& b : boxBuffer) {
		float dimX = b->max.x - b->min.x;
		float dimY = b->max.y - b->min.y;
		float dimZ = b->max.z - b->min.z;
		PrimitiveBox box = PrimitiveBox(dimX, dimY, dimZ, 1, 1, 1);
		Vector3 t = b->min;
		box.applyMatrix(Matrix4().makeTranslation(t.x, t.y, t.z));
		geoms->push_back(box);
	}
}

void Octree::GenerateFullOctreeBoxVisualisation(VTKExporter& e)
{
	std::cout << "--------------------------------------------" << std::endl;
	std::vector<Geometry> boxGeoms = {};
	std::vector<OctreeNode> nodeBuffer = {};
	std::cout << "obtaining Octree nodes..." << std::endl;
	this->getAllNodes(&nodeBuffer);
	std::cout << nodeBuffer.size() << "Octree nodes retrieved" << std::endl;

	std::cout << "generating box geometries..." << std::endl;
	for (auto&& n : nodeBuffer) {
		float dimX = n.box.max.x - n.box.min.x;
		float dimY = n.box.max.y - n.box.min.y;
		float dimZ = n.box.max.z - n.box.min.z;
		PrimitiveBox box = PrimitiveBox(dimX, dimY, dimZ, 1, 1, 1);
		Vector3 t = n.box.min;
		box.applyMatrix(Matrix4().makeTranslation(t.x, t.y, t.z));
		boxGeoms.push_back(box);
	}

	Geometry resultGeom = mergeGeometries(boxGeoms);
	e.initExport(&resultGeom, this->aabbTree->geom->name + "_Octree_allBoxes");
	std::cout << "Octree box geometries retrieved and exported" << std::endl;
}

void Octree::GenerateLeafCellVisualisation(VTKExporter& e, bool visualizeCentroids)
{
	std::cout << "--------------------------------------------" << std::endl;
	std::vector<Geometry> boxGeoms = {};
	std::cout << "obtaining Octree leaf boxes..." << std::endl;
	this->getLeafBoxGeoms(&boxGeoms);
	std::cout << boxGeoms.size() << " Octree leaf boxes retrieved" << std::endl;
	Geometry resultGeom = mergeGeometries(boxGeoms);
	e.initExport(&resultGeom, this->aabbTree->geom->name + "_Octree_leafBoxes");
	std::cout << "Octree leaf boxes exported" << std::endl;

	if (visualizeCentroids) {
		std::cout << "obtaining Octree leaf box centroids..." << std::endl;
		std::vector<Vector3> centroids = {};
		for (auto&& b : boxGeoms) {
			centroids.push_back(0.5 * (b.uniqueVertices[0] + b.uniqueVertices[6]));
		}
		e.exportPointData(centroids, this->aabbTree->geom->name + "_Octree_leafCentroids");
		std::cout << centroids.size() << " Octree leaf box centroids retrieved and exported" << std::endl;
	}
}

void Octree::setLeafValueToScalarGrid(Grid* grid)
{
	auto startGetLeaves = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	std::vector<Box3*> boxBuffer = {};
	std::vector<float> valueBuffer = {};
	this->root->getLeafBoxesAndValues(&boxBuffer, &valueBuffer);
	// === Timed code ============
	auto endGetLeaves = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedAABBLeaves = (endGetLeaves - startGetLeaves);
	this->leaf_retrieve_time = elapsedAABBLeaves.count();

	size_t NLeaves = boxBuffer.size();
	uint Nx = grid->Nx, Ny = grid->Ny, Nz = grid->Nz;
	float scaleX = grid->scale.x, scaleY = grid->scale.y, scaleZ = grid->scale.z;
	float gMinX = grid->cubeBox.min.x, gMinY = grid->cubeBox.min.y, gMinZ = grid->cubeBox.min.z;
	
	uint ix, iy, iz, gridPos;

	for (uint i = 0; i < NLeaves; i++) {
		// transform from real space to grid index space
		ix = (uint)std::floor((0.5 * (boxBuffer[i]->min.x + boxBuffer[i]->max.x) - gMinX) * Nx / scaleX);
		iy = (uint)std::floor((0.5 * (boxBuffer[i]->min.y + boxBuffer[i]->max.y) - gMinY) * Ny / scaleY);
		iz = (uint)std::floor((0.5 * (boxBuffer[i]->min.z + boxBuffer[i]->max.z) - gMinZ) * Nz / scaleZ);

		gridPos = Nx * Ny * iz + Nx * iy + ix;
		grid->field[gridPos] = valueBuffer[i];
		grid->frozenCells[gridPos] = true; // freeze initial condition
	}

	Box3 targetBox = bbox;
	grid->clip(targetBox);
}

void Octree::setConstantValueToScalarGrid(Grid* grid, float value)
{
	auto startGetLeaves = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	std::vector<Box3*> boxBuffer = {};
	this->root->getLeafBoxes(&boxBuffer);
	// === Timed code ============
	auto endGetLeaves = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedAABBLeaves = (endGetLeaves - startGetLeaves);
	this->leaf_retrieve_time = elapsedAABBLeaves.count();

	uint Nx = grid->Nx, Ny = grid->Ny, Nz = grid->Nz;
	float scaleX = grid->scale.x, scaleY = grid->scale.y, scaleZ = grid->scale.z;
	float gMinX = grid->cubeBox.min.x, gMinY = grid->cubeBox.min.y, gMinZ = grid->cubeBox.min.z;

	grid->max = grid->max < value ? value + 1 : grid->max;

	uint ix, iy, iz, gridPos;

	for (auto&& b : boxBuffer) {
		// transform from real space to grid index space
		ix = (uint)std::round((b->min.x - gMinX) * Nx / scaleX);
		iy = (uint)std::round((b->min.y - gMinY) * Ny / scaleY);
		iz = (uint)std::round((b->min.z - gMinZ) * Nz / scaleZ);

		gridPos = Nx * Ny * iz + Nx * iy + ix;
		grid->field[gridPos] = value;
		grid->frozenCells[gridPos] = true; // freeze initial condition
	}

	Box3 targetBox = bbox;
	grid->clip(targetBox);
}

void Octree::applyMatrix(Matrix4& m)
{
	cubeBox.min.applyMatrix4(m);
	cubeBox.max.applyMatrix4(m);

	std::stack<OctreeNode*> stack = {};
	stack.push(this->root);
	uint i; bool largerThanLeaf;
	float size;

	while (stack.size()) {
		OctreeNode* item = stack.top();
		stack.pop();

		item->applyMatrix(m);
		size = GET_CUBE_SIZE(item->box);

		largerThanLeaf = item->isLargerThanLeaf(&size);
		if (largerThanLeaf && item->depthLeft > 0 && item->children == nullptr) {
			item->children = new OctreeNode * [8];
			uint chId = 0; Box3 childBox = Box3();
			int i, j, k;
			for (i = 0; i < 2; i++) {
				for (j = 0; j < 2; j++) {
					for (k = 0; k < 2; k++) {
						SET_BOX_MIN_COORD(childBox, item->box, i, j, k);
						SET_BOX_MAX_COORD(childBox, item->box, i, j, k);
						if (item->intersectsPrimitives(&childBox)) {
							item->children[chId] = new OctreeNode(this, childBox, item, item->depthLeft - 1);
							stack.push(item->children[chId]);
						}
						else {
							item->children[chId] = nullptr;
						}
						chId++;
					}
				}
			}
		} else if (!largerThanLeaf && item->children != nullptr) { // turn into a leaf
			delete[] item->children;
			item->children = nullptr;
		}
		else {
			for (i = 0; i < 8; i++) {
				if (item->children[i]) stack.push(item->children[i]);
			}
		}
	}
}
