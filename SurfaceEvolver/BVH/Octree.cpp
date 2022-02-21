#include "Octree.h"

Octree::OctreeNode::OctreeNode(const OctreeNode& other)
{
	tree = other.tree;
	box = other.box;
	centroidDistance = other.centroidDistance;
	children = other.children;
}

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

// ======================== Helper macros for Octree subdivision ======================

#define SET_BOX_MIN_COORD(b, B, i, j, k)									\
		b.min.set(															\
				B.min.x + ((double)i / 2.0) * size,							\
				B.min.y + ((double)j / 2.0) * size,							\
				B.min.z + ((double)k / 2.0) * size							\
		);

#define SET_BOX_MAX_COORD(b, B, i, j, k)									\
		b.max.set(															\
				B.min.x + ((double)(i + 1) / 2.0) * size,					\
				B.min.y + ((double)(j + 1) / 2.0) * size,					\
				B.min.z + ((double)(k + 1) / 2.0) * size					\
		);

#define GET_CUBE_SIZE(b)													\
		b.max.x - b.min.x;

Octree::OctreeNode::OctreeNode(Octree* tree, Box3 box, uint depthLeft)
{
	this->tree = tree; // so it knows what tree it belongs to
	this->box = box;
	double size = GET_CUBE_SIZE(this->box);
	this->depthLeft == depthLeft;

	if (this->isLargerThanLeaf(size) && depthLeft > 0) {
		this->children = {};
		children.reserve(8);
		Box3 childBox{};
		int i, j, k;
		for (i = 0; i < 2; i++) {
			for (j = 0; j < 2; j++) {
				for (k = 0; k < 2; k++) {
					SET_BOX_MIN_COORD(childBox, this->box, i, j, k);
					SET_BOX_MAX_COORD(childBox, this->box, i, j, k);

					if (intersectsPrimitives(childBox)) {
						this->children.emplace_back(std::make_shared<OctreeNode>(
							OctreeNode(this->tree, childBox, depthLeft - 1)));
					}
				}
			}
		}

		children.shrink_to_fit();
	}
	else {
		// complete leaf construction by computing the mesh distance
		std::vector<uint> intersectedTriangleIds = {};
		this->tree->aabbTree->getPrimitivesInABox(&this->box, &intersectedTriangleIds);
		double distSq, resultDistSq = FLT_MAX;
		Vector3 center = this->box.getCenter();
		uint id = 0;

		for (auto& ti : intersectedTriangleIds) {
			distSq = getDistanceToAPrimitiveSq(this->tree->aabbTree->primitives[ti], center);
			resultDistSq = distSq < resultDistSq ? distSq : resultDistSq;
		}
		this->centroidDistance = sqrt(resultDistSq);
	}

	this->tree->nodeCount++;
}

bool Octree::OctreeNode::intersectsPrimitives(Box3& box)
{
	return this->tree->aabbTree->boxIntersectsAPrimitive(&box);
}

bool Octree::OctreeNode::isLargerThanLeaf(double size)
{
	return (size > this->tree->leafSize); // they're all cubes
}

bool Octree::OctreeNode::isALeaf()
{
	double size = this->box.max.x - this->box.min.x;
	return (this->children.empty() && !this->isLargerThanLeaf(size));
}

using Leaf = Octree::OctreeNode;
void Octree::OctreeNode::getLeafNodes(std::vector<Leaf*>& leafBuffer)
{
	std::stack<Leaf*> stack = {};
	stack.push(this);

	while (!stack.empty()) {
		Leaf* item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
			leafBuffer.push_back(item);
		} else {
			for (auto& child : item->children)
			{
				stack.push(child.get());
			}
		}
	}
}

void Octree::OctreeNode::getLeafBoxes(std::vector<Box3*>& boxBuffer)
{
	std::stack<Leaf*> stack = {};
	stack.push(this);

	while (!stack.empty()) {
		Leaf* item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
			boxBuffer.push_back(&item->box);
		}
		else {
			for (auto& child : item->children)
			{
				stack.push(child.get());
			}
		}
	}
}

void Octree::OctreeNode::getLeafBoxesAndValues(std::vector<Box3*>& boxBuffer, std::vector<double>& valueBuffer)
{
	std::stack<Leaf*> stack = {};
	stack.push(this);

	while (!stack.empty()) {
		Leaf* item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
			boxBuffer.push_back(&item->box);
			valueBuffer.push_back(item->centroidDistance);
		}
		else {
			for (auto& child : item->children)
			{
				stack.push(child.get());
			}
		}
	}
}

void Octree::OctreeNode::applyMatrix(Matrix4& m)
{
	box.min.applyMatrix4(m);
	box.max.applyMatrix4(m);
}

Octree::Octree(const Octree& other)
{
	root = other.root;
	aabbTree = other.aabbTree;
	bbox = other.bbox;
	cubeBox = other.cubeBox;
	depth = other.depth;
	leafSize = other.leafSize;

	leaf_retrieve_time = other.leaf_retrieve_time;
	nodeCount = other.nodeCount;
}

Octree::Octree(const std::shared_ptr<AABBTree>& aabbTree, const Box3& bbox, uint resolution) {

	/*
	this->bbox = bbox;
	Vector3 boxCenter = bbox.getCenter();

	Vector3 gridCenterMin = Vector3(std::floor(boxCenter.x / leafSize - 0.5), std::floor(boxCenter.y / leafSize - 0.5), std::floor(boxCenter.z / leafSize - 0.5));
	Vector3 gridCenterMax = Vector3(std::ceil(boxCenter.x / leafSize - 0.5), std::ceil(boxCenter.y / leafSize - 0.5), std::ceil(boxCenter.z / leafSize - 0.5));
	Vector3 newGridCenter = (gridCenterMin + gridCenterMax) * leafSize / 2.0;

	Vector3 boxSize = bbox.getSize();
	double maxDim = std::max({ boxSize.x, boxSize.y, boxSize.z });
	depth = std::floor(log2(maxDim / leafSize));
	double boxHalfDim = pow(2, depth) * leafSize;

	cubeBox = Box3(newGridCenter.clone().subScalar(boxHalfDim), newGridCenter.clone().addScalar(boxHalfDim));

	this->leafSize = leafSize;
	this->aabbTree = aabbTree;

	root = new OctreeNode(this, cubeBox);*/

	this->bbox = bbox;
	Vector3 size = bbox.getSize();
	double maxDim = std::max({ size.x, size.y, size.z });

	// this cube box will be subdivided
	Box3 cubeBox = Box3(bbox.min, bbox.min + Vector3(maxDim, maxDim, maxDim));
	cubeBox.expandByFactor(1.1);
	this->bbox.expandByFactor(1.1);

	this->leafSize = maxDim / resolution;
	this->cubeBox = cubeBox;
	this->aabbTree = aabbTree; // for fast lookup

	this->root = std::make_shared<OctreeNode>(OctreeNode(this, this->cubeBox));
}

void Octree::getAllNodes(std::vector<OctreeNode>& nodeBuffer)
{
	std::stack<OctreeNode*> stack = {};
	stack.push(this->root.get());

	while (!stack.empty()) {
		OctreeNode* item = stack.top();
		stack.pop();

		nodeBuffer.push_back(*item);

		for (auto& child : item->children)
		{
			stack.push(child.get());
		}
		delete item;
	}
}

void Octree::getLeafBoxGeoms(std::vector<Geometry>& geoms)
{
	auto startGetLeaves = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	std::vector<Box3*> boxBuffer = {};
	std::vector<double> valueBuffer = {};
	this->root->getLeafBoxes(boxBuffer);
	// === Timed code ============
	auto endGetLeaves = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsedAABBLeaves = (endGetLeaves - startGetLeaves);
	std::cout << "Octree leaf nodes retrieved after " << elapsedAABBLeaves.count() << " seconds" << std::endl;

	geoms.reserve(boxBuffer.size());
	for (auto&& b : boxBuffer) {
		double dimX = b->max.x - b->min.x;
		double dimY = b->max.y - b->min.y;
		double dimZ = b->max.z - b->min.z;
		PrimitiveBox box = PrimitiveBox(dimX, dimY, dimZ, 1, 1, 1);
		Vector3 t = b->min;
		box.applyMatrix(Matrix4().makeTranslation(t.x, t.y, t.z));
		geoms.push_back(box);
	}
}

void Octree::GenerateFullOctreeBoxVisualisation(VTKExporter& e)
{
	std::cout << "--------------------------------------------" << std::endl;
	std::vector<Geometry> boxGeoms = {};
	std::vector<OctreeNode> nodeBuffer = {};
	std::cout << "obtaining Octree nodes..." << std::endl;
	this->getAllNodes(nodeBuffer);
	std::cout << nodeBuffer.size() << "Octree nodes retrieved" << std::endl;

	std::cout << "generating box geometries..." << std::endl;
	for (auto&& n : nodeBuffer) {
		double dimX = n.box.max.x - n.box.min.x;
		double dimY = n.box.max.y - n.box.min.y;
		double dimZ = n.box.max.z - n.box.min.z;
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
	this->getLeafBoxGeoms(boxGeoms);
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

void Octree::setLeafValueToScalarGrid(Grid& grid)
{
	auto startGetLeaves = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	std::vector<Box3*> boxBuffer = {};
	std::vector<double> valueBuffer = {};
	this->root->getLeafBoxesAndValues(boxBuffer, valueBuffer);
	// === Timed code ============
	auto endGetLeaves = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsedAABBLeaves = (endGetLeaves - startGetLeaves);
	this->leaf_retrieve_time = elapsedAABBLeaves.count();

	size_t NLeaves = boxBuffer.size();
	uint Nx = grid.Nx, Ny = grid.Ny, Nz = grid.Nz;
	double scaleX = grid.scale.x, scaleY = grid.scale.y, scaleZ = grid.scale.z;
	double gMinX = grid.cubeBox.min.x, gMinY = grid.cubeBox.min.y, gMinZ = grid.cubeBox.min.z;
	
	uint ix, iy, iz, gridPos;

	for (uint i = 0; i < NLeaves; i++) {
		// transform from real space to grid index space
		ix = (uint)std::floor((0.5 * (boxBuffer[i]->min.x + boxBuffer[i]->max.x) - gMinX) * Nx / scaleX);
		iy = (uint)std::floor((0.5 * (boxBuffer[i]->min.y + boxBuffer[i]->max.y) - gMinY) * Ny / scaleY);
		iz = (uint)std::floor((0.5 * (boxBuffer[i]->min.z + boxBuffer[i]->max.z) - gMinZ) * Nz / scaleZ);

		gridPos = Nx * Ny * iz + Nx * iy + ix;
		grid.field[gridPos] = valueBuffer[i];
		grid.frozenCells[gridPos] = true; // freeze initial condition
	}

	Box3 targetBox = bbox;
	grid.clip(targetBox);
}

void Octree::setConstantValueToScalarGrid(Grid& grid, double value)
{
	auto startGetLeaves = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	std::vector<Box3*> boxBuffer = {};
	this->root->getLeafBoxes(boxBuffer);
	// === Timed code ============
	auto endGetLeaves = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsedAABBLeaves = (endGetLeaves - startGetLeaves);
	this->leaf_retrieve_time = elapsedAABBLeaves.count();

	uint Nx = grid.Nx, Ny = grid.Ny, Nz = grid.Nz;
	double scaleX = grid.scale.x, scaleY = grid.scale.y, scaleZ = grid.scale.z;
	double gMinX = grid.cubeBox.min.x, gMinY = grid.cubeBox.min.y, gMinZ = grid.cubeBox.min.z;

	grid.max = grid.max < value ? value + 1 : grid.max;

	uint ix, iy, iz, gridPos;

	for (auto&& b : boxBuffer) {
		// transform from real space to grid index space
		ix = (uint)std::round((b->min.x - gMinX) * Nx / scaleX);
		iy = (uint)std::round((b->min.y - gMinY) * Ny / scaleY);
		iz = (uint)std::round((b->min.z - gMinZ) * Nz / scaleZ);

		gridPos = Nx * Ny * iz + Nx * iy + ix;
		grid.field[gridPos] = value;
		grid.frozenCells[gridPos] = true; // freeze initial condition
	}

	Box3 targetBox = bbox;
	grid.clip(targetBox);
}

void Octree::applyMatrix(Matrix4& m)
{
	/*cubeBox.min.applyMatrix4(m);
	cubeBox.max.applyMatrix4(m);

	std::stack<OctreeNode*> stack = {};
	stack.push(this->root.get());
	uint i; bool largerThanLeaf;
	double size;

	while (!stack.empty()) {
		OctreeNode* item = stack.top();
		stack.pop();

		item->applyMatrix(m);
		size = GET_CUBE_SIZE(item->box);

		largerThanLeaf = item->isLargerThanLeaf(size);
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
		}
		else {
			for (auto& child : item->children)
			{
				stack.push(child.get());
			}
		}
	}*/
}
