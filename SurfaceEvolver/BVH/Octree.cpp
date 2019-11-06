#include "Octree.h"

Octree::OctreeNode::OctreeNode()
{
}

Octree::OctreeNode::OctreeNode(Octree* tree, Box3 box, OctreeNode* parent, uint depthLeft)
{
	this->tree = tree; // so it knows what tree it belongs to
	this->parent = parent;
	this->box = box;
	Vector3 size = this->box.getSize();
	this->tree->nodeCount++;

	if (this->isLargerThanLeaf(&size) && depthLeft > 0) {
		for (uint i = 0; i < 2; i++) {
			for (uint j = 0; j < 2; j++) {
				for (uint k = 0; k < 2; k++) {
					Vector3 offset0 = multiply(Vector3((float)i / 2.0f, (float)j / 2.0f, (float)k / 2.0f), size);
					Vector3 offset1 = multiply(Vector3((float)(i + 1) / 2.0f, (float)(j + 1) / 2.0f, (float)(k + 1) / 2.0f), size);
					Box3 box = Box3(this->box.min + offset0, this->box.min + offset1);
					if (intersectsTriangles(&box)) {
						this->children.push_back(new OctreeNode(this->tree, box, this, depthLeft - 1));
					}
				}
			}
		}
	}
}

bool Octree::OctreeNode::intersectsTriangles(Box3* box)
{
	return this->tree->aabbTree->boxIntersectsATriangle(box);
}

bool Octree::OctreeNode::isLargerThanLeaf(Vector3* size)
{
	return (size->x > this->tree->leafSize); // they're all cubes
}

bool Octree::OctreeNode::isALeaf()
{
	return this->children.empty() && !this->isLargerThanLeaf(&this->box.getSize());
}

using Leaf = Octree::OctreeNode;
void Octree::OctreeNode::getLeafNodes(std::vector<Leaf*>* leafBuffer)
{
	std::stack<Leaf*> stack = {};
	stack.push(this);

	while (stack.size()) {
		Leaf* item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
			leafBuffer->push_back(item);
		} else {
			for (uint i = 0; i < item->children.size(); i++) {
				stack.push(item->children[i]);
			}
		}
	}
}

void Octree::OctreeNode::getLeafBoxes(std::vector<Box3*>* boxBuffer)
{
	std::stack<Leaf*> stack = {};
	stack.push(this);

	while (stack.size()) {
		Leaf* item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
			boxBuffer->push_back(&item->box);
		}
		else {
			for (uint i = 0; i < item->children.size(); i++) {
				stack.push(item->children[i]);
			}
		}
	}
}

Octree::Octree()
{
}

Octree::Octree(AABBTree* aabbTree, Box3 bbox, uint resolution)
{
	Vector3 size = bbox.getSize();
	float maxDim = std::max({ size.x, size.y, size.z });

	// this cube box will be subdivided
	Box3 cubeBox = Box3(bbox.min, bbox.min + Vector3(maxDim, maxDim, maxDim));

	this->leafSize = maxDim / resolution;
	this->cubeBox = cubeBox;
	this->aabbTree = aabbTree; // for fast lookup

	this->root = new OctreeNode(this, this->cubeBox);
}

Octree::~Octree()
{
}

void Octree::getLeafBoxGeoms(std::vector<Geometry>* geoms)
{
	auto startGetLeaves = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	std::vector<Box3*> boxBuffer = {};
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

void Octree::setLeafValueToScalarGrid(Grid* grid, float value)
{
	auto startGetLeaves = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	std::vector<Box3*> boxBuffer = {};
	this->root->getLeafBoxes(&boxBuffer);
	// === Timed code ============
	auto endGetLeaves = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedAABBLeaves = (endGetLeaves - startGetLeaves);
	std::cout << "Octree leaf nodes retrieved after " << elapsedAABBLeaves.count() << " seconds" << std::endl;

	uint Nx = grid->Nx, Ny = grid->Ny, Nz = grid->Nz;
	float scaleX = grid->scale.x, scaleY = grid->scale.y, scaleZ = grid->scale.z;
	float gMinX = grid->bbox.min.x, gMinY = grid->bbox.min.y, gMinZ = grid->bbox.min.z;
	
	grid->max = grid->max < value ? value + 1 : grid->max;

	uint ix, iy, iz, gridPos;

	for (auto&& b : boxBuffer) {
		// transform from real space to grid index space
		ix = (uint)std::round((b->min.x - gMinX) * Nx / scaleX);
		iy = (uint)std::round((b->min.y - gMinY) * Ny / scaleY);
		iz = (uint)std::round((b->min.z - gMinZ) * Nz / scaleZ);

		gridPos = Nx * Ny * iz + Nx * iy + ix;
		grid->field[gridPos] = value;
		grid->frozenCells[gridPos] = true; // freeze cell;
	}
}
