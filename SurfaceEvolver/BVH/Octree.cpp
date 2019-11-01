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
		std::vector<Box3> boxes = getOctantBoxes(&size);
		for (uint i = 0; i < boxes.size(); i++) {
			this->children.push_back(new OctreeNode(this->tree, boxes[i], this, depthLeft - 1));
		}
	}
}

bool Octree::OctreeNode::intersectsTriangles(Box3* box)
{
	Vector3 center = box->getCenter();
	Vector3 halfSize = box->getSize();
	halfSize = 0.5 * halfSize;

	// this also filters out boxes that do not intersect with the root AABB
	/* std::vector<Tri> triangles = this->tree->aabbTree->getTrianglesInABox(box);
	for (auto&& t : triangles) {
		if (getTriangleBoundingBoxIntersection(t, center, halfSize, 0.0f)) {
			return true;
		}
	} */
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

std::vector<Box3> Octree::OctreeNode::getOctantBoxes(Vector3* size)
{
	std::vector<Box3> result = {};
	for (uint i = 0; i < 2; i++) {
		for (uint j = 0; j < 2; j++) {
			for (uint k = 0; k < 2; k++) {
				Vector3 offset0 = multiply(Vector3((float)i / 2.0f, (float)j / 2.0f, (float)k / 2.0f), *size);
				Vector3 offset1 = multiply(Vector3((float)(i + 1) / 2.0f, (float)(j + 1) / 2.0f, (float)(k + 1) / 2.0f), *size);
				Box3 box = Box3(this->box.min + offset0, this->box.min + offset1);
				if (intersectsTriangles(&box)) {
					result.push_back(box);
				}				
			}
		}
	}
	return result;
}

using Leaf = Octree::OctreeNode;
void Octree::OctreeNode::getLeafNodes(std::vector<Leaf>* leafBuffer)
{
	std::stack<Leaf*> stack = {};
	stack.push(this);

	while (stack.size()) {
		Leaf item = *stack.top();
		stack.pop();

		if (item.isALeaf()) {
			leafBuffer->push_back(item);
		} else {
			for (uint i = 0; i < item.children.size(); i++) {
				stack.push(item.children[i]);
			}
		}
	}
}

void Octree::OctreeNode::getLeafBoxes(std::vector<Box3>* boxBuffer)
{
	std::stack<Leaf*> stack = {};
	stack.push(this);

	while (stack.size()) {
		Leaf item = *stack.top();
		stack.pop();

		if (item.isALeaf()) {
			boxBuffer->push_back(item.box);
		}
		else {
			for (uint i = 0; i < item.children.size(); i++) {
				stack.push(item.children[i]);
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
	std::vector<Box3> boxBuffer = {};
	this->root->getLeafBoxes(&boxBuffer);
	// === Timed code ============
	auto endGetLeaves = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedAABBLeaves = (endGetLeaves - startGetLeaves);
	std::cout << "Octree leaf nodes retrieved after " << elapsedAABBLeaves.count() << " seconds" << std::endl;
	
	for (auto&& b : boxBuffer) {
		float dimX = b.max.x - b.min.x;
		float dimY = b.max.y - b.min.y;
		float dimZ = b.max.z - b.min.z;
		PrimitiveBox box = PrimitiveBox(dimX, dimY, dimZ, 1, 1, 1);
		Vector3 t = b.min;
		box.applyMatrix(Matrix4().makeTranslation(t.x, t.y, t.z));
		geoms->push_back(box);
	}
}
