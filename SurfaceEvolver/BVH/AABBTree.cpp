#include "AABBTree.h"

AABBTree::AABBTree()
{
}

AABBTree::AABBTree(const AABBTree& other)
{
	triangles = std::vector<Tri>(other.triangles);
	bbox = other.bbox;
	geom = other.geom;
	depth = other.depth;
	root = other.root;
}

AABBTree::AABBTree(Geometry* geom)
{
	this->triangles = geom->getTriangles();
	this->bbox = geom->getBoundingBox();
	this->bbox.min.addScalar(-0.001);
	this->bbox.max.addScalar(0.001);

	this->geom = geom;

	// generate index array to triangles
	std::vector<uint> triIds;
	triIds.reserve(this->triangles.size());
	uint n(0); std::generate_n(std::back_inserter(triIds), this->triangles.size(), [n]() mutable { return n++;  });

	// generate root and all that follow
	this->root = new AABBNode(&triIds, this->bbox, this);
}

AABBTree::~AABBTree()
{
}

bool AABBTree::hasTriangles()
{
	return this->triangles.size() > 0;
}

bool AABBTree::boxIntersectsATriangle(Box3* box)
{
	std::stack<AABBNode*> stack = {};
	stack.push(this->root);

	while (stack.size()) {
		AABBNode* item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
			Vector3 center = box->getCenter();
			Vector3 halfSize = box->getSize();
			halfSize = 0.5 * halfSize;
			for (auto&& t:item->triangles) {
				if (getTriangleBoundingBoxIntersection(&this->triangles.at(t), center, halfSize, 0.0f)) {
					return true;
				}
			}
		}
		else {
			if (item->left && box->intersectsBox(item->left->bbox)) {
				stack.push(item->left);
			}
			if (item->right && box->intersectsBox(item->right->bbox)) {
				stack.push(item->right);
			}
		}
	}

	return false;
}

float AABBTree::boxIntersectsATriangleAtDistance(Box3* box)
{
	std::stack<AABBNode*> stack = {};
	stack.push(this->root);

	while (stack.size()) {
		AABBNode* item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
			Vector3 center = box->getCenter();
			Vector3 halfSize = 0.5 * box->getSize();

			for (auto&& t : item->triangles) {
				if (getTriangleBoundingBoxIntersection(&this->triangles.at(t), center, halfSize, 0.0f)) {
					return getDistanceToATriangleSq2(&this->triangles.at(t), center);
				}
			}
		}
		else {
			if (item->left && box->intersectsBox(item->left->bbox)) {
				stack.push(item->left);
			}
			if (item->right && box->intersectsBox(item->right->bbox)) {
				stack.push(item->right);
			}
		}
	}

	return INFINITY;
}

std::vector<AABBTree::AABBNode> AABBTree::flatten()
{
	std::vector<AABBNode> resultArray = {};

	std::stack<AABBNode*> nodeStack = {};
	nodeStack.push(this->root);

	while (nodeStack.size()) {
		AABBNode* item = nodeStack.top();
		nodeStack.pop();

		if (item->isALeaf()) {
			resultArray.push_back(*item);
		}
		else {
			if (item->left) {
				nodeStack.push(item->left);
			}
			if (item->right) {
				nodeStack.push(item->right);
			}
		}
	}

	return resultArray;
}

std::vector<AABBTree::AABBNode> AABBTree::flattenToDepth(uint depth)
{
	std::vector<AABBNode> resultArray = {};

	std::stack<AABBNode*> nodeStack;
	nodeStack.push(this->root);

	while (nodeStack.size()) {
		AABBNode* item = nodeStack.top();
		nodeStack.pop();

		if (item->isALeaf() && item->depth == depth) {// is leaf?
			resultArray.push_back(*item);
		}
		else {
			if (item->left) {
				nodeStack.push(item->left);
			}
			if (item->right) {
				nodeStack.push(item->right);
			}
		}
	}

	return resultArray;
}

// returns all triangles that are intersecting overlapping AABB leaves
void AABBTree::getTrianglesInABox(Box3* box, std::vector<uint>* triIdBuffer)
{
	if (!triIdBuffer->empty()) triIdBuffer->clear();

	std::stack<AABBNode*> stack = {};
	stack.push(this->root);

	while (stack.size()) {
		AABBNode* item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
			for (uint i = 0; i < item->triangles.size(); i++) {
				triIdBuffer->push_back(item->triangles[i]);
			}
		} else {
			if (item->left && box->intersectsBox(item->left->bbox)) {
				stack.push(item->left);
			}
			if (item->right && box->intersectsBox(item->right->bbox)) {
				stack.push(item->right);
			}
		}
	}
}

AABBTree::AABBNode* AABBTree::getClosestNode(Vector3& point)
{
	float t_left, t_right;
	Vector3 leftCenter, rightCenter, leftSize, rightSize;
	Vector3 invLeftCenter = Vector3(), invRightCenter = Vector3();

	std::stack<AABBNode*> stack = {};
	stack.push(this->root);

	while (stack.size()) {
		AABBNode* item = stack.top();
		stack.pop();

		if (item->left || item->right) {
			leftCenter = item->left->bbox.getCenter() - point;
			invLeftCenter.set(1.0f / leftCenter.x, 1.0f / leftCenter.y, 1.0f / leftCenter.z);
			leftSize = 0.5 * item->left->bbox.getSize();

			rightCenter = item->right->bbox.getCenter() - point;
			invLeftCenter.set(1.0f / rightCenter.x, 1.0f / rightCenter.y, 1.0f / rightCenter.z);
			rightSize = 0.5 * item->right->bbox.getSize();

			t_left = fminf(dot((item->left->bbox.min - point), invLeftCenter), dot((item->left->bbox.max - point), invLeftCenter));
			t_right = fminf(dot((item->right->bbox.min - point), invRightCenter), dot((item->right->bbox.max - point), invRightCenter));

			bool leftIsNear = leftCenter.lengthSq() < rightCenter.lengthSq() && t_left < t_right;

			if (leftIsNear) {
				stack.push(item->left);
			}
			else {
				stack.push(item->right);
			}
		}
		else {
			return item;
		}
	}

	return nullptr;
}

int AABBTree::getClosestTriangleId(Vector3& point)
{
	AABBNode* closestNode = this->getClosestNode(point);
	if (!closestNode) {
		return -1;
	}
	return closestNode->triangles[0];
}

std::vector<Geometry> AABBTree::getAABBGeomsOfDepth(uint depth)
{
	std::vector<AABBNode> leafNodes = this->flattenToDepth(depth);
	std::vector<Geometry> boxGeoms = {};

	for (auto&& leaf : leafNodes) {
		float dimX = leaf.bbox.max.x - leaf.bbox.min.x;
		float dimY = leaf.bbox.max.y - leaf.bbox.min.y;
		float dimZ = leaf.bbox.max.z - leaf.bbox.min.z;
		PrimitiveBox box = PrimitiveBox(dimX, dimY, dimZ, 1, 1, 1);
		Vector3 t = leaf.bbox.min;
		box.applyMatrix(Matrix4().makeTranslation(t.x, t.y, t.z));
		boxGeoms.push_back(box);
	}

	return boxGeoms;
}

std::vector<Geometry> AABBTree::getAABBLeafGeoms()
{
	std::vector<AABBNode> leafNodes = this->flatten();
	std::vector<Geometry> boxGeoms = {};

	for (auto&& leaf : leafNodes) {
		float dimX = leaf.bbox.max.x - leaf.bbox.min.x;
		float dimY = leaf.bbox.max.y - leaf.bbox.min.y;
		float dimZ = leaf.bbox.max.z - leaf.bbox.min.z;
		PrimitiveBox box = PrimitiveBox(dimX, dimY, dimZ, 1, 1, 1);
		Vector3 t = leaf.bbox.min;
		box.applyMatrix(Matrix4().makeTranslation(t.x, t.y, t.z));
		boxGeoms.push_back(box);
	}

	return boxGeoms;
}

std::vector<Geometry> AABBTree::getAABBTrianglesOfDepth(uint depth)
{
	std::vector<AABBNode> leafNodes = this->flattenToDepth(depth);
	std::vector<Geometry> triGeoms = {};

	for (auto&& leaf : leafNodes) {
		std::vector<Geometry> leafTriGeoms = {};

		for (uint i = 0; i < leaf.triangles.size(); i++) {
			Geometry triGeom = Geometry();
			triGeom.uniqueVertices = std::vector<Vector3>(*this->triangles[leaf.triangles[i]].begin(), *this->triangles[leaf.triangles[i]].end());

			for (uint j = 0; j < 3; j++) {
				triGeom.vertices.push_back(this->triangles[leaf.triangles[i]][j]->x);
				triGeom.vertices.push_back(this->triangles[leaf.triangles[i]][j]->y);
				triGeom.vertices.push_back(this->triangles[leaf.triangles[i]][j]->z);

				// normals will be computed when merging geoms

				triGeom.vertexIndices.push_back(j);
			}
			triGeoms.push_back(triGeom);
		}
	}

	return triGeoms;
}

uint depth(AABBTree* tree)
{
	if (tree->root == nullptr) {
		return 0;
	}

	std::list<AABBTree::AABBNode*> queue;
	queue.push_back(tree->root);

	AABBTree::AABBNode* front = nullptr;
	uint depth = 0;

	while (!queue.empty()) {
		uint size = queue.size();

		while (size--) {
			front = queue.front();
			queue.pop_front();

			if (front->left != nullptr) {
				queue.push_back(front->left);
			}
			if (front->right != nullptr) {
				queue.push_back(front->right);
			}
		}

		depth++;
	}

	return depth;
}

AABBTree::AABBNode::AABBNode()
{
}

AABBTree::AABBNode::AABBNode(const AABBNode& other)
{
	parent = other.parent;
	left = other.left;
	right = other.right;
	tree = other.tree;
	bbox = other.bbox;
	axis = other.axis;
	depth = other.depth;
	splitPosition = other.splitPosition;
	triangles = std::vector<uint>(other.triangles);
}

AABBTree::AABBNode::AABBNode(std::vector<uint>* triangles, Box3& bbox, AABBTree* tree, uint depthLeft, AABBNode* parent)
{
	this->tree = tree;
	this->parent = parent;
	this->bbox = bbox;
	this->depth = MAX_DEPTH - depthLeft;

	this->construct(triangles, depthLeft);
}

AABBTree::AABBNode::~AABBNode()
{
}

bool AABBTree::AABBNode::isALeaf()
{
	return (!this->left && !this->right);
}

bool AABBTree::AABBNode::isALeafWithTriangles()
{
	return this->isALeaf() && !this->triangles.empty();
}

void AABBTree::AABBNode::construct(std::vector<uint>* triangles, uint depthLeft)
{
	if (depthLeft == 0 || triangles->size() <= 2) {
		this->triangles = *triangles;
		this->filterTriangles();
		return;
	}

	// choose axis with maximum extent to split by. this appears to be a bit better than to cycle axes
	Vector3 bboxSize = this->bbox.getSize();
	if (bboxSize.x > bboxSize.y&& bboxSize.x > bboxSize.z) {
		this->axis = 0;
	}
	else if (bboxSize.y > bboxSize.z) {
		this->axis = 1;
	}
	else {
		this->axis = 2;
	}

	std::vector<uint> leftTriangles = {};
	std::vector<uint> rightTriangles = {};

	this->splitPosition = this->getSplitPosition(*triangles, &leftTriangles, &rightTriangles); // classify triangles

	const bool stopBranching = !this->hasEnoughBranching(leftTriangles.size(), rightTriangles.size(), triangles->size());
	if (stopBranching) {
		this->triangles = *triangles;
		this->filterTriangles();
	}
	else {
		if (leftTriangles.size() > 0) {
			Box3 bboxLeft = this->bbox;
			bboxLeft.max.setCoordById(this->splitPosition, this->axis);
			this->left = new AABBNode(&leftTriangles, bboxLeft, this->tree, depthLeft - 1, this);

			if (this->left->triangles.empty() && this->left->left == nullptr && this->left->right == nullptr) {
				this->left = nullptr;
			}
		}
		if (rightTriangles.size() > 0) {
			Box3 bboxRight = this->bbox;
			bboxRight.min.setCoordById(this->splitPosition, this->axis);
			this->right = new AABBNode(&rightTriangles, bboxRight, this->tree, depthLeft - 1, this);

			if (this->right->triangles.empty() && this->right->left == nullptr && this->right->right == nullptr) {
				this->right = nullptr;
			}
		}
	}
}

float AABBTree::AABBNode::getSplitPosition(std::vector<uint>& triangles, std::vector<uint>* out_left, std::vector<uint>* out_right)
{
	std::vector<float> splitPositions = {};
	std::vector<float> splitsLeft = {};
	std::vector<float> splitsRight = {};

	const uint CUTS = 8;
	for (uint i = 1; i <= CUTS; i++) {
		splitPositions.push_back(bbox.min.getCoordById(axis) * (1.0f - ((float)i / (CUTS + 1.0f))) + bbox.max.getCoordById(axis) * ((float)i / (CUTS + 1.0f)));
		splitsLeft.push_back(0.0f);
		splitsRight.push_back(0.0f);
	}

	// count triangles for each split, we don't want to fill arrays here (lot of wasted cycles and memory ops)
	for (uint i = 0; i < triangles.size(); i++) {
		float min = std::fminf(
			this->tree->triangles[triangles[i]][0]->getCoordById(axis), 
			std::fminf(
				this->tree->triangles[triangles[i]][1]->getCoordById(axis),
				this->tree->triangles[triangles[i]][2]->getCoordById(axis)
			)
		);		
		float max = std::fmaxf(
			this->tree->triangles[triangles[i]][0]->getCoordById(axis),
			std::fmaxf(
				this->tree->triangles[triangles[i]][1]->getCoordById(axis),
				this->tree->triangles[triangles[i]][2]->getCoordById(axis)
			)
		);

		for (uint j = 0; j < splitPositions.size(); j++) {
			if (min <= splitPositions[j]) {
				++splitsLeft[j];
			}
			if (max >= splitPositions[j]) {
				++splitsRight[j];
			}
		}
	}

	// pick the best split location
	float bestCost = FLT_MAX;
	float bestSplitPosition = bbox.min.getCoordById(axis) * 0.5f + bbox.max.getCoordById(axis) * 0.5f;
	for (uint i = 0; i < splitPositions.size(); i++) {
		float cost = getCostEstimate(splitPositions[i], splitsLeft[i], splitsRight[i]);
		if (cost < bestCost) {
			bestCost = cost;
			bestSplitPosition = splitPositions[i];
		}
	}

	// fill left and right arrays now that best split position is known
	for (uint i = 0; i < triangles.size(); i++) {
		float min = std::fminf(
			this->tree->triangles[triangles[i]][0]->getCoordById(axis),
			std::fminf(
				this->tree->triangles[triangles[i]][1]->getCoordById(axis),
				this->tree->triangles[triangles[i]][2]->getCoordById(axis)
			)
		);
		float max = std::fmaxf(
			this->tree->triangles[triangles[i]][0]->getCoordById(axis),
			std::fmaxf(
				this->tree->triangles[triangles[i]][1]->getCoordById(axis),
				this->tree->triangles[triangles[i]][2]->getCoordById(axis)
			)
		);

		if (min <= bestSplitPosition) {
			out_left->push_back(triangles[i]);
		}
		if (max >= bestSplitPosition) {
			out_right->push_back(triangles[i]);
		}
	}

	return bestSplitPosition;
}

float AABBTree::AABBNode::getCostEstimate(float splitPos, uint nLeft, uint nRight)
{
	Box3 l_bbox = this->bbox;
	l_bbox.max.setCoordById(splitPos, this->axis);
	Box3 r_bbox = this->bbox;
	r_bbox.min.setCoordById(splitPos, this->axis);
	Vector3 ls = l_bbox.getSize() / 1000.0f;
	Vector3 rs = r_bbox.getSize() / 1000.0f;
	float leftArea = ls.x * ls.y * 2.0f + ls.y * ls.z * 2.0f + ls.z * ls.x * 2.0f;
	float rightArea = rs.x * rs.y * 2.0f + rs.y * rs.z * 2.0f + rs.z * rs.x * 2.0f;

	return leftArea * nLeft + rightArea * nRight;
}

bool AABBTree::AABBNode::hasEnoughBranching(size_t nLeftTris, size_t nRightTris, size_t nTris)
{
	return nLeftTris + nRightTris < 1.5f * nTris;
}

void AABBTree::AABBNode::filterTriangles()
{
	Vector3 bboxCenter = bbox.getCenter();
	Vector3 bboxHalfSize = 0.5f * bbox.getSize();
	std::vector<uint> newTriangles = {};
	for (uint i = 0; i < this->triangles.size(); i++) {
		if (getTriangleBoundingBoxIntersection(&this->tree->triangles.at(this->triangles[i]), bboxCenter, bboxHalfSize)) {
			newTriangles.push_back(this->triangles[i]);
		}
	}
	this->triangles = newTriangles;
}