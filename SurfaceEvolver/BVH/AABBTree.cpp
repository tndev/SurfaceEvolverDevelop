#include "AABBTree.h"

AABBTree::AABBTree()
{
}

AABBTree::AABBTree(const AABBTree& other)
{
	primitives = std::vector<Primitive>(other.primitives);
	bbox = other.bbox;
	geom = other.geom;
	depth = other.depth;
	root = other.root;
}

AABBTree::AABBTree(Geometry* geom, PrimitiveType type)
{
	this->type = type;
	this->primitives = geom->getPrimitives(type);
	this->bbox = geom->getBoundingBox();
	this->bbox.min.addScalar(-0.001);
	this->bbox.max.addScalar(0.001);

	this->geom = geom;

	// generate index array to primitives
	std::vector<uint> primitiveIds;
	primitiveIds.reserve(this->primitives.size());
	uint n(0); std::generate_n(std::back_inserter(primitiveIds), this->primitives.size(), [n]() mutable { return n++;  });

	// generate root and all that follow
	this->root = new AABBNode(&primitiveIds, this->bbox, this);
}

AABBTree::~AABBTree()
{
}

bool AABBTree::hasPrimitives()
{
	return this->primitives.size() > 0;
}

bool AABBTree::boxIntersectsAPrimitive(Box3* box)
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
			for (auto&& pId:item->primitiveIds) {
				if (getPrimitiveBoxIntersection(this->primitives.at(pId), &center, &box->min, &box->max, &halfSize, 0.0f)) {
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

float AABBTree::boxIntersectsAPrimitiveAtDistance(Box3* box)
{
	std::stack<AABBNode*> stack = {};
	stack.push(this->root);

	while (stack.size()) {
		AABBNode* item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
			Vector3 center = box->getCenter();
			Vector3 halfSize = 0.5 * box->getSize();

			for (auto&& pId : item->primitiveIds) {
				if (getPrimitiveBoxIntersection(this->primitives.at(pId), &center, &box->min, &box->max, &halfSize, 0.0f)) {
					return getDistanceToAPrimitiveSq(this->primitives.at(pId), center);
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

// returns all primitives that are intersecting overlapping AABB leaves
void AABBTree::getPrimitivesInABox(Box3* box, std::vector<uint>* primIdBuffer)
{
	if (!primIdBuffer->empty()) primIdBuffer->clear();

	std::stack<AABBNode*> stack = {};
	stack.push(this->root);

	while (stack.size()) {
		AABBNode* item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
			for (uint i = 0; i < item->primitiveIds.size(); i++) {
				primIdBuffer->push_back(item->primitiveIds[i]);
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
			bool leftIsNear = point.getCoordById(item->axis) < item->splitPosition;
			AABBNode* nearNode = item->left;
			AABBNode* farNode = item->right;
			if (!leftIsNear) {
				nearNode = item->right;
				farNode = item->left;
			}

			if (nearNode) {
				stack.push(nearNode);
			}			
		}
		else {
			return item;
		}
	}

	return nullptr;
}

int AABBTree::getClosestPrimitiveId(Vector3& point)
{
	AABBNode* closestNode = this->getClosestNode(point);
	if (!closestNode) {
		return -1;
	}

	int id; float distSq; float result_distSq = FLT_MAX;
	for (auto&& i : closestNode->primitiveIds) {
		distSq = getDistanceToAPrimitiveSq(this->primitives[i], point);

		if (distSq < result_distSq) {
			result_distSq = distSq;
			id = i;
		}
	}
	return id;
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

std::vector<Geometry> AABBTree::getAABBPrimitivesOfDepth(uint depth)
{
	std::vector<AABBNode> leafNodes = this->flattenToDepth(depth);
	std::vector<Geometry> primGeoms = {};

	for (auto&& leaf : leafNodes) {
		std::vector<Geometry> leafTriGeoms = {};

		for (uint i = 0; i < leaf.primitiveIds.size(); i++) {
			Geometry primGeom = Geometry();
			primGeom.uniqueVertices = std::vector<Vector3>(
				*this->primitives[leaf.primitiveIds[i]].vertices.begin(),
				*this->primitives[leaf.primitiveIds[i]].vertices.end()
			);

			for (uint j = 0; j < 3; j++) {
				primGeom.vertices.push_back(this->primitives[leaf.primitiveIds[i]].vertices[j]->x);
				primGeom.vertices.push_back(this->primitives[leaf.primitiveIds[i]].vertices[j]->y);
				primGeom.vertices.push_back(this->primitives[leaf.primitiveIds[i]].vertices[j]->z);

				// normals will be computed when merging geoms

				primGeom.vertexIndices.push_back(j);
			}
			primGeoms.push_back(primGeom);
		}
	}

	return primGeoms;
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
	primitiveIds = std::vector<uint>(other.primitiveIds);
}

AABBTree::AABBNode::AABBNode(std::vector<uint>* primitiveIds, Box3& bbox, AABBTree* tree, uint depthLeft, AABBNode* parent)
{
	this->tree = tree;
	this->parent = parent;
	this->bbox = bbox;
	this->depth = MAX_DEPTH - depthLeft;

	this->construct(primitiveIds, depthLeft);
}

AABBTree::AABBNode::~AABBNode()
{
}

bool AABBTree::AABBNode::isALeaf()
{
	return (!this->left && !this->right);
}

bool AABBTree::AABBNode::isALeafWithPrimitives()
{
	return this->isALeaf() && !this->primitiveIds.empty();
}

void AABBTree::AABBNode::construct(std::vector<uint>* primitiveIds, uint depthLeft)
{
	if (depthLeft == 0 || primitiveIds->size() <= 2) {
		this->primitiveIds = *primitiveIds;
		this->filterPrimitives();
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

	std::vector<uint> leftPrimitiveIds = {};
	std::vector<uint> rightPrimitiveIds = {};

	this->splitPosition = this->getSplitPosition(*primitiveIds, &leftPrimitiveIds, &rightPrimitiveIds); // classify primitives

	const bool stopBranching = !this->hasEnoughBranching(leftPrimitiveIds.size(), rightPrimitiveIds.size(), primitiveIds->size());
	if (stopBranching) {
		this->primitiveIds = *primitiveIds;
		this->filterPrimitives();
	}
	else {
		if (leftPrimitiveIds.size() > 0) {
			Box3 bboxLeft = this->bbox;
			bboxLeft.max.setCoordById(this->splitPosition, this->axis);
			this->left = new AABBNode(&leftPrimitiveIds, bboxLeft, this->tree, depthLeft - 1, this);

			if (this->left->primitiveIds.empty() && this->left->left == nullptr && this->left->right == nullptr) {
				this->left = nullptr;
			}
		}
		if (rightPrimitiveIds.size() > 0) {
			Box3 bboxRight = this->bbox;
			bboxRight.min.setCoordById(this->splitPosition, this->axis);
			this->right = new AABBNode(&rightPrimitiveIds, bboxRight, this->tree, depthLeft - 1, this);

			if (this->right->primitiveIds.empty() && this->right->left == nullptr && this->right->right == nullptr) {
				this->right = nullptr;
			}
		}
	}
}

float AABBTree::AABBNode::getSplitPosition(std::vector<uint>& primitiveIds, std::vector<uint>* out_left, std::vector<uint>* out_right)
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
	for (uint i = 0; i < primitiveIds.size(); i++) {
		float min = this->tree->primitives[primitiveIds[i]].getMinById(axis);
		float max = this->tree->primitives[primitiveIds[i]].getMaxById(axis);

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
	for (uint i = 0; i < primitiveIds.size(); i++) {
		float min = this->tree->primitives[primitiveIds[i]].getMinById(axis);
		float max = this->tree->primitives[primitiveIds[i]].getMaxById(axis);

		if (min <= bestSplitPosition) {
			out_left->push_back(primitiveIds[i]);
		}
		if (max >= bestSplitPosition) {
			out_right->push_back(primitiveIds[i]);
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

bool AABBTree::AABBNode::hasEnoughBranching(size_t nLeftPrims, size_t nRightPrims, size_t nPrims)
{
	return nLeftPrims + nRightPrims < 1.5f * nPrims;
}

void AABBTree::AABBNode::filterPrimitives()
{
	Vector3 bboxCenter = bbox.getCenter();
	Vector3 bboxHalfSize = 0.5f * bbox.getSize();
	std::vector<uint> newPrimitiveIds = {};
	for (uint i = 0; i < this->primitiveIds.size(); i++) {
		if (getPrimitiveBoxIntersection(this->tree->primitives.at(this->primitiveIds[i]), &bboxCenter, &bbox.min, &bbox.max, &bboxHalfSize)) {
			newPrimitiveIds.push_back(this->primitiveIds[i]);
		}
	}
	this->primitiveIds = newPrimitiveIds;
}