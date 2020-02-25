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
	uint n(0); std::generate_n(std::back_inserter(primitiveIds), this->primitives.size(), [n]() mutable { return n++; });

	// generate root and all that follow
	this->root = new AABBNode(&primitiveIds, this->bbox, this);
}

AABBTree::~AABBTree()
{
	delete this->root;
	this->primitives.clear();
}

uint AABBTree::getDepth()
{
	if (this->root == nullptr) {
		return 0;
	}

	std::list<AABBTree::AABBNode*> queue;
	queue.push_back(this->root);

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

bool AABBTree::hasPrimitives()
{
	return this->primitives.size() > 0;
}

bool AABBTree::boxIntersectsAPrimitive(Box3* box)
{
	std::stack<AABBNode*> stack = {};
	stack.push(this->root);
	AABBNode* item;
	Vector3 center; Vector3 halfSize;
	box->setToCenter(&center);
	box->setToHalfSize(&halfSize);

	while (stack.size()) {
		item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
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
	AABBNode* item;
	Vector3 center; Vector3 halfSize;
	box->setToCenter(&center);
	box->setToHalfSize(&halfSize);

	while (stack.size()) {
		item = stack.top();
		stack.pop();

		if (item->isALeaf()) {
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

uint AABBTree::rayIntersect(Vector3& rayOrigin, Vector3& rayDirection, float rayMinParam, float rayMaxParam)
{
	uint hitCount = 0;
	Vector3 rayStart = rayOrigin;
	Vector3 rayInvDirection = Vector3(1.0f / rayDirection.x, 1.0f / rayDirection.y, 1.0f / rayDirection.z);
	AABBRay ray = AABBRay(&rayOrigin, &rayDirection);

	float hitParam;
	bool treeBoxHit = getRayBoxIntersection(ray, &bbox.min, &bbox.max, &hitParam);
	if (!treeBoxHit) {
		return hitCount;
	}

	AABBNode* item;	AABBNode* nearNode;	AABBNode* farNode;
	float t, t_split;
	bool leftIsNear, force_near, force_far, t_splitAtInfinity;


	std::stack<AABBNode*> stack = {};
	stack.push(this->root);
	std::set<float> t_set = {};

	while (stack.size()) {
		item = stack.top();
		stack.pop();


		if (!item->isALeaf()) {
			leftIsNear = rayOrigin.getCoordById(item->axis) < item->splitPosition;
			nearNode = item->left; farNode = item->right;
			if (!leftIsNear) {
				nearNode = item->right;
				farNode = item->left;
			}

			t_split = (item->splitPosition - rayStart.getCoordById(item->axis)) * rayInvDirection.getCoordById(item->axis);
			force_near = false; force_far = false;
			t_splitAtInfinity = t_split > FLT_MAX || t_split < -FLT_MAX;
			if (t_splitAtInfinity) {
				if (rayStart.getCoordById(item->axis) <= item->splitPosition &&
					rayStart.getCoordById(item->axis) >= item->bbox.min.getCoordById(item->axis)) {
					if (leftIsNear) force_near = true;
					else force_far = true;
				}
				if (rayStart.getCoordById(item->axis) >= item->splitPosition &&
					rayStart.getCoordById(item->axis) <= item->bbox.max.getCoordById(item->axis)) {
					if (leftIsNear) force_near = true;
					else force_far = true;
				}
			}

			if (farNode && ((!t_splitAtInfinity && rayMaxParam >= t_split) || force_far)) {
				stack.push(farNode);
			}
			if (nearNode && ((!t_splitAtInfinity && rayMinParam <= t_split) || force_near)) {
				stack.push(nearNode);
			}
		}
		else {
			for (auto&& pId : item->primitiveIds) {
				/*
				if (pId == 10) {
					t = 0.0f;
				}*/
				t = getRayTriangleIntersection(rayStart, rayDirection, &primitives[pId].vertices, rayMinParam, rayMaxParam);
				if (t > 0) {
					hitCount++;
					// t += (bbox.max.x - bbox.min.x) / 1000.0f;
					// rayStart.set(rayStart.x + t * rayDirection.x, rayStart.y + t * rayDirection.y, rayStart.z + t * rayDirection.z);
					t_set.insert(t);
				}
			}

		}
	}

	return t_set.size();
}

bool AABBTree::boolRayIntersect(Vector3& rayOrigin, Vector3& rayDirection, float rayMinParam, float rayMaxParam)
{
	Vector3 rayStart = rayOrigin;
	Vector3 rayInvDirection = Vector3(1.0f / rayDirection.x, 1.0f / rayDirection.y, 1.0f / rayDirection.z);
	AABBRay ray = AABBRay(&rayOrigin, &rayDirection);

	float hitParam;
	bool treeBoxHit = getRayBoxIntersection(ray, &bbox.min, &bbox.max, &hitParam);
	if (!treeBoxHit) {
		return false;
	}

	AABBNode* item;	AABBNode* nearNode;	AABBNode* farNode;
	float t, t_split;
	bool leftIsNear, force_near, force_far, t_splitAtInfinity;


	std::stack<AABBNode*> stack = {};
	stack.push(this->root);
	std::set<float> t_set = {};

	while (stack.size()) {
		item = stack.top();
		stack.pop();


		if (!item->isALeaf()) {
			leftIsNear = rayOrigin.getCoordById(item->axis) < item->splitPosition;
			nearNode = item->left; farNode = item->right;
			if (!leftIsNear) {
				nearNode = item->right;
				farNode = item->left;
			}

			t_split = (item->splitPosition - rayStart.getCoordById(item->axis)) * rayInvDirection.getCoordById(item->axis);
			force_near = false; force_far = false;
			t_splitAtInfinity = t_split > FLT_MAX || t_split < -FLT_MAX;
			if (t_splitAtInfinity) {
				if (rayStart.getCoordById(item->axis) <= item->splitPosition &&
					rayStart.getCoordById(item->axis) >= item->bbox.min.getCoordById(item->axis)) {
					if (leftIsNear) force_near = true;
					else force_far = true;
				}
				if (rayStart.getCoordById(item->axis) >= item->splitPosition &&
					rayStart.getCoordById(item->axis) <= item->bbox.max.getCoordById(item->axis)) {
					if (leftIsNear) force_near = true;
					else force_far = true;
				}
			}

			if (farNode && ((!t_splitAtInfinity && rayMaxParam >= t_split) || force_far)) {
				stack.push(farNode);
			}
			if (nearNode && ((!t_splitAtInfinity && rayMinParam <= t_split) || force_near)) {
				stack.push(nearNode);
			}
		}
		else {

			for (auto&& pId : item->primitiveIds) {
				t = getRayTriangleIntersection(rayStart, rayDirection, &primitives[pId].vertices, rayMinParam, rayMaxParam);
				if (t > 0) {
					return true;
				}
			}
		}
	}

	return false;
}

std::vector<AABBTree::AABBNode> AABBTree::flatten()
{
	std::vector<AABBNode> resultArray = {};

	std::stack<AABBNode*> nodeStack = {};
	nodeStack.push(this->root);
	AABBNode* item;

	while (nodeStack.size()) {
		item = nodeStack.top();
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

void AABBTree::getAllNodes(std::vector<AABBNode>* buffer)
{
	std::stack<AABBNode*> nodeStack;
	nodeStack.push(this->root);

	while (nodeStack.size()) {
		AABBNode* item = nodeStack.top();
		nodeStack.pop();

		buffer->push_back(*item);

		if (item->left) {
			nodeStack.push(item->left);
		}
		if (item->right) {
			nodeStack.push(item->right);
		}
	}
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
	AABBNode* item;

	while (stack.size()) {
		item = stack.top();
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
	AABBNode* item; AABBNode* nearNode; AABBNode* farNode;
	bool leftIsNear;

	while (stack.size()) {
		item = stack.top();
		stack.pop();

		if (item->left || item->right) {
			leftIsNear = point.getCoordById(item->axis) < item->splitPosition;
			nearNode = item->left;
			farNode = item->right;
			if (!leftIsNear) {
				nearNode = item->right;
				farNode = item->left;
			}

			if (nearNode) {
				stack.push(nearNode);
			}
			else if (farNode) {
				stack.push(farNode);
			}
		}
		else {
			return item;
		}
	}

	return nullptr;
}

int AABBTree::getClosestPrimitiveIdAndDist(Vector3& point, float* result)
{
	AABBNode* closestNode = this->getClosestNode(point);
	if (!closestNode) {
		return -1;
	}

	int id = -1; float distSq; *result = FLT_MAX;
	for (auto&& i : closestNode->primitiveIds) {
		distSq = getDistanceToAPrimitiveSq(this->primitives[i], point);

		if (distSq < *result) {
			*result = distSq;
			id = i;
		}
	}

	return id;
}

void AABBTree::applyMatrix(Matrix4& m)
{
	// this->geom->applyMatrix(m);
	this->bbox.min.applyMatrix4(m);
	this->bbox.max.applyMatrix4(m);

	std::stack<AABBNode*> nodeStack;
	nodeStack.push(this->root);
	AABBNode* item;

	while (nodeStack.size()) {
		item = nodeStack.top();
		nodeStack.pop();

		item->applyMatrix(m);

		if (item->left) {
			nodeStack.push(item->left);
		}
		if (item->right) {
			nodeStack.push(item->right);
		}
	}
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
			uint Nverts = this->primitives[leaf.primitiveIds[i]].vertices.size();
			for (uint k = 0; k < Nverts; k++) {
				primGeom.uniqueVertices.push_back(*this->primitives[leaf.primitiveIds[i]].vertices[k]);

				primGeom.vertices.push_back(this->primitives[leaf.primitiveIds[i]].vertices[k]->x);
				primGeom.vertices.push_back(this->primitives[leaf.primitiveIds[i]].vertices[k]->y);
				primGeom.vertices.push_back(this->primitives[leaf.primitiveIds[i]].vertices[k]->z);

				primGeom.vertexIndices.push_back(k);
			}
			primGeoms.push_back(primGeom);
		}
	}

	return primGeoms;
}

void AABBTree::GenerateFullTreeBoxVisualisation(VTKExporter& e)
{
	std::cout << "--------------------------------------------" << std::endl;
	std::vector<Geometry> boxGeoms = {};
	std::vector<AABBNode> nodeBuffer = {};
	std::cout << "obtaining AABB nodes..." << std::endl;
	this->getAllNodes(&nodeBuffer);
	std::cout << nodeBuffer.size() << " AABB nodes retrieved" << std::endl;

	std::cout << "generating box geometries..." << std::endl;
	for (auto&& n : nodeBuffer) {
		float dimX = n.bbox.max.x - n.bbox.min.x;
		float dimY = n.bbox.max.y - n.bbox.min.y;
		float dimZ = n.bbox.max.z - n.bbox.min.z;
		PrimitiveBox box = PrimitiveBox(dimX, dimY, dimZ, 1, 1, 1);
		Vector3 t = n.bbox.min;
		box.applyMatrix(Matrix4().makeTranslation(t.x, t.y, t.z));
		boxGeoms.push_back(box);
	}

	Geometry resultGeom = mergeGeometries(boxGeoms);
	e.initExport(&resultGeom, geom->name + "_AABB_allBoxes");
	std::cout << "AABB box geometries retrieved and exported" << std::endl;
}

void AABBTree::GenerateFullLeafBoxVisualisation(VTKExporter& e)
{
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "obtaining AABB leaf boxes..." << std::endl;
	std::vector<Geometry> boxGeoms = getAABBLeafGeoms();
	std::cout << boxGeoms.size() << " AABB leaf boxes retrieved" << std::endl;
	Geometry resultGeom = mergeGeometries(boxGeoms);
	e.initExport(&resultGeom, geom->name + "_AABB_leafBoxes");
	std::cout << "AABB leaf boxes exported" << std::endl;
}

uint depth(AABBTree* tree)
{
	return tree->getDepth();
}

bool getRayBoxIntersection(const AABBTree::AABBRay& r, Vector3* boxMin, Vector3* boxMax, float* hitParam)
{
	float tmin, tmax, tymin, tymax, tzmin, tzmax;
	bool signX = r.invDirection.x < 0;
	bool signY = r.invDirection.y < 0;
	bool signZ = r.invDirection.z < 0;
	float b0X = (signX ? boxMax->x : boxMin->x);
	float b0Y = (signY ? boxMax->y : boxMin->y);
	float b0Z = (signZ ? boxMax->z : boxMin->z);
	float b1X = (signX ? boxMin->x : boxMax->x);
	float b1Y = (signY ? boxMin->y : boxMax->y);
	float b1Z = (signZ ? boxMin->z : boxMax->z);
	
	tmin = (b0X - r.start.x) * r.invDirection.x;
	tmax = (b1X - r.start.x) * r.invDirection.x;
	tymin = (b0Y - r.start.y) * r.invDirection.y;
	tymax = (b1Y - r.start.y) * r.invDirection.y;

	if ((tmin > tymax) || (tymin > tmax)) return false;
	if (tymin > tmin) tmin = tymin;
	if (tymax < tmax) tmax = tymax;

	tzmin = (b0Z - r.start.z) * r.invDirection.z;
	tzmax = (b1Z - r.start.z) * r.invDirection.z;

	if ((tmin > tzmax) || (tzmin > tmax)) return false;
	if (tzmin > tmin) tmin = tzmin;
	if (tzmax < tmax) tmax = tzmax;

	if ((tmin < r.maxParam) && (tmax > r.minParam)) {
		*hitParam = tmax;
		return true;
	}
	return false;
}

void AABBTree::GenerateStepwiseLeafBoxVisualisation(VTKExporter& e)
{
	std::cout << "--------------------------------------------" << std::endl;
	uint t_depth = this->getDepth();
	std::cout << "generating AABB leaf visualisation of depth = 0, ..., " << t_depth << std::endl;

	for (uint d = 0; d < t_depth; d++) {
		std::cout << "d = " << d << "..." << std::endl;
		std::vector<Geometry> boxGeoms = this->getAABBGeomsOfDepth(d);
		std::vector<Geometry> triGeoms = this->getAABBPrimitivesOfDepth(d);

		Geometry resultBoxGeom = mergeGeometries(boxGeoms);
		Geometry resultTriGeom = mergeGeometries(triGeoms);

		e.initExport(&resultBoxGeom, geom->name + "_AABB_leafBB_step-" + std::to_string(d));
		e.initExport(&resultTriGeom, geom->name + "_AABB_leafTris_step-" + std::to_string(d));
	}
	std::cout << "AABB leaf visualisations finished" << std::endl;
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

	// choosing the longest split axis seems to be faster
	if ((bbox.max.x - bbox.min.x) > (bbox.max.y - bbox.min.y) && (bbox.max.x - bbox.min.x) > (bbox.max.z - bbox.min.z)) {
		this->axis = 0;
	}
	else if ((bbox.max.y - bbox.min.y) > (bbox.max.z - bbox.min.z)) {
		this->axis = 1;
	}
	else {
		this->axis = 2;
	}

	std::vector<uint> leftPrimitiveIds = {};
	std::vector<uint> rightPrimitiveIds = {};

	this->splitPosition = this->getAdaptivelyResampledSplitPosition(*primitiveIds, &leftPrimitiveIds, &rightPrimitiveIds); // classify primitives

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
	const uint CUTS = 8; uint i, j;
	float* splitPositions = new float[CUTS];
	float* splitsLeft = new float[CUTS];
	float* splitsRight = new float[CUTS];

	for (i = 0; i < CUTS; i++) {
		splitPositions[i] = 
			bbox.min.getCoordById(axis) * (1.0f - ((float)(i + 1) / (float)(CUTS + 1.0f))) +
			bbox.max.getCoordById(axis) * ((float)(i + 1) / (float)(CUTS + 1.0f));
		splitsLeft[i] = 0.0f;
		splitsRight[i] = 0.0f;
	}

	float min, max;

	// count triangles for each split, we don't want to fill arrays here (lot of wasted cycles and memory ops)
	for (i = 0; i < primitiveIds.size(); i++) {
		min = this->tree->primitives[primitiveIds[i]].getMinById(axis);
		max = this->tree->primitives[primitiveIds[i]].getMaxById(axis);

		for (j = 0; j < CUTS; j++) {
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
	for (i = 0; i < CUTS; i++) {
		float cost = getCostEstimate(splitPositions[i], splitsLeft[i], splitsRight[i]);
		if (cost < bestCost) {
			bestCost = cost;
			bestSplitPosition = splitPositions[i];
		}
	}

	// fill left and right arrays now that best split position is known
	for (i = 0; i < primitiveIds.size(); i++) {
		float min = this->tree->primitives[primitiveIds[i]].getMinById(axis);
		float max = this->tree->primitives[primitiveIds[i]].getMaxById(axis);

		if (min <= bestSplitPosition) {
			out_left->push_back(primitiveIds[i]);
		}
		if (max >= bestSplitPosition) {
			out_right->push_back(primitiveIds[i]);
		}
	}

	delete[] splitsLeft;
	delete[] splitsRight;
	delete[] splitPositions;

	return bestSplitPosition;
}

// ============ Adaptive resampling helper macros ============

#define SHIFT_m256_LEFT(source, target)												\
		t0 = _mm256_permute_ps(source, 0x39);										\
		t1 = _mm256_permute2f128_ps(t0, t0, 0x01);									\
		target = _mm256_blend_ps(t0, t1, 0x88);

// ---------------------------------------------------------

#define SHIFT_m256_RIGHT(source, target)											\
		t0 = _mm256_permute_ps(source, _MM_SHUFFLE(2, 1, 0, 3));					\
		t1 = _mm256_permute2f128_ps(t0, t0, 41);									\
		target = _mm256_blend_ps(t0, t1, 0x11);

// --------------------------------------------------------

#define BOX_AREA(b)																	\
		2.0f * (b.max.x - b.min.x) +												\
		2.0f * (b.max.y - b.min.y) +												\
		2.0f * (b.max.z - b.min.z);

// ---------------------------------------------------------

#define GET_QUAD_MIN_BETWEEN_PTS(x0, y0, x1, y1, x2, y2, target)					\
		det = (x0 - x1) * (x0 - x2) * (x1 - x2);									\
		det0 = x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1);					\
		det1 = x2 * x2 * (y0 - y1) + x0 * x0 * (y1 - y2) + x1 * x1 * (y2 - y0);		\
		a2 = det0 / det; a1 = det1 / det;											\
		target = -a1 / (2.0f * a2);

// --------------------------------------------------------

// Fast kd-tree Construction with an Adaptive Error-Bounded Heuristic (Hunt, Mark, Stoll)
//
float AABBTree::AABBNode::getAdaptivelyResampledSplitPosition(std::vector<uint>& primitiveIds, std::vector<uint>* out_left, std::vector<uint>* out_right)
{
	const float tot_BoxArea = BOX_AREA(bbox);
	const uint CUTS = 4; uint i, j;

	// === Stage 1: Initial sampling of C_L(x) and C_R(x) ========================================================

	// split interval limits:
	const float a = bbox.min.getCoordById(axis);
	const float b = bbox.max.getCoordById(axis);

	// splits:
	union { __m256 B_min_R; float BminR[6]; };
	B_min_R = _mm256_setr_ps(
		a,
		a * (1.0f - (1.0f / (CUTS + 1))) + b * (1.0f / (CUTS + 1)),
		a * (1.0f - (2.0f / (CUTS + 1))) + b * (2.0f / (CUTS + 1)),
		a * (1.0f - (3.0f / (CUTS + 1))) + b * (3.0f / (CUTS + 1)),
		a * (1.0f - (4.0f / (CUTS + 1))) + b * (4.0f / (CUTS + 1)),
		b,
		0.0f, 0.0f
	);

	// set C_L(x) = 0, C_R(x) = 0
	union { __m128 C_L; float cl[4]; };
	union { __m128 C_R; float cr[4]; };
	C_L = _mm_setzero_ps();
	C_R = _mm_setzero_ps();

	// masks for C_L and C_R increments:
	union { __m128 mask_L; uint m_Li[4]; };
	union { __m128 mask_R; uint m_Ri[4]; };
	__m128 r_L, r_R;
	// shift B_min_R to the left to avoid applying mask to split at x = a:
	union { __m128 B_min_shift; float bmn[4]; };
	__m256 s_res; __m256 t0, t1;
	SHIFT_m256_LEFT(B_min_R, s_res);
	B_min_shift = _mm256_castps256_ps128(s_res);

	uint N_primitives = primitiveIds.size();
	float* primMins = new float[N_primitives];
	float* primMaxes = new float[N_primitives];

	for (i = 0; i < N_primitives; i++) {
		primMins[i] = this->tree->primitives[primitiveIds[i]].getMinById(axis);
		primMaxes[i] = this->tree->primitives[primitiveIds[i]].getMaxById(axis);

		mask_L = _mm_cmple_ps(_mm_set1_ps(primMins[i]), B_min_shift);
		mask_R = _mm_cmpge_ps(_mm_set1_ps(primMaxes[i]), B_min_shift);
		r_L = _mm_setr_ps((m_Li[0] > 0) * 1.0f,	(m_Li[1] > 0) * 1.0f, (m_Li[2] > 0) * 1.0f,	(m_Li[3] > 0) * 1.0f);
		r_R = _mm_setr_ps((m_Ri[0] > 0) * 1.0f,	(m_Ri[1] > 0) * 1.0f, (m_Ri[2] > 0) * 1.0f,	(m_Ri[3] > 0) * 1.0f);
		C_L = _mm_add_ps(C_L, r_L);
		C_R = _mm_add_ps(C_R, r_R);
	}

	// ===== Stage 2: Sample range [0, N_primitives] uniformly & count the number of samples within each segment ======

	// range [0, N_primitives] sampling
	union { __m256i S_L; uint sl[5]; };
	union { __m256i S_R; uint sr[5]; };
	float ran_s, cMin, cMax;

	// segment samples have to be counted here:
	// extended ranges
	union { __m256 C_L_ext; float cl_ext[6]; };
	union { __m256 C_R_ext; float cr_ext[6]; };
	C_L_ext = _mm256_castps128_ps256(C_L);
	SHIFT_m256_RIGHT(C_L_ext, C_L_ext); cl_ext[5] = 1.0f * N_primitives;
	C_R_ext = _mm256_castps128_ps256(C_R);
	SHIFT_m256_RIGHT(C_R_ext, C_R_ext); cr_ext[0] = 1.0f * N_primitives;
	// store in reverse since C_R(x) is non-increasing
	C_R_ext = _mm256_setr_ps(cr_ext[5], cr_ext[4], cr_ext[3], cr_ext[2], cr_ext[1], cr_ext[0], 0.0f, 0.0f);

	S_L = _mm256_setzero_si256();
	S_R = _mm256_setzero_si256();

	// compare masks for sampling interval bounds
	union { __m256 cmpMin_SL; uint cmSL[6]; };
	union { __m256 cmpMin_SR; uint cmSR[6]; };
	union { __m256 cmpMax_SL; uint cmxSL[6]; };
	union { __m256 cmpMax_SR; uint cmxSR[6]; };

	for (i = 0; i < CUTS; i++) {
		ran_s = (float)(i + 1) / (float)(CUTS + 1) * N_primitives;

		t0 = _mm256_set1_ps(ran_s);
		cmpMin_SL = _mm256_cmp_ps(t0, C_L_ext, _CMP_GT_OS);
		cmpMin_SR = _mm256_cmp_ps(t0, C_R_ext, _CMP_GT_OS);
		cmpMax_SL = _mm256_cmp_ps(t0, C_L_ext, _CMP_LT_OS);
		cmpMax_SR = _mm256_cmp_ps(t0, C_R_ext, _CMP_LT_OS);
		S_L = _mm256_add_epi32(S_L,
			_mm256_setr_epi32(
				(cmSL[0] > 0 && cmxSL[1] > 0) * 1,
				(cmSL[1] > 0 && cmxSL[2] > 0) * 1,
				(cmSL[2] > 0 && cmxSL[3] > 0) * 1,
				(cmSL[3] > 0 && cmxSL[4] > 0) * 1,
				(cmSL[4] > 0 && cmxSL[5] > 0) * 1,
				0, 0, 0));

		S_R = _mm256_add_epi32(S_R,
			_mm256_setr_epi32(
				(cmSR[0] > 0 && cmxSR[1] > 0) * 1,
				(cmSR[1] > 0 && cmxSR[2] > 0) * 1,
				(cmSR[2] > 0 && cmxSR[3] > 0) * 1,
				(cmSR[3] > 0 && cmxSR[4] > 0) * 1,
				(cmSR[4] > 0 && cmxSR[5] > 0) * 1,
				0, 0, 0));
	}
	S_R = _mm256_setr_epi32(sr[4], sr[3], sr[2], sr[1], sr[0], 0, 0, 0);

	// ==== Stage 3: add more sampling positions to subdivided segments ===========================================

	union { __m256 all_samplePos_L; float allspl[8]; };
	union { __m256 all_samplePos_R; float allspr[8]; };
	uint nSeg_L = 0, nSeg_R = 0;
	float segLen = (b - a) / 5.0f;

	for (i = 0; i <= CUTS; i++) {
		if (i > 0) {
			allspl[nSeg_L++] = BminR[i];
			allspr[nSeg_R++] = BminR[i];
		}
		for (j = 0; j < sl[i]; j++) {
			allspl[nSeg_L++] = BminR[i] + (float)(j + 1) / (float)(sl[i] + 1) * segLen;
		}
		for (j = 0; j < sr[i]; j++) {
			allspr[nSeg_R++] = BminR[i] + (float)(j + 1) / (float)(sr[i] + 1) * segLen;
		}
	}

	// Compute surface area heuristic SAH:
	// remaining two dimensions of the child box candidates
	float boxDim0 = this->bbox.max.getCoordById((axis + 1) % 3) - this->bbox.min.getCoordById((axis + 1) % 3);
	float boxDim1 = this->bbox.max.getCoordById((axis + 2) % 3) - this->bbox.min.getCoordById((axis + 2) % 3);
	__m256 SA_L, SA_R;

	// SA_L(x) = (boxDim_L(x) + boxDim0 + boxDim1) * 2.0f / tot_BoxArea
	// SA_R(x) = (boxDim_R(x) + boxDim0 + boxDim1) * 2.0f / tot_BoxArea
	SA_L = _mm256_mul_ps(_mm256_add_ps(_mm256_sub_ps(all_samplePos_L, _mm256_set1_ps(a)), _mm256_add_ps(_mm256_set1_ps(boxDim0), _mm256_set1_ps(boxDim1))), _mm256_set1_ps(2.0f / tot_BoxArea));
	SA_R = _mm256_mul_ps(_mm256_add_ps(_mm256_sub_ps(_mm256_set1_ps(b), all_samplePos_R), _mm256_add_ps(_mm256_set1_ps(boxDim0), _mm256_set1_ps(boxDim1))), _mm256_set1_ps(2.0f / tot_BoxArea));

	// ==== Stage 4: RESAMPLE C_L and C_R on all sample points & construct a piecewise quadratic approximation of cost(x) to minimize

	union { __m256 C_L_final; float cl_fin[8]; };
	union { __m256 C_R_final; float cr_fin[8]; };
	C_L_final = _mm256_setzero_ps();
	C_R_final = _mm256_setzero_ps();
	// masks for C_L and C_R increments:
	union { __m256 mask_CL; uint mcl[8]; };
	union { __m256 mask_CR; uint mcr[8]; };
	__m256 rs_L, rs_R;
	float min, max;

	for (i = 0; i < N_primitives; i++) {
		min = primMins[i];
		max = primMaxes[i];

		mask_CL = _mm256_cmp_ps(_mm256_set1_ps(min), all_samplePos_L, _CMP_LT_OS);
		mask_CR = _mm256_cmp_ps(_mm256_set1_ps(max), all_samplePos_R, _CMP_GT_OS);
		rs_L = _mm256_setr_ps(
			(mcl[0] > 0) * 1.0f, (mcl[1] > 0) * 1.0f, (mcl[2] > 0) * 1.0f, (mcl[3] > 0) * 1.0f,
			(mcl[4] > 0) * 1.0f, (mcl[5] > 0) * 1.0f, (mcl[6] > 0) * 1.0f, (mcl[7] > 0) * 1.0f
		);
		rs_R = _mm256_setr_ps(
			(mcr[0] > 0) * 1.0f, (mcr[1] > 0) * 1.0f, (mcr[2] > 0) * 1.0f, (mcr[3] > 0) * 1.0f,
			(mcr[4] > 0) * 1.0f, (mcr[5] > 0) * 1.0f, (mcr[6] > 0) * 1.0f, (mcr[7] > 0) * 1.0f
		);
		C_L_final = _mm256_add_ps(C_L_final, rs_L);
		C_R_final = _mm256_add_ps(C_R_final, rs_R);
	}

	// cost(x) = C_L(x) * SA_L(x) + C_R(x) * SA_R(x):
	union { __m256 COST; float cost[8]; };
	COST = _mm256_add_ps(_mm256_mul_ps(C_L_final, SA_L), _mm256_mul_ps(C_R_final, SA_R));

	// ==== Stage 5: Minimize cost(x) & classify primitives  =============================================================

	/*
	// Alternative I: Quad interpolate argmin of cost(x) between triplets of sampled positions when when triplets form a convex parabolic arc
	//
	// extended discretisation domain including a and b
	float dD[10] = {a, allspl[0], allspl[1], allspl[2], allspl[3], allspl[4], allspl[5], allspl[6], allspl[7], b};
	// extended cost sample vector also evaluated at a and b	
	float cost_ab[10] = {(float)N_primitives, cost[0], cost[1], cost[2], cost[3], cost[4], cost[5], cost[6], cost[7], (float)N_primitives};

	float bestSplit = a;
	float det, det0, det1, a1, a2; // used in quad interpolation

	for (i = 1; i <= 2 * CUTS; i++) {
		if ((cost_ab[i - 1] >= cost_ab[i]) && (cost_ab[i] <= cost_ab[i + 1])) {
			// interpolating min iff the leading pts form a concave parabola
			GET_QUAD_MIN_BETWEEN_PTS(
				dD[i - 1], cost_ab[i - 1],
				dD[i], cost_ab[i],
				dD[i + 1], cost_ab[i + 1], bestSplit
			);
			break;
		}
	}*/
	// Alternative I.a: implement sort order method for __m265 vector to obtain the indices of minimal splits
	// then, interpolate a parabola using the ids of 3 minimal indices

	/**/
	// Alternative II: simply find min cost by comparing vals - less exact, but
	// faster than previous by ~20%
	
	float bestSplit, minCost = FLT_MAX;
	for (i = 1; i < 2 * CUTS; i++) {
		if (cost[i] < minCost) {
			minCost = cost[i];
			bestSplit = allspl[i];
		}
	}
	

	// fill left and right arrays now that best split position is known:
	for (i = 0; i < N_primitives; i++) {
		min = primMins[i];
		max = primMaxes[i];

		if (min <= bestSplit) {
			out_left->push_back(primitiveIds[i]);
		}
		if (max >= bestSplit) {
			out_right->push_back(primitiveIds[i]);
		}
	}

	delete[] primMins;
	delete[] primMaxes;

	return bestSplit;
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

void AABBTree::AABBNode::applyMatrix(Matrix4& m)
{
	bbox.min.applyMatrix4(m);
	bbox.max.applyMatrix4(m);
}

AABBTree::AABBRay::AABBRay()
{
}

AABBTree::AABBRay::AABBRay(Vector3* start, Vector3* dir)
{
	this->start.set(start->x, start->y, start->z);
	if (dir) {
		this->direction.set(dir->x, dir->y, dir->z);
	}	
}
