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
	e.initExport(resultGeom, geom->name + "_AABB_allBoxes");
	std::cout << "AABB box geometries retrieved and exported" << std::endl;
}

void AABBTree::GenerateFullLeafBoxVisualisation(VTKExporter& e)
{
	std::cout << "--------------------------------------------" << std::endl;
	std::cout << "obtaining AABB leaf boxes..." << std::endl;
	std::vector<Geometry> boxGeoms = getAABBLeafGeoms();
	std::cout << boxGeoms.size() << " AABB leaf boxes retrieved" << std::endl;
	Geometry resultGeom = mergeGeometries(boxGeoms);
	e.initExport(resultGeom, geom->name + "_AABB_leafBoxes");
	std::cout << "AABB leaf boxes exported" << std::endl;
}

uint depth(AABBTree* tree)
{
	return tree->getDepth();
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

		e.initExport(resultBoxGeom, geom->name + "_AABB_leafBB_step-" + std::to_string(d));
		e.initExport(resultTriGeom, geom->name + "_AABB_leafTris_step-" + std::to_string(d));
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
	Vector3 bboxSize = this->bbox.getSize();
	if (bboxSize.x > bboxSize.y && bboxSize.x > bboxSize.z) {
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

#define FLERP(y0, y1, x0, x1, x)					\
		y0 + (x - x0) * (y1 - y0) / (x1 - x0);

float AABBTree::AABBNode::getAdaptivelyResampledSplitPosition(std::vector<uint>& primitiveIds, std::vector<uint>* out_left, std::vector<uint>* out_right)
{
	const uint CUTS = 4; uint i, j;
	float min, max;

	// === Stage 1: Initial sampling of C_L(x) and C_R(x)
	// split interval limits:
	float a = bbox.min.getCoordById(axis);
	float b = bbox.max.getCoordById(axis);

	// splits:
	union { __m256 B_min_R; float BminR[6]; };
	B_min_R = _mm256_setzero_ps();
	for (i = 0; i <= CUTS + 1; i++) {
		BminR[i] = a * (1.0f - ((float)i/ (float)(CUTS + 1.0f))) + b * ((float)i / (float)(CUTS + 1.0f));
	}
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
	__m256 t0 = _mm256_permute_ps(B_min_R, 0x39);
	__m256 t1 = _mm256_permute2f128_ps(t0, t0, 0x01);
	__m256 s_res = _mm256_blend_ps(t0, t1, 0x88);
	__m128 B_min_shift = _mm256_castps256_ps128(s_res);

	for (i = 0; i < primitiveIds.size(); i++) {
		min = this->tree->primitives[primitiveIds[i]].getMinById(axis);
		max = this->tree->primitives[primitiveIds[i]].getMaxById(axis);

		mask_L = _mm_cmple_ps(_mm_set1_ps(min), B_min_shift);
		mask_R = _mm_cmpge_ps(_mm_set1_ps(max), B_min_shift);
		r_L = _mm_setr_ps(
			(m_Li[0] > 0) * 1.0f,
			(m_Li[1] > 0) * 1.0f,
			(m_Li[2] > 0) * 1.0f,
			(m_Li[3] > 0) * 1.0f
		);
		r_R = _mm_setr_ps(
			(m_Ri[0] > 0) * 1.0f,
			(m_Ri[1] > 0) * 1.0f,
			(m_Ri[2] > 0) * 1.0f,
			(m_Ri[3] > 0) * 1.0f
		);
		C_L = _mm_add_ps(C_L, r_L);
		C_R = _mm_add_ps(C_R, r_R);
	}

	// ===== Stage 2: Sample range [0, N_primitives] & lerp sample inputs S_L, and S_R ======

	uint N_primitives = primitiveIds.size();
	// range [0, N_primitives] sampling
	union { __m128 S_L; float sl[4]; };
	union { __m128 S_R; float sr[4]; };
	float ran_s, cMax;
	for (i = 0; i < CUTS; i++) {
		ran_s = (float)(i + 1) / (float)(CUTS + 1) * N_primitives;

		cMax = (i + 1 < CUTS ? cl[i + 1] : 1.0f * N_primitives);
		sl[i] = FLERP(BminR[i], BminR[i + 1], cl[i], cMax, ran_s);
		cMax = (i + 1 < CUTS ? cr[i + 1] : 0.0f);
		sr[i] = FLERP(BminR[i], BminR[i + 1], cr[i], cMax, ran_s);
	}

	// ===== Stage 3: Counting samples in sample regions:
	// N_L += (BminR < S_L ? 1 : 0)
	// N_R += (BminR > S_R ? 1 : 0)
	__m128 N_L = _mm_setzero_ps();
	__m128 N_R = _mm_setzero_ps();

	mask_L = _mm_cmplt_ps(B_min_shift, S_L);
	mask_R = _mm_cmpgt_ps(B_min_shift, S_R);
	r_L = _mm_setr_ps(
		(m_Li[0] > 0) * 1.0f,
		(m_Li[1] > 0) * 1.0f,
		(m_Li[2] > 0) * 1.0f,
		(m_Li[3] > 0) * 1.0f
	);
	r_R = _mm_setr_ps(
		(m_Ri[0] > 0) * 1.0f,
		(m_Ri[1] > 0) * 1.0f,
		(m_Ri[2] > 0) * 1.0f,
		(m_Ri[3] > 0) * 1.0f
	);
	N_L = _mm_add_ps(N_L, r_L);
	N_R = _mm_add_ps(N_R, r_R);

	return 0.0f;
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
