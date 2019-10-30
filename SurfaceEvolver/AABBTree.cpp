#include "AABBTree.h"

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

bool getTriangleBoundingBoxIntersection(Tri& vertices, Vector3& bboxCenter, Vector3& bboxHalfSize, Vector3* optTriNormal = nullptr) {
	bboxHalfSize.addScalar(0.0001);

	std::vector<Vector3> verts = { vertices[0] - bboxCenter, vertices[1] - bboxCenter, vertices[2] - bboxCenter };

	Vector3 t_min = verts[0];
	t_min.min(verts[1]);
	t_min.min(verts[2]);

	Vector3 t_max = verts[0];
	t_max.max(verts[1]);
	t_max.max(verts[2]);

	bool aabb_overlap = (
		!(t_max.x < -bboxHalfSize.x || t_min.x > bboxHalfSize.x) &&
		!(t_max.y < -bboxHalfSize.y || t_min.y > bboxHalfSize.y) &&
		!(t_max.z < -bboxHalfSize.z || t_min.z > bboxHalfSize.z)
		);

	if (!aabb_overlap) {
		return false;
	}

	Vector3 n = optTriNormal != nullptr ? *optTriNormal : normalize(cross(verts[1] - verts[0], verts[2] - verts[0]));

	// plane-bbox intersection
	Vector3 nf_mask = Vector3((n.x > 0 ? 1 : -1), (n.y > 0 ? 1 : -1), (n.z > 0 ? 1 : -1));
	Vector3 near_corner = multiply(bboxHalfSize, nf_mask);
	Vector3 far_corner = -1.0f * multiply(bboxHalfSize, nf_mask);
	float dist_near_s = sgn(dot(n, near_corner) - dot(n, verts[0]));
	float dist_far_s = sgn(dot(n, far_corner) - dot(n, verts[0]));
	if (fabs(dist_near_s - dist_far_s) < FLT_EPSILON) {
		return false;
	}

	Tri f = { verts[1] - verts[0], verts[2] - verts[1], verts[0] - verts[2] };
	Tri axes = { Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1) };

	for (uint i = 0; i < 3; ++i) {
		for (uint j = 0; j < 3; ++j) {
			Vector3 a = cross(axes[i], f[j]);
			float p0 = a.dot(verts[0]);
			float p1 = a.dot(verts[1]);
			float p2 = a.dot(verts[2]);
			float min = std::fminf(p0, std::fminf(p1, p2));
			float max = std::fmaxf(p0, std::fmaxf(p1, p2));

			float r = bboxHalfSize.dot(Vector3(fabs(a.x), fabs(a.y), fabs(a.z)));
			if (min > r || max < -r) {
				return false;
			}
		}
	}

	return true;
}


AABBTree::AABBTree()
{
}

AABBTree::AABBTree(std::vector<Tri>& triangles, Box3 bbox, uint depthLeft, AABBTree* parent)
{
	this->bbox = bbox;
	this->bbox.min.addScalar(-0.001);
	this->bbox.max.addScalar(0.001);
	this->parent = parent;
	this->depth = MAX_DEPTH - depthLeft;

	this->construct(triangles, depthLeft);
}

AABBTree::~AABBTree()
{
}

bool AABBTree::isLeaf()
{
	return (this->left == nullptr && this->right == nullptr);
}

bool AABBTree::isLeafWithTriangles()
{
	if (!this->isLeaf()) return false;

	return this->triangles.size() > 0;
}

bool AABBTree::hasTriangles()
{
	return this->triangles.size() > 0;
}

void AABBTree::construct(std::vector<Tri>& triangles, uint depthLeft)
{
	if (depthLeft == 0 || triangles.size() <= 2) {
		this->triangles = triangles;
		this->filterTriangles();
		return;
	}

	// choose axis with maximum extent to split by. this appears to be a bit better than to cycle axes
	Vector3 bboxSize = this->bbox.getSize();
	if (bboxSize.x > bboxSize.y && bboxSize.x > bboxSize.z) {
		this->axis = 0;
	} else if (bboxSize.y > bboxSize.z) {
		this->axis = 1;
	} else {
		this->axis = 2;
	}

	std::vector<Tri> leftTriangles = {};
	std::vector<Tri> rightTriangles = {};

	this->splitPosition = this->getSplitPosition(triangles, &leftTriangles, &rightTriangles); // classify triangles

	const bool stopBranching = !this->hasEnoughBranching(leftTriangles.size(), rightTriangles.size(), triangles.size());
	if (stopBranching) {
		this->triangles = triangles;
		this->filterTriangles();
	}
	else {
		if (leftTriangles.size() > 0) {
			Box3 bboxLeft = this->bbox;
			bboxLeft.max.setCoordById(this->splitPosition, this->axis);
			this->left = new AABBTree(leftTriangles, bboxLeft, depthLeft - 1, this);
			if (this->left->triangles.empty() && this->left->left == nullptr && this->left->right == nullptr) {
				this->left = nullptr;
			}
		}
		if (rightTriangles.size() > 0) {
			Box3 bboxRight = this->bbox;
			bboxRight.min.setCoordById(this->splitPosition, this->axis);
			this->right = new AABBTree(rightTriangles, bboxRight, depthLeft - 1, this);
			if (this->right->triangles.empty() && this->right->left == nullptr && this->right->right == nullptr) {
				this->right = nullptr;
			}
		}
	}
}

std::vector<AABBTree> AABBTree::flatten()
{
	std::vector<AABBTree> resultArray = {};

	std::stack<AABBTree> nodeStack;
	nodeStack.push(*this);

	while (nodeStack.size()) {
		AABBTree item = nodeStack.top();
		nodeStack.pop();

		if (!item.left && !item.right) {// is leaf?
			resultArray.push_back(item);
		}
		else {
			if (item.left) {
				nodeStack.push(*item.left);
			}
			if (item.right) {
				nodeStack.push(*item.right);
			}
		}
	}

	return resultArray;
}

std::vector<AABBTree> AABBTree::flattenToDepth(uint depth)
{
	std::vector<AABBTree> resultArray = {};

	std::stack<AABBTree> nodeStack;
	nodeStack.push(*this);

	while (nodeStack.size()) {
		AABBTree item = nodeStack.top();
		nodeStack.pop();

		if (!item.left && !item.right && item.depth == depth) {// is leaf?
			resultArray.push_back(item);
		}
		else {
			if (item.left) {
				nodeStack.push(*item.left);
			}
			if (item.right) {
				nodeStack.push(*item.right);
			}
		}
	}

	return resultArray;
}

std::vector<Geometry> AABBTree::getAABBGeomsOfDepth(uint depth)
{
	std::vector<AABBTree> leafNodes = this->flattenToDepth(depth);
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
	std::vector<AABBTree> leafNodes = this->flatten();
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
	std::vector<AABBTree> leafNodes = this->flattenToDepth(depth);
	std::vector<Geometry> triGeoms = {};

	for (auto&& leaf : leafNodes) {
		std::vector<Geometry> leafTriGeoms = {};

		for (uint i = 0; i < leaf.triangles.size(); i++) {
			Geometry triGeom = Geometry();
			triGeom.uniqueVertices = leaf.triangles[i];

			for (uint j = 0; j < 3; j++) {
				triGeom.vertices.push_back(leaf.triangles[i][j].x);
				triGeom.vertices.push_back(leaf.triangles[i][j].y);
				triGeom.vertices.push_back(leaf.triangles[i][j].z);

				// normals will be computed when merging geoms

				triGeom.vertexIndices.push_back(j);
			}
			triGeoms.push_back(triGeom);
		}
	}

	return triGeoms;
}

float AABBTree::getCostEstimate(float splitPos, uint nLeft, uint nRight)
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

bool AABBTree::hasEnoughBranching(size_t nLeftTris, size_t nRightTris, size_t nTris)
{
	return nLeftTris + nRightTris < 1.5f * nTris;
}

void AABBTree::filterTriangles()
{
	Vector3 bboxCenter = bbox.getCenter();
	Vector3 bboxHalfSize = 0.5f * bbox.getSize();
	std::vector<Tri> newTriangles = {};
	for (uint i = 0; i < this->triangles.size(); i++) {
		// TODO: Check if this evaluates correctly
		if (getTriangleBoundingBoxIntersection(this->triangles[i], bboxCenter, bboxHalfSize)) {
			newTriangles.push_back(this->triangles[i]);
		}
	}
	this->triangles = newTriangles;
}

// TODO: Improve by adaptive sampling
float AABBTree::getSplitPosition(std::vector<Tri>& triangles, std::vector<Tri>* out_left, std::vector<Tri>* out_right)
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
		float min = std::fminf(triangles[i][0].getCoordById(axis), std::fminf(triangles[i][1].getCoordById(axis), triangles[i][2].getCoordById(axis)));
		float max = std::fmaxf(triangles[i][0].getCoordById(axis), std::fmaxf(triangles[i][1].getCoordById(axis), triangles[i][2].getCoordById(axis)));

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
		float min = std::fminf(triangles[i][0].getCoordById(axis), std::fminf(triangles[i][1].getCoordById(axis), triangles[i][2].getCoordById(axis)));
		float max = std::fmaxf(triangles[i][0].getCoordById(axis), std::fmaxf(triangles[i][1].getCoordById(axis), triangles[i][2].getCoordById(axis)));

		if (min <= bestSplitPosition) {
			out_left->push_back(triangles[i]);
		}
		if (max >= bestSplitPosition) {
			out_right->push_back(triangles[i]);
		}
	}

	return bestSplitPosition;
}

uint depth(AABBTree* root)
{
	if (root == nullptr) {
		return 0;
	}

	std::list<AABBTree*> queue;
	queue.push_back(root);

	AABBTree* front = nullptr;
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
