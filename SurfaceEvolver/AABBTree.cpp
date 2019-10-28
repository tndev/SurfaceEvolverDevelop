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

	Vector3 axis_a0_f0 = cross(axes[0], f[0]);
	Vector3 axis_a0_f1 = cross(axes[0], f[1]);
	Vector3 axis_a0_f2 = cross(axes[0], f[2]);

	Vector3 axis_a1_f0 = cross(axes[1], f[0]);
	Vector3 axis_a1_f1 = cross(axes[1], f[1]);
	Vector3 axis_a1_f2 = cross(axes[1], f[2]);

	Vector3 axis_a2_f0 = cross(axes[2], f[0]);
	Vector3 axis_a2_f1 = cross(axes[2], f[1]);
	Vector3 axis_a2_f2 = cross(axes[2], f[2]);

	auto testSeparationWithAxis = [&](Vector3& axis) -> bool {
		float p0 = dot(verts[0], axis);
		float p1 = dot(verts[1], axis);
		float p2 = dot(verts[2], axis);

		float r =
			bboxHalfSize.x * fabs(dot(axes[0], axis)) +
			bboxHalfSize.y * fabs(dot(axes[1], axis)) +
			bboxHalfSize.z * fabs(dot(axes[2], axis));

		if (std::fmaxf(-std::max({ p0, p1, p2 }), std::min({ p0, p1, p2 })) > r) {
			return true;
		}

		return false;
	};

	// test 0:

	if (testSeparationWithAxis(axis_a0_f0)) {
		return false;
	}

	// test 1:

	if (testSeparationWithAxis(axis_a0_f1)) {
		return false;
	}

	// test 2:

	if (testSeparationWithAxis(axis_a0_f2)) {
		return false;
	}

	// -----------

	// test 3:

	if (testSeparationWithAxis(axis_a1_f0)) {
		return false;
	}

	// test 4:

	if (testSeparationWithAxis(axis_a1_f1)) {
		return false;
	}

	// test 5:

	if (testSeparationWithAxis(axis_a1_f2)) {
		return false;
	}

	// -----------

	// test 6:

	if (testSeparationWithAxis(axis_a2_f0)) {
		return false;
	}

	// test 7:

	if (testSeparationWithAxis(axis_a2_f1)) {
		return false;
	}

	// test 8:

	if (testSeparationWithAxis(axis_a2_f2)) {
		return false;
	}

	// box normal tests:

	if (testSeparationWithAxis(axes[0])) {
		return false;
	}

	if (testSeparationWithAxis(axes[1])) {
		return false;
	}

	if (testSeparationWithAxis(axes[2])) {
		return false;
	}

	// triangle normal test:

	if (testSeparationWithAxis(n)) {
		return false;
	}

	// passed all 13 tests

	return true;
}


AABBTree::AABBTree()
{
}

AABBTree::AABBTree(std::vector<Tri>& triangles, Box3 bbox, uint depthLeft)
{
	this->bbox = bbox;
	this->bbox.min.addScalar(-0.001);
	this->bbox.max.addScalar(0.001);
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
		if (leftTriangles.size()) {
			Box3 bboxLeft = this->bbox;
			bboxLeft.max.setCoordById(this->splitPosition, this->axis);
			this->left = new AABBTree(leftTriangles, bboxLeft, depthLeft - 1);
			if (!this->left->triangles.size() && this->left->left == nullptr && this->left->right == nullptr) {
				this->left = nullptr;
			}
		}
		if (rightTriangles.size()) {
			Box3 bboxRight = this->bbox;
			bboxRight.min.setCoordById(this->splitPosition, this->axis);
			this->right = new AABBTree(rightTriangles, bboxRight, depthLeft - 1);
			if (!this->right->triangles.size() && this->right->left == nullptr && this->right->right == nullptr) {
				this->right = nullptr;
			}
		}
	}
}

// crashes on depth >= 1, this->depth = NULL!
std::vector<Geometry> AABBTree::getAABBGeomsOfDepth(uint depth)
{
	if (this->depth < depth) {
		std::vector<Geometry> boxesLeft = {};
		std::vector<Geometry> boxesRight = {};
		if (left != nullptr) {
			boxesLeft = left->getAABBGeomsOfDepth(depth - 1);
		}
		if (right != nullptr) {
			boxesRight = right->getAABBGeomsOfDepth(depth - 1);
		}
		boxesLeft.insert(boxesLeft.end(), boxesRight.begin(), boxesRight.end());

		return boxesLeft;
	}

	float dimX = bbox.max.x - bbox.min.x;
	float dimY = bbox.max.y - bbox.min.y;
	float dimZ = bbox.max.z - bbox.min.z;
	PrimitiveBox box = PrimitiveBox(dimX, dimY, dimZ, 1, 1, 1);
	Vector3 t = bbox.min;
	box.applyMatrix(Matrix4().makeTranslation(t.x, t.y, t.z));

	return { box };
}

std::vector<Geometry> AABBTree::getAABBLeafGeoms()
{
	if (!this->isLeaf()) {
		std::vector<Geometry> boxesLeft = {};
		std::vector<Geometry> boxesRight = {};
		if (left != nullptr) {
			boxesLeft = left->getAABBLeafGeoms();
		}
		if (right != nullptr) {
			boxesRight = right->getAABBLeafGeoms();
		}
		boxesLeft.insert(boxesLeft.end(), boxesRight.begin(), boxesRight.end());

		return boxesLeft;
	}

	float dimX = bbox.max.x - bbox.min.x;
	float dimY = bbox.max.y - bbox.min.y;
	float dimZ = bbox.max.z - bbox.min.z;
	PrimitiveBox box = PrimitiveBox(dimX, dimY, dimZ, 1, 1, 1);
	Vector3 t = bbox.min;
	box.applyMatrix(Matrix4().makeTranslation(t.x, t.y, t.z));

	return { box };
}

float AABBTree::getCostEstimate(float splitPos, uint nLeft, uint nRight)
{
	Box3 l_bbox = this->bbox;
	l_bbox.max.setCoordById(this->axis, splitPos);
	Box3 r_bbox = this->bbox;
	r_bbox.min.setCoordById(this->axis, splitPos);
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
