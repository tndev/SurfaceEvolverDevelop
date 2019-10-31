#include "Geometry.h"

#define failedTriang 0

Geometry::Geometry()
{
}

Geometry::~Geometry()
{
	clear();
}

Geometry::Geometry(const Geometry& other)
{
	normals = std::vector<float>(other.normals);
	uniqueVertices = std::vector<Vector3>(other.uniqueVertices);
	vertices = std::vector<float>(other.vertices);
	vertexIndices = std::vector<unsigned int>(other.vertexIndices);
	triangulations = std::vector<std::vector<unsigned int>>(other.triangulations);
}

bool Geometry::hasVertices()
{
	return this->vertices.size() > 0;
}

bool Geometry::hasVertexIndices()
{
	return this->vertexIndices.size() > 0;
}

bool Geometry::hasNormals()
{
	return this->normals.size() > 0;
}

bool Geometry::hasTriangulations()
{
	return this->triangulations.size() > 0;
}

Box3 Geometry::getBoundingBox(Box3 bbox, Matrix4 matrix)
{
	Vector3 helperVector = Vector3();
	for (unsigned int i = 0; i < vertices.size(); i += 3) {
		helperVector.set(vertices[i], vertices[i + 1], vertices[i + 2]);
		helperVector.applyMatrix4(matrix);
		bbox.expandByPoint(helperVector);
	}
	return bbox;
}

void Geometry::computeNormals()
{
	StructGeom::Triangle faceVerts = { Vector3(), Vector3(), Vector3() };

	Vector3 normal = Vector3();
	for (unsigned int i = 0; i < this->vertices.size(); i += 9) {

		for (unsigned int j = 0; j < 3; j++) {
			faceVerts[j].x = this->vertices[i + j * 3];
			faceVerts[j].y = this->vertices[i + j * 3 + 1];
			faceVerts[j].z = this->vertices[i + j * 3 + 2];
		}

		normal = getTriangleNormal(faceVerts, normal);

		for (unsigned int j = 0; j < 3; j++) {
			this->normals.push_back(normal.x);
			this->normals.push_back(normal.y);
			this->normals.push_back(normal.z);
		}
	}
}

void Geometry::computeTriangulations()
{
	if (this->hasVertexIndices() || this->hasVertices()) {
		unsigned int N = (this->hasVertexIndices() && !this->hasVertices()) ? this->vertexIndices.size() : this->vertices.size() / 3;
		for (unsigned int i = 0; i < N; i += 3) {
			this->triangulations.push_back({ i / 3 });
		}
	}
	// otherwise the geometry is just a point cloud and has to be triangulated separately
}

std::vector<unsigned int> Geometry::getPolygonIndicesFromTriangulation(BufferGeom::Triangulation t)
{
	Geometry* g = this;
	std::vector<BufferGeom::Triangle> triangles = std::vector<BufferGeom::Triangle>();
	for (unsigned int k = 0; k < t.size(); k++) {
		triangles.push_back({ g->vertexIndices[3 * t[k]], g->vertexIndices[3 * t[k] + 1], g->vertexIndices[3 * t[k] + 2] });
	}
	std::vector<unsigned int> polygonIds = g->getPolygonIndicesFromTriangles(triangles);
	return polygonIds;
}

std::vector<unsigned int> Geometry::getPolygonIndicesFromTriangles(std::vector<BufferGeom::Triangle> triangles)
{
	if (triangles.size() > 1) {
		std::vector<std::vector<unsigned int>> edges = std::vector<std::vector<unsigned int>>();

		for (auto&& t : triangles) {
			edges.push_back({ t[0], t[1] });
			edges.push_back({ t[1], t[2] });
			edges.push_back({ t[2], t[0] });
		}
		for (unsigned int i = 0; i < edges.size() - 1; i++) {
			for (unsigned int j = i + 1; j < edges.size(); j++) {
				if ((edges[i][0] == edges[j][0] && edges[i][1] == edges[j][1]) || (edges[i][0] == edges[j][1] && edges[i][1] == edges[j][0])) {
					edges.erase(edges.begin() + j);
					edges.erase(edges.begin() + i);
					i--;
					break;
				}
			}
		}

		std::vector<unsigned int> faceVertIndices = std::vector<unsigned int>();
		faceVertIndices.push_back(edges[0][1]);

		std::vector<std::vector<unsigned int>> edgesSet = std::vector<std::vector<unsigned int>>(edges);
		edgesSet.erase(edgesSet.begin());

		for (unsigned int i = 1; i < edges.size(); i++) {
			for (unsigned int j = 0; j < edgesSet.size(); j++) {
				if (edgesSet[j][0] == faceVertIndices[faceVertIndices.size() - 1]) {
					faceVertIndices.push_back(edgesSet[j][1]);
					edgesSet.erase(edgesSet.begin() + j);
					break;
				}
			}
		}

		return faceVertIndices;
	}
	else if (triangles.size() == 1) {
		return triangles[0];
	}

	return std::vector<unsigned int>();
}

std::vector<Vector3> Geometry::getVertices()
{
	std::vector<Vector3> result = std::vector<Vector3>();

	for (unsigned int i = 0; i < this->vertices.size(); i += 3) {
		result.push_back(Vector3(this->vertices[i], this->vertices[i + 1], this->vertices[i + 2]));
	}
	return result;
}

std::vector<Vector3> Geometry::getUniqueVertices()
{
	return this->uniqueVertices;
}

std::vector<Vector3> Geometry::getProjectionsAlongNormal(BufferGeom::Face& vertices)
{
	Vector3 normal = getNormal(vertices); // directions
	Vector3 referencePoint = vertices[0];
	Vector3 up = Vector3(0, 0, 1);

	// ------ -----------
	// TODO: create Quaternion

	if (fabs(dot(normal, up)) >= 0.999) {
		if (fabs(normal.x) < fabs(normal.y) && fabs(normal.x) < fabs(normal.z)) {
			up = normalize(Vector3(1, 0, 0).cross(normal));
		}
		else if (fabs(normal.y) < fabs(normal.z)) {
			up = normalize(Vector3(0, 1, 0).cross(normal));
		}
		else {
			up = normalize(Vector3(0, 0, 1).cross(normal));
		}
	}

	Vector3 crossProduct = normalize(cross(up, normal));
	Vector3 pseudoUp = normalize(cross(normal, crossProduct));

	Matrix4 transform = Matrix4(
		normal.x, normal.y, normal.z, 0,
		crossProduct.x, crossProduct.y, crossProduct.z, 0,
		pseudoUp.x, pseudoUp.y, pseudoUp.z, 0,
		0, 0, 0, 1
	);
	transform.transpose();

	Vector3 projX = transform * Vector3(0, 1, 0);
	Vector3 projY = transform * Vector3(0, 0, 1);

	std::vector<Vector3> projections = std::vector<Vector3>();
	for (unsigned int k = 0; k < vertices.size(); k++) {
		Vector3 vec = vertices[k] - referencePoint;
		float px = dot(vec, projX);
		float py = dot(vec, projY);
		projections.push_back(Vector3(px, py, 0.0f));
	}

	return projections;
}

std::vector<std::vector<unsigned int>> Geometry::getTriangulatedIndices(BufferGeom::Face& vertices)
{
	std::vector<std::vector<unsigned int>> faces = std::vector<std::vector<unsigned int>>();
	if (vertices.size() < 3) {
	}
	else if (vertices.size() == 3) {
		faces = { {0, 1, 2} };
	}
	else if (vertices.size() == 4) {
		faces = { {0, 1, 2}, {0, 2, 3} };
		Vector3 e2 = vertices[2] - vertices[0];
		Vector3 e1 = vertices[1] - vertices[0];
		Vector3 e3 = vertices[3] - vertices[0];
		Vector3 c21 = cross(e2, e1);
		Vector3 c23 = cross(e2, e3);
		if (dot(c21, c23) > 0.0) {
			faces = { {0, 1, 3}, {1, 2, 3} };
		}
	}

	else {

		std::vector<Vector3> projections = getProjectionsAlongNormal(vertices);

		p2t::CDT* cdt = NULL;

		// poly2tri doesn't work well with duplicate or collinear points. it is possible that it will fail on certain inputs
		// if it fails, we retry with increasingly higher jitter. it usually works on first retry
		// if it doesn't work with 42 retries, there is no hope
		std::vector<p2t::Point> points = {};

		unsigned int MAX_TRIES = 42;
		unsigned int tries_utilized = 0;
		srand((unsigned)time(NULL));
		for (tries_utilized = 0; tries_utilized < MAX_TRIES; tries_utilized++) {
			//// poly2tri triangulation, seems to work well with holes
			points.clear();
			for (unsigned int i = 0; i < projections.size(); i++) {
				p2t::Point pt = p2t::Point(
					projections[i].x + (tries_utilized * (rand() / RAND_MAX * 0.001 - 0.0005)),
					projections[i].y + (tries_utilized * (rand() / RAND_MAX * 0.001 - 0.0005))
				);
				pt._vIdx = i;
				points.push_back(pt);
			}
			
			try {
				std::vector<p2t::Point*> contour = {};
				for (unsigned int i = 0; i < points.size(); i++) {
					contour.push_back(&points[i]);
				}
				cdt = new p2t::CDT(contour);
				cdt->Triangulate();
				if (cdt == NULL) {
					throw(failedTriang);
				}
			}
			catch (int e) {
				continue;
			}
			break;// assuming it all goes well
		}
		if (cdt == NULL) {
			std::cout << "ERROR! Triangulation was not able to finish with" << MAX_TRIES << "tries" << std::endl;
			return {};
		}

		std::vector<p2t::Triangle*> triangles = cdt->GetTriangles();

		for (unsigned int i = 0; i < triangles.size(); i++) {
			faces.push_back({
				triangles[i]->GetPoint(0)->_vIdx,
				triangles[i]->GetPoint(1)->_vIdx,
				triangles[i]->GetPoint(2)->_vIdx
			});
		}
	}

	return faces;
}

template<class T>
void flipArray(std::vector<T> arr, unsigned int elementSize) {
	if (arr.size() && elementSize) {
		unsigned int iStep = 3 * elementSize;
		float tmp[3] = { 0.0f, 0.0f, 0.0f };
		for (unsigned int i = 0; i < arr.size(); i += iStep) {
			for (unsigned int j = 0; j < elementSize; j++) {
				const unsigned int i0 = i + j;
				const unsigned int i1 = i + j + 2 * elementSize;
				tmp[j] = arr[i0];
				arr[i0] = arr[i1];
				arr[i1] = tmp[j];
			}
		}
	}	
}

void Geometry::flipFaceOrientation()
{
	flipArray(vertices, 3);
	if (hasNormals()) {
		flipArray(normals, 3);
	}
	flipArray(vertexIndices, 1);
}

Matrix3 getSubMatrix3(Matrix4 m) {
	float* e = m.elements;
	float resultElems[9] = {
		e[0],	e[4],	e[8],
		e[1],	e[5],	e[9],
		e[2],	e[6],	e[10]
	};
	return Matrix3(resultElems);
}

void Geometry::applyMatrix(Matrix4 m)
{
	if (m.isIdentity()) {
		return;
	}

	const bool hasNormals = this->hasNormals();
	const bool changeOrientation = (m.determinant() < 0);
	Matrix4 mInverse = Matrix4().getInverse(m);
	Matrix4 tmp = transpose(mInverse);
	Matrix3 matrix3 = m.getSubMatrix3();
	Matrix3 matrix3InverseTransposed = getSubMatrix3(tmp);

	Vector3 helperVector = Vector3();

	for (int i = 0; i < vertices.size(); i += 3) {
		const unsigned int i0 = i, i1 = i + 1, i2 = i + 2;

		helperVector.set(vertices[i0], vertices[i1], vertices[i2]);
		helperVector.applyMatrix4(m);

		vertices[i0] = helperVector.x;
		vertices[i1] = helperVector.y;
		vertices[i2] = helperVector.z;

		this->uniqueVertices[this->vertexIndices[i / 3]] = helperVector;

		if (hasNormals) {
			helperVector.set(normals[i0], normals[i1], normals[i2]);
			helperVector.applyMatrix3(matrix3InverseTransposed);
			const float length = helperVector.length();

			normals[i0] = helperVector.x / length;
			normals[i1] = helperVector.y / length;
			normals[i2] = helperVector.z / length;
		}
	}

	if (changeOrientation) {
		flipFaceOrientation();
	}
}

Vector3 Geometry::getNormal(BufferGeom::Face f)
{
	Vector3 normal = Vector3();
	for (unsigned int i = 0; i < f.size(); i++) {
		Vector3 fromNext = f[i] - f[(i + 1) % f.size()];
		Vector3 plusNext = f[i] + f[(i + 1) % f.size()];
		normal.x = normal.x + fromNext.y * plusNext.z;
		normal.y = normal.y + fromNext.z * plusNext.x;
		normal.z = normal.z + fromNext.x * plusNext.y;
	}
	return normalize(normal);
}

void Geometry::clear()
{
	uniqueVertices.clear();
	vertices.clear();
	normals.clear();
	vertexIndices.clear();
	triangulations.clear();

	triangles.clear();
	edges.clear();
}

std::pair<std::vector<BufferGeom::Triangulation>, std::vector<size_t>> Geometry::getSortedPolygonTriangulationsAndSizes()
{
	std::vector<BufferGeom::Triangulation> T = triangulations;
	std::vector<size_t> sizes = std::vector<size_t>();
	std::sort(T.begin(), T.end(), [](const BufferGeom::Triangulation& a, const BufferGeom::Triangulation& b) { return a.size() < b.size(); });
	unsigned int tCount = 0; // triangulation count
	size_t currentSize = T.begin()->size();

	for (std::vector<BufferGeom::Triangulation>::iterator it = T.begin(); it < T.end(); ++it) {
		if (it->size() > currentSize) {
			currentSize = it->size();
			sizes.push_back(tCount);
			tCount = 0;
			continue;
		}
		tCount++;
	}
	sizes.push_back(tCount);

	std::pair<std::vector<BufferGeom::Triangulation>, std::vector<size_t>> result = { T, sizes };
	return result;
}

std::vector<StructGeom::Triangle> Geometry::getTriangles()
{
	if (this->triangles.empty()) {
		for (unsigned int i = 0; i < this->vertexIndices.size(); i += 3) {
			Vector3 v0 = this->uniqueVertices[this->vertexIndices[i]];
			Vector3 v1 = this->uniqueVertices[this->vertexIndices[i + 1]];
			Vector3 v2 = this->uniqueVertices[this->vertexIndices[i + 2]];

			StructGeom::Triangle T = { v0, v1, v2 };

			this->triangles.push_back(T);
		}
	}

	return this->triangles;
}

// TODO: Use sets to get unique edges
std::vector<StructGeom::Edge> Geometry::getEdges()
{
	if (this->edges.empty()) {
		for (unsigned int i = 0; i < this->vertexIndices.size(); i += 3) {
			Vector3 v0 = this->uniqueVertices[this->vertexIndices[i]];
			Vector3 v1 = this->uniqueVertices[this->vertexIndices[i + 1]];
			Vector3 v2 = this->uniqueVertices[this->vertexIndices[i + 2]];

			StructGeom::Edge e0 = { v0, v1 };
			StructGeom::Edge e1 = { v1, v2 };
			StructGeom::Edge e2 = { v2, v0 };

			this->edges.push_back(e0);
			this->edges.push_back(e1);
			this->edges.push_back(e2);
		}
	}

	return this->edges;
}

Geometry mergeGeometries(std::vector<Geometry>& geometries)
{
	Geometry result = Geometry();
	unsigned int vertexCount = 0;
	unsigned int uniqueVertCount = 0;

	for (unsigned int i = 0; i < geometries.size(); i++) {
		if (!geometries[i].hasNormals()) {
			geometries[i].computeNormals();
		}
		if (!geometries[i].hasTriangulations()) {
			geometries[i].computeTriangulations();
		}
		vertexCount += geometries[i].vertices.size() / 3;
		uniqueVertCount += geometries[i].uniqueVertices.size();
	};

	result.uniqueVertices = {};
	result.vertices = {};
	result.normals = {};
	result.vertexIndices = std::vector<unsigned int>(vertexCount);
	result.triangulations = {};

	unsigned int idx = 0;
	unsigned int uIdx = 0;
	for (unsigned int i = 0; i < geometries.size(); i++) {
		const unsigned int currVertexCount = geometries[i].vertices.size() / 3;
		const unsigned int currUniqueCount = geometries[i].uniqueVertices.size();

		result.uniqueVertices.insert(result.uniqueVertices.end(), geometries[i].uniqueVertices.begin(), geometries[i].uniqueVertices.end());
		result.vertices.insert(result.vertices.end(), geometries[i].vertices.begin(), geometries[i].vertices.end());
		result.normals.insert(result.normals.end(), geometries[i].normals.begin(), geometries[i].normals.end());

		for (unsigned int j = 0; j < currVertexCount; j++) {
			const unsigned int vertexIdx = geometries[i].hasVertexIndices() ? geometries[i].vertexIndices[j] : j;
			result.vertexIndices[j + idx] = vertexIdx + uIdx;
		}

		if (geometries[i].hasTriangulations()) {
			for (unsigned int j = 0; j < geometries[i].triangulations.size(); j++) {
				BufferGeom::Triangulation triangulationCopy = geometries[i].triangulations[j];
				for (unsigned int k = 0; k < triangulationCopy.size(); k++) {
					triangulationCopy[k] += idx / 3;
				}
				result.triangulations.push_back(triangulationCopy);
			}

		}

		idx += currVertexCount;
		uIdx += currUniqueCount;
	}
	return result;
}

Vector3 getTriangleNormal(StructGeom::Triangle triangle, Vector3 resultNormal)
{
	Vector3 u = triangle[1] - triangle[0];
	Vector3 v = triangle[2] - triangle[0];

	resultNormal.set(
		(u.y * v.z) - (u.z * v.y),
		(u.z * v.x) - (u.x * v.z),
		(u.x * v.y) - (u.y * v.x)
	);
	return resultNormal;
}

using Tri = StructGeom::Triangle;
using uint = unsigned int;
template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

bool getTriangleBoundingBoxIntersection(Tri& vertices, Vector3& bboxCenter, Vector3& bboxHalfSize, Vector3* optTriNormal)
{
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
