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

bool Geometry::hasNormals()
{
	return normals.size();
}

bool Geometry::hasTriangulations()
{
	return triangulations.size();
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

std::vector<unsigned int> Geometry::getPolygonIndicesFromTriangulation(Triangulation t)
{
	Geometry* g = this;
	std::vector<Triangle> triangles = std::vector<Triangle>();
	for (unsigned int k = 0; k < t.size(); k++) {
		triangles.push_back({ g->vertexIndices[3 * t[k]], g->vertexIndices[3 * t[k] + 1], g->vertexIndices[3 * t[k] + 2] });
	}
	std::vector<unsigned int> polygonIds = g->getPolygonIndicesFromTriangles(triangles);
	return polygonIds;
}

std::vector<unsigned int> Geometry::getPolygonIndicesFromTriangles(std::vector<Triangle> triangles)
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

std::vector<Vector3> Geometry::getProjectionsAlongNormal(Face& vertices)
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

std::vector<std::vector<unsigned int>> Geometry::getTriangulatedIndices(Face& vertices)
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

Vector3 Geometry::getNormal(Face f)
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
}
