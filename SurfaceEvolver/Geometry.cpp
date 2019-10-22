#include "Geometry.h"

Geometry::Geometry()
{
}

Geometry::~Geometry()
{
	clear();
}

bool Geometry::hasNormals()
{
	return normals.size();
}

bool Geometry::hasTriangulations()
{
	return triangulations.size();
}

void Geometry::copy(Geometry other)
{
	if (other.hasNormals()) {
		normals = std::vector<float>(other.normals);
	}

	vertices = std::vector<float>(other.vertices);
	vertexIndices = std::vector<unsigned int>(other.vertexIndices);

	if (other.hasTriangulations()) {
		triangulations = std::vector<std::vector<unsigned int>>(other.triangulations);
	}
}

Geometry Geometry::clone()
{
	Geometry result = Geometry();
	result.copy(*this);
	return result;
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

std::vector<unsigned int> Geometry::getPolygonVerticesFromTriangulation(std::vector<std::vector<unsigned int>> triangles)
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

std::vector<Vector3> Geometry::getUniqueVertices()
{
	std::vector<Vector3> vertices = std::vector<Vector3>();
	std::set<Vector3> vertexSet = std::set<Vector3>();

	for (unsigned int i = 0; i < this->vertices.size(); i += 9) {

		for (unsigned int j = 0; j < 3; j++) {
			unsigned int vId = i + 3 * j;
			Vector3 v = Vector3(this->vertices[vId], this->vertices[vId + 1], this->vertices[vId + 2]);
			std::pair<std::set<Vector3>::iterator, bool> unique = vertexSet.insert(v);

			if (unique.second) {
				vertices.push_back(v);
			}
		}
	}
	return vertices;
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

void Geometry::clear()
{
	vertices.clear();
	normals.clear();
	vertexIndices.clear();
	triangulations.clear();
}
