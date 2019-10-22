#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <vector>
#include <set>
#include "Box3.h"
#include "Matrix4.h"
#include "Vector3.h"

class Geometry
{
public:
	std::vector<float> vertices; // duplicated (for each triangle)
	std::vector<float> normals; // vertex normals (for each triangle, i.e.: there's as many as there are vertex indices)
	std::vector<unsigned int> vertexIndices; // [0, 1, 2, 0, 2, 3, ... ] (e.g.: quads are made of 2 consecutive triplets of vert indices)
	std::vector<std::vector<unsigned int>> triangulations; // each contains ids of triangles inside a polygon

	Geometry();
	~Geometry();

	bool hasNormals();
	bool hasTriangulations();
	void copy(Geometry other);
	Geometry clone();

	Box3 getBoundingBox(Box3 bbox = Box3(), Matrix4 matrix = Matrix4());

	std::vector<unsigned int> getPolygonVerticesFromTriangulation(std::vector<std::vector<unsigned int>> triangles);
	std::vector<Vector3> getUniqueVertices();

	void flipFaceOrientation();
	void applyMatrix(Matrix4 m);

protected:
	void clear();
};

#endif