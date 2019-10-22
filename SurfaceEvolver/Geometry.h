#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <vector>
#include "Matrix4.h"
#include "Vector3.h"

class Geometry
{
public:
	std::vector<float> vertices; // unique! no copies based on vertex indices
	std::vector<float> normals; // vertex normals (for each triangle, i.e.: there's as many as there are vertex indices)
	std::vector<unsigned int> vertexIndices; // [0, 1, 2, 0, 2, 3, ... ] (e.g.: quads are made of 2 consecutive triplets of vert indices)
	std::vector<std::vector<unsigned int>> triangulations; // each contains ids of triangles inside a polygon

	Geometry();
	~Geometry();

	bool hasNormals();
	bool hasTriangulations();
	void copy(Geometry other);
	// std::vector<unsigned int> getIndicesForPolygons(std::vector<unsigned int>& indices);
	Geometry clone();

	void flipFaceOrientation();
	void applyMatrix(Matrix4 m);

protected:
	void clear();
};

#endif