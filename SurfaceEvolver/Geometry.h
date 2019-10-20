#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <vector>
#include "Matrix4.h"
#include "Vector3.h"

class Geometry
{
public:
	std::vector<float> vertices; 
	std::vector<float> normals; // vertex normals
	std::vector<unsigned int> vertexIndices; // [0, 1, 2, 0, 2, 3, ... ]
	std::vector<std::vector<unsigned int>> triangulations;

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