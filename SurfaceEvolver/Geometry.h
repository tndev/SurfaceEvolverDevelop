#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <vector>
#include "Matrix4.h"

class Geometry
{
public:
	std::vector<float> vertices; 
	std::vector<float> normals; // vertex normals
	std::vector<unsigned int> vertexIndices; // [0, 1, 2, 0, 2, 3, ... ]
	bool quadified = false; // if true then vertexIndices are iterated with an increment of 6 (two triangles), otherwise with 3

	Geometry();
	~Geometry();

	bool hasNormals();
	void copy(Geometry other);
	Geometry clone();

	void applyMatrix(Matrix4 m);
protected:
	void clear();
};

#endif