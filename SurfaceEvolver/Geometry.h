#ifndef GEOMETRY_H_
#define GEOMETRY_H_

class Geometry
{
public:
	unsigned int nVerts = 0;
	unsigned int nTris = 0;
	float* vertices = nullptr; // [v0x, v0y, v0z, v1x, v1y, v1z, ... ]
	float* normals = nullptr; // vertex normals
	int* vertexIndices = nullptr; // [0, 1, 2, 0, 2, 3, ... ]
	bool quadified = false; // if true then vertexIndices are iterated with an increment of 6 (two triangles), otherwise with 3

	Geometry();
	~Geometry();

	bool hasNormals();
	void copy(Geometry other);
	Geometry clone();
	void clear();
};

#endif