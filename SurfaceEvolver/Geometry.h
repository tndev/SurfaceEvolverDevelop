#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <vector>
#include <algorithm>
#include "Box3.h"
#include "Matrix4.h"
#include "Vector3.h"

using Triangle = std::vector<unsigned int>;
using Triangulation = std::vector<unsigned int>;

class Geometry
{
public:
	// every time a geometry is created, it is preferrable to keep an additional array of vertices to, which
	// vertexIndices actually point
	std::vector<Vector3> uniqueVertices;
	std::vector<float> vertices; // duplicated (for each triangle)
	std::vector<float> normals; // vertex normals (for each triangle, i.e.: there's as many as there are vertex indices)
	// [0, 1, 2, 0, 2, 3, ... ] (e.g.: quads are made of 2 consecutive triplets of vert indices)
	std::vector<unsigned int> vertexIndices; // values correspond to the positions in uniqueVertices array;
	std::vector<Triangulation> triangulations; // each contains ids of triangles inside a polygon

	Geometry();
	~Geometry();
	Geometry(const Geometry& other);

	bool hasNormals();
	bool hasTriangulations();

	Box3 getBoundingBox(Box3 bbox = Box3(), Matrix4 matrix = Matrix4());

	std::vector<unsigned int> getPolygonIndicesFromTriangulation(Triangulation t);	
	std::vector<Vector3> getVertices();

	void flipFaceOrientation();
	void applyMatrix(Matrix4 m);

protected:
	void clear();
private:
	std::vector<unsigned int> getPolygonIndicesFromTriangles(std::vector<Triangle> triangles);
};

#endif