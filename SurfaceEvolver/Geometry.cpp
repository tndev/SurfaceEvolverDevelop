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
	return (normals != nullptr);
}

void Geometry::copy(Geometry other)
{
	nVerts = other.nVerts; nTris = other.nTris; quadified = other.quadified;
	vertices = new float[3 * nVerts];
	if (other.hasNormals()) {
		normals = new float[3 * nVerts];

		for (int i = 0; i < 3 * nVerts; i++) {
			normals[i] = other.normals[i];
		}
	}

	for (int i = 0; i < 3 * nVerts; i++) {
		vertices[i] = other.vertices[i];		
	}

	for (int i = 0; i < 3 * nTris; i++) {
		vertexIndices[i] = other.vertexIndices[i];
	}
}

Geometry Geometry::clone()
{
	Geometry *result = new Geometry();
	result->copy(*this);
	return *result;
}

void Geometry::clear()
{
	delete vertices;
	delete normals;
	delete vertexIndices;
	nVerts = 0; nTris = 0; quadified = false;
}
