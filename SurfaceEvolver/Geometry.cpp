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

void Geometry::copy(Geometry other)
{
	quadified = other.quadified;
	if (other.hasNormals()) {
		for (unsigned int i = 0; i < other.normals.size(); i++) {
			normals.push_back(other.normals[i]);
		}
	}

	for (unsigned int i = 0; i < other.vertices.size(); i++) {
		vertices.push_back(other.vertices[i]);		
	}

	for (unsigned int i = 0; i < other.vertexIndices.size(); i++) {
		vertexIndices.push_back(other.vertexIndices[i]);
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
	vertices.clear();
	normals.clear();
	vertexIndices.clear();
	quadified = false;
}
