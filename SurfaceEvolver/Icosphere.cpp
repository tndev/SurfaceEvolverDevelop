#include "Icosphere.h"
#include "Vector3.h"
#include <vector>
#include <map>

Icosphere::Icosphere()
{
}

Icosphere::Icosphere(unsigned int detail, float radius)
{
	this->detail = detail; this->radius = radius;
	build();
}

Icosphere::~Icosphere()
{
}

void Icosphere::copy(Icosphere other)
{
	Geometry::copy(other);
	detail = other.detail;
	radius = other.radius;
}

Icosphere Icosphere::clone()
{
	Icosphere* result = new Icosphere();
	result->copy(*this);
	return *result;
}

struct Triangle {
	unsigned int vertex[3];
};

using TriangleList = std::vector<Triangle>;
using VertexList = std::vector<Vector3>;

namespace Icosahedron {
	const float t = (1.f + sqrt(5.f)) / 2.f;
	const float norm = sqrt(1.f + t * t);

	static const VertexList vertices = {
		   {-1.f / norm, t / norm, 0.f}, {1.f / norm, t / norm, 0.f},   {-1.f / norm, -t / norm,  0.f},    {1.f / norm, -t / norm, 0.f},
		   {0.f, -1.f / norm, t / norm}, {0.f, 1.f / norm, t / norm},    {0.f, -1.f / norm, -t / norm},    {0.f, 1.f / norm, -t / norm},
		   {t / norm, 0.f, -1.f / norm}, {t / norm, 0.f, 1.f / norm},    {-t / norm, 0.f, -1.f / norm},    {-t / norm, 0.f, 1.f / norm}
	};

	static const TriangleList triangles = {
		{0, 11, 5},    {0, 5, 1},     {0, 1, 7},    {0, 7, 10},    {0, 10, 11},
		{1, 5, 9},    {5, 11, 4},   {11, 10, 2},    {10, 7, 6},      {7, 1, 8},
		{3, 9, 4},     {3, 4, 2},     {3, 2, 6},     {3, 6, 8},      {3, 8, 9},
		{4, 9, 5},    {2, 4, 11},    {6, 2, 10},     {8, 6, 7},      {9, 8, 1}
	};
};


using Lookup = std::map<std::pair<unsigned int, unsigned int>, unsigned int>;

unsigned int midpointId(Lookup& lookup, VertexList& vertices, unsigned int first, unsigned int second)
{
	Lookup::key_type key(first, second);

	if (key.first > key.second) {
		std::swap(key.first, key.second);
	}
	auto inserted = lookup.insert({ key, vertices.size() });
	if (inserted.second) {
		auto& edge0 = vertices[first];
		auto& edge1 = vertices[second];
		auto midPt = normalize(edge0 + edge1);
		vertices.push_back(midPt);
	}

	return inserted.first->second;
};

// using IndexedMesh = std::pair<VertexList, TriangleList>;


void Icosphere::build()
{
	VertexList vertices = Icosahedron::vertices;
	TriangleList triangles = Icosahedron::triangles;

	for (unsigned int i = 0; i < detail; i++) {

	}


}
