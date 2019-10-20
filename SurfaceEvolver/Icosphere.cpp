#include "Icosphere.h"
#include "Vector3.h"
#include <vector>
#include <array>
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

// ===============================================
// ======== Pre-requisites for Icosphere =========

struct Triangle {
	unsigned int vertex[3];
};

using TriangleList = std::vector<Triangle>;
using VertexList = std::vector<Vector3>;

namespace Icosahedron {
	const float t = (1.0f + sqrt(5.f)) / 2.0f;
	const float norm = sqrt(1.0f + t * t);

	static const VertexList vertices = {
		{-1.0f / norm, t / norm, 0.0f}, {1.0f / norm, t / norm, 0.0f},   {-1.0f / norm, -t / norm,  0.0f},    {1.0f / norm, -t / norm, 0.0f},
		{0.0f, -1.0f / norm, t / norm}, {0.0f, 1.0f / norm, t / norm},    {0.0f, -1.0f / norm, -t / norm},    {0.0f, 1.0f / norm, -t / norm},
		{t / norm, 0.0f, -1.0f / norm}, {t / norm, 0.0f, 1.0f / norm},    {-t / norm, 0.0f, -1.0f / norm},    {-t / norm, 0.0f, 1.0f / norm}
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

TriangleList subdivide(VertexList& vertices, TriangleList triangles) {
	Lookup lookup;
	TriangleList result;

	for (auto&& each:triangles) {
		std::array<unsigned int, 3> mid;
		for (int edge = 0; edge < 3; ++edge) {
			mid[edge] = midpointId(lookup, vertices, each.vertex[edge], each.vertex[(edge + 1) % 3]);
		}

		result.push_back({ each.vertex[0], mid[0], mid[2] });
		result.push_back({ each.vertex[1], mid[1], mid[0] });
		result.push_back({ each.vertex[2], mid[2], mid[1] });
		result.push_back({ mid[0], mid[1], mid[2] });
	}

	return result;
}


void Icosphere::build()
{
	VertexList vertices = Icosahedron::vertices;
	TriangleList triangles = Icosahedron::triangles;

	for (unsigned int i = 0; i < detail; i++) {
		triangles = subdivide(vertices, triangles);
	}

	for (unsigned int i = 0; i < triangles.size(); i++) {
		for (unsigned int j = 0; j < 3; j++) {
			Vector3& v = vertices[triangles[i].vertex[j]];
			this->normals.push_back(v.x);
			this->normals.push_back(v.y);
			this->normals.push_back(v.z);

			this->vertices.push_back(radius * v.x);
			this->vertices.push_back(radius * v.y);
			this->vertices.push_back(radius * v.z + radius);
			this->vertexIndices.push_back(triangles[i].vertex[j]);
		}
	}
}
