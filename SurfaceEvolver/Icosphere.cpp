#include "Icosphere.h"
#include "Vector3.h"
#include <vector>
#include <array>
#include <map>

IcoSphere::IcoSphere()
{
}

IcoSphere::IcoSphere(unsigned int detail, float radius)
{
	this->detail = detail; this->radius = radius;
	build();
}

IcoSphere::~IcoSphere()
{
}

void IcoSphere::copy(IcoSphere other)
{
	Geometry::copy(other);
	detail = other.detail;
	radius = other.radius;
}

IcoSphere IcoSphere::clone()
{
	IcoSphere result = IcoSphere();
	result.copy(*this);
	return result;
}

// ===============================================
// ======== Pre-requisites for Icosphere =========

struct Triangle {
	unsigned int vertex[3];
};

using TriangleList = std::vector<Triangle>;
using VertexList = std::vector<Vector3>;

namespace Icosahedron {
	const float t = (1.0f + sqrt(5.0f)) / 2.0f;
	const float norm = sqrt(1.0f + t * t);

	static const VertexList vertices = {
		{-1.0f / norm, t / norm, 0.0f},    {1.0f / norm, t / norm, 0.0f},   {-1.0f / norm, -t / norm,  0.0f},    {1.0f / norm, -t / norm, 0.0f},
		{0.0f, -1.0f / norm, t / norm},    {0.0f, 1.0f / norm, t / norm},    {0.0f, -1.0f / norm, -t / norm},    {0.0f, 1.0f / norm, -t / norm},
		{t / norm, 0.0f, -1.0f / norm},    {t / norm, 0.0f, 1.0f / norm},    {-t / norm, 0.0f, -1.0f / norm},    {-t / norm, 0.0f, 1.0f / norm}
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

TriangleList subdivide(VertexList& vertices, TriangleList& triangles) {
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


void IcoSphere::build()
{
	VertexList vertices = Icosahedron::vertices;
	TriangleList triangles = Icosahedron::triangles;

	for (unsigned int i = 0; i < detail; i++) {
		triangles = subdivide(vertices, triangles);
	}

	for (unsigned int i = 0; i < vertices.size(); i++) {
		this->vertices.push_back(vertices[i].x);
		this->vertices.push_back(vertices[i].y);
		this->vertices.push_back(vertices[i].z);

		this->uniqueVertices.push_back(radius * vertices[i]);
	}

	for (unsigned int i = 0; i < triangles.size(); i++) {
		for (unsigned int j = 0; j < 3; j++) {
			this->vertexIndices.push_back(triangles[i].vertex[j]);
		}
	}

	std::vector<float> geometryVertices = std::vector<float>(3 * this->vertexIndices.size());
	this->normals = std::vector<float>(3 * this->vertexIndices.size());


	// apply radius to duplicate vertices
	for (unsigned int i = 0; i < this->vertexIndices.size(); i++) {
		geometryVertices[i * 3] = radius * this->vertices[this->vertexIndices[i] * 3];
		geometryVertices[i * 3 + 1] = radius * this->vertices[this->vertexIndices[i] * 3 + 1];
		geometryVertices[i * 3 + 2] = radius * this->vertices[this->vertexIndices[i] * 3 + 2];
	}

	unsigned int triId = 0;

	auto bufferTriangleFromIds = [&](
		unsigned int i0, unsigned int i1, unsigned int i2
	) {
			this->vertexIndices[3 * triId] = i0;
			this->vertexIndices[3 * triId + 1] = i1;
			this->vertexIndices[3 * triId + 2] = i2;

			if (detail > 0) {
				this->normals[9 * triId] = this->vertices[3 * i0];
				this->normals[9 * triId + 1] = this->vertices[3 * i0 + 1];
				this->normals[9 * triId + 2] = this->vertices[3 * i0 + 2];

				this->normals[9 * triId + 3] = this->vertices[3 * i1];
				this->normals[9 * triId + 4] = this->vertices[3 * i1 + 1];
				this->normals[9 * triId + 5] = this->vertices[3 * i1 + 2];

				this->normals[9 * triId + 6] = this->vertices[3 * i2];
				this->normals[9 * triId + 7] = this->vertices[3 * i2 + 1];
				this->normals[9 * triId + 8] = this->vertices[3 * i2 + 2];
			}
			else {
				Vector3 pA = Vector3(this->vertices[3 * i0], this->vertices[3 * i0 + 1], this->vertices[3 * i0 + 2]);
				Vector3 pB = Vector3(this->vertices[3 * i1], this->vertices[3 * i1 + 1], this->vertices[3 * i1 + 2]);
				Vector3 pC = Vector3(this->vertices[3 * i2], this->vertices[3 * i2 + 1], this->vertices[3 * i2 + 2]);
				Vector3 centroid = Vector3();

				centroid = (pA + pB + pC) / 3.0f;

				float norm = centroid.length();

				this->normals[9 * triId] = centroid.x / norm;
				this->normals[9 * triId + 1] = centroid.y / norm;
				this->normals[9 * triId + 2] = centroid.z / norm;

				this->normals[9 * triId + 3] = centroid.x / norm;
				this->normals[9 * triId + 4] = centroid.y / norm;
				this->normals[9 * triId + 5] = centroid.z / norm;

				this->normals[9 * triId + 6] = centroid.x / norm;
				this->normals[9 * triId + 7] = centroid.y / norm;
				this->normals[9 * triId + 8] = centroid.z / norm;
			}

			this->triangulations.push_back({ triId });
	};

	// buffer geom postprocessing
	for (unsigned int i = 0; i < this->vertexIndices.size(); i += 3) {
		bufferTriangleFromIds(this->vertexIndices[i], this->vertexIndices[i + 1], this->vertexIndices[i + 2]);
		triId++;
	}

	this->vertices = std::vector<float>(geometryVertices);
}
