#include "Icosphere.h"
#include "../Geometry/Vector3.h"
#include <vector>
#include <array>
#include <map>

IcoSphere::IcoSphere()
{
}

IcoSphere::IcoSphere(const IcoSphere& other)
{
	Geometry::Geometry(other);
	detail = other.detail;
	radius = other.radius;
}

IcoSphere::IcoSphere(unsigned int detail, float radius, std::string name)
{
	this->detail = detail; this->radius = radius;
	if (name.empty()) {
		this->name = "IcoSphere, detail: " + std::to_string(this->detail) + ", radius: " + std::to_string(this->radius);
	}
	build();
}

IcoSphere::~IcoSphere()
{
}

// ===============================================
// https://schneide.blog/2016/07/15/generating-an-icosphere-in-c/
// ======== Pre-requisites for Icosphere =========

using TriangleList = std::vector<BufferGeom::Triangle>;
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
			mid[edge] = midpointId(lookup, vertices, each[edge], each[(edge + 1) % 3]);
		}

		result.push_back({ each[0], mid[0], mid[2] });
		result.push_back({ each[1], mid[1], mid[0] });
		result.push_back({ each[2], mid[2], mid[1] });
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
			this->vertexIndices.push_back(triangles[i][j]);
		}
	}

	std::vector<float> geometryVertices = std::vector<float>(3 * this->vertexIndices.size());
	this->normals = std::vector<float>(3 * this->vertexIndices.size());


	// apply radius to duplicate vertices
	for (unsigned int i = 0; i < this->vertexIndices.size(); i++) {
		geometryVertices[(size_t)i * 3] = radius * this->vertices[(size_t)this->vertexIndices[i] * 3];
		geometryVertices[(size_t)i * 3 + 1] = radius * this->vertices[(size_t)this->vertexIndices[i] * 3 + 1];
		geometryVertices[(size_t)i * 3 + 2] = radius * this->vertices[(size_t)this->vertexIndices[i] * 3 + 2];
	}

	unsigned int triId = 0;

	auto bufferTriangleFromIds = [&](
		unsigned int i0, unsigned int i1, unsigned int i2
	) {
			this->vertexIndices[(size_t)3 * triId] = i0;
			this->vertexIndices[(size_t)3 * triId + 1] = i1;
			this->vertexIndices[(size_t)3 * triId + 2] = i2;

			if (detail > 0) {
				this->normals[(size_t)9 * triId] = this->vertices[(size_t)3 * i0];
				this->normals[(size_t)9 * triId + 1] = this->vertices[(size_t)3 * i0 + 1];
				this->normals[(size_t)9 * triId + 2] = this->vertices[(size_t)3 * i0 + 2];

				this->normals[(size_t)9 * triId + 3] = this->vertices[(size_t)3 * i1];
				this->normals[(size_t)9 * triId + 4] = this->vertices[(size_t)3 * i1 + 1];
				this->normals[(size_t)9 * triId + 5] = this->vertices[(size_t)3 * i1 + 2];

				this->normals[(size_t)9 * triId + 6] = this->vertices[(size_t)3 * i2];
				this->normals[(size_t)9 * triId + 7] = this->vertices[(size_t)3 * i2 + 1];
				this->normals[(size_t)9 * triId + 8] = this->vertices[(size_t)3 * i2 + 2];
			}
			else {
				Vector3 pA = Vector3(this->vertices[(size_t)3 * i0], this->vertices[(size_t)3 * i0 + 1], this->vertices[(size_t)3 * i0 + 2]);
				Vector3 pB = Vector3(this->vertices[(size_t)3 * i1], this->vertices[(size_t)3 * i1 + 1], this->vertices[(size_t)3 * i1 + 2]);
				Vector3 pC = Vector3(this->vertices[(size_t)3 * i2], this->vertices[(size_t)3 * i2 + 1], this->vertices[(size_t)3 * i2 + 2]);
				Vector3 centroid = Vector3();

				centroid = (pA + pB + pC) / 3.0f;

				float norm = centroid.length();

				this->normals[(size_t)9 * triId] = centroid.x / norm;
				this->normals[(size_t)9 * triId + 1] = centroid.y / norm;
				this->normals[(size_t)9 * triId + 2] = centroid.z / norm;

				this->normals[(size_t)9 * triId + 3] = centroid.x / norm;
				this->normals[(size_t)9 * triId + 4] = centroid.y / norm;
				this->normals[(size_t)9 * triId + 5] = centroid.z / norm;

				this->normals[(size_t)9 * triId + 6] = centroid.x / norm;
				this->normals[(size_t)9 * triId + 7] = centroid.y / norm;
				this->normals[(size_t)9 * triId + 8] = centroid.z / norm;
			}

			this->triangulations.push_back({ triId });
	};

	// buffer geom postprocessing
	for (unsigned int i = 0; i < this->vertexIndices.size(); i += 3) {
		bufferTriangleFromIds(this->vertexIndices[(size_t)i], this->vertexIndices[(size_t)i + 1], this->vertexIndices[(size_t)i + 2]);
		triId++;
	}

	this->vertices = std::vector<float>(geometryVertices);
}
