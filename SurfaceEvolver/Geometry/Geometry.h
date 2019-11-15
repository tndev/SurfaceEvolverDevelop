#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <time.h>
#include "Box3.h"
#include "Matrix4.h"
#include "Vector3.h"
#include "../poly2tri/poly2tri.h"

#define uint unsigned int

// TODO: types in this namespace should be classes containing all necessary topological info
namespace StructGeom {
	using Triangle = std::vector<Vector3*>;
	using Edge = std::pair<Vector3*, Vector3*>;
};

namespace BufferGeom {
	using Triangle = std::vector<uint>;
	using TriWithMarkedVertex = std::pair<Triangle, uint>;
	using Face = std::vector<Vector3*>;
	using Triangulation = std::vector<uint>;
};

using Tri = StructGeom::Triangle;

class Geometry
{
private:
	std::vector<StructGeom::Triangle> triangles = {};
	std::vector<StructGeom::Edge> edges = {};
	std::multimap<Vector3, BufferGeom::TriWithMarkedVertex> vertexToTriangles;
public:
	std::string name = "Geometry - Object";

	// every time a geometry is created, it is preferrable to keep an additional array of vertices to, which
	// vertexIndices actually point
	std::vector<Vector3> uniqueVertices;
	std::vector<float> vertices; // duplicated (for each triangle)
	std::vector<float> normals; // vertex normals (for each triangle, i.e.: there's as many as there are vertex indices)
	// [0, 1, 2, 0, 2, 3, ... ] (e.g.: quads are made of 2 consecutive triplets of vert indices)
	std::vector<unsigned int> vertexIndices; // values correspond to the positions in uniqueVertices array;
	std::vector<BufferGeom::Triangulation> triangulations; // each contains ids of triangles inside a polygon

	// TODO: tangents and uvs

	Geometry(std::string name = "Geometry - Object");
	~Geometry();
	Geometry(const Geometry& other);

	bool hasVertices();
	bool hasVertexIndices();
	bool hasNormals();
	bool hasTriangulations();
	bool hasVertexToTrianglesMap();

	Box3 getBoundingBox(Box3 bbox = Box3(), Matrix4 matrix = Matrix4());
	void computeNormals();
	void computeTriangulations();

	std::vector<uint> getPolygonIndicesFromTriangulation(BufferGeom::Triangulation t);
	std::vector<Vector3> getProjectionsAlongNormal(BufferGeom::Face& vertices); // TODO: use Vector2
	std::vector<std::vector<uint>> getTriangulatedIndices(BufferGeom::Face& vertices);
	std::pair<std::vector<BufferGeom::Triangulation>, std::vector<size_t>> getSortedPolygonTriangulationsAndSizes();
	std::vector<StructGeom::Triangle> getTriangles();
	std::vector<StructGeom::Edge> getEdges();
	std::vector<Vector3> getVertices();
	std::vector<Vector3> getUniqueVertices();

	void getVertexToTriangleMap(std::multimap<Vector3, BufferGeom::TriWithMarkedVertex>* buffer);
	std::vector<Vector3> getAngleWeightedVertexPseudoNormals();
	Vector3 getAngleWeightedPseudonormalToVertex(uint vId);

	void applyMatrix(Matrix4 m);
	Vector3 getNormal(BufferGeom::Face f);

protected:
	void clear();
private:
	void flipFaceOrientation();
	std::vector<uint> getPolygonIndicesFromTriangles(std::vector<BufferGeom::Triangle> triangles);
};


Geometry mergeGeometries(std::vector<Geometry>& geometries);
Vector3 getTriangleNormal(StructGeom::Triangle triangle, Vector3 resultNormal);
bool getTriangleBoundingBoxIntersection(Tri* vertices, Vector3& bboxCenter, Vector3& bboxHalfSize, float offset = 0.0001f, Vector3* optTriNormal = nullptr);
float getDistanceToATriangleSq(Tri* vertices, Vector3& point);
float getDistanceToATriangleSq2(Tri* vertices, Vector3& point);
float getRayTriangleIntersection(Vector3& rayStart, Vector3& rayDirection, Tri* tri, float minParam, float maxParam);

#endif