#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#define _USE_MATH_DEFINES

#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <set>
#include <stack>
#include <time.h>
#include "Box3.h"
#include "Matrix4.h"
#include "Vector3.h"
#include "../poly2tri/poly2tri.h"

#define uint unsigned int

namespace StructGeom {
	using Vertex = Vector3*;
	using Triangle = std::vector<Vector3*>;
	using Edge = std::vector<Vector3*>;

};

namespace BufferGeom {
	using Triangle = std::vector<uint>;
	using TriWithMarkedVertex = std::pair<Triangle, uint>;

	using Face = std::vector<Vector3*>;
	using Triangulation = std::vector<uint>;
};

using Vertex = StructGeom::Vertex;
using Tri = StructGeom::Triangle;
using Edge = StructGeom::Edge;

/*
bool operator<(const Edge& left, const Edge& right) {
	return (*left[0] < *right[0] && *left[1] < *right[1]);
}*/

enum class PrimitiveType {
	vert = 1,
	edge = 2,
	tri = 3
};

struct Primitive {
	// TODO: Test if all methods comparing size of vertices array get faster by using PrimitiveType enum
	std::vector<Vector3*> vertices = {};

	Primitive(std::vector<Vector3*> verts);
	~Primitive();

	float getMinById(uint id); // returns coord of the lowest feature vertex
	float getMaxById(uint id); // returns coord of the highest feature vertex
};

struct VertexScalarData {
	std::string name = "scalar_data";
	std::vector<float> data = {};

	VertexScalarData(std::vector<float>* data, std::string name = "scalar_data");
	~VertexScalarData();

	float& operator[](int i);
};

class Geometry
{
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
	std::vector<VertexScalarData> scalarTables = {}; // tables containing scalar data;


	// TODO: tangents and uvs

	Geometry(std::string name = "Geometry - Object");
	~Geometry();
	Geometry(const Geometry& other);

	bool hasVertices();
	bool hasUniqueVertices();
	bool hasVertexIndices();
	bool hasNormals();
	bool hasTriangulations();
	bool hasScalarData();

	Box3 getBoundingBox(Box3 bbox = Box3(), Matrix4 matrix = Matrix4());
	void computeNormals();
	void computeTriangulations();
	void fillVerticesFromUniqueVertices();

	void setScalarData(std::vector<float>* data, std::string name = "scalar_data");
	void clearScalarData();

	// getters
	std::vector<uint> getPolygonIndicesFromTriangulation(BufferGeom::Triangulation t);
	std::vector<Vector3> getProjectionsAlongNormal(BufferGeom::Face& vertices); // TODO: use Vector2
	std::vector<std::vector<uint>> getTriangulatedIndices(BufferGeom::Face& vertices);
	std::pair<std::vector<BufferGeom::Triangulation>, std::vector<size_t>> getSortedPolygonTriangulationsAndSizes();
	void getTriangles(std::vector<Tri>* trianglesBuffer);
	void getEdgesSet(std::set<Edge>* edgesSet);  // TODO: "<" operator for an edge
	std::vector<Vector3> getVertices();
	std::vector<Vector3> getUniqueVertices();
	std::vector<Primitive> getPrimitives(PrimitiveType type);

	void getVertexToTriangleMap(std::multimap<Vector3, BufferGeom::TriWithMarkedVertex>* buffer);
	void getEdgeToTriangleMap(std::multimap<Edge, BufferGeom::Triangle>* buffer);
	std::vector<Vector3> getAngleWeightedVertexPseudoNormals();
	std::vector<Vector3> getTriangleNormals();

	void getVertexToPolygonMap(std::multimap<Vector3, BufferGeom::TriWithMarkedVertex>* buffer);
	// returns vertex ring corresponding to finite volume partitions of a ring of adjacent polygons
	void getVertexFiniteVolumes(std::vector<std::vector<Vector3>>* vVolVerts, std::vector<std::vector<std::vector<uint>>>* adjacentPolyIds);
	Vector3 getNormal(BufferGeom::Face f);

	void applyMatrix(Matrix4 m);
protected:
	void clear();
private:
	void flipFaceOrientation();
	std::vector<uint> getPolygonIndicesFromTriangles(std::vector<BufferGeom::Triangle> triangles);
};

// ==== External Methods ============
// geom:
Geometry mergeGeometries(std::vector<Geometry>& geometries);
Vector3 getTriangleNormal(StructGeom::Triangle triangle, Vector3& resultNormal);

// intersections:
bool getPlaneBoxIntersection(Vector3* normal, Vector3* vert, Vector3* boxMax);
bool getTriangleBoxIntersection(Vector3** T, Vector3* boxCenter, Vector3* boxHalfSize);
bool getEdgeBoxIntersection(Edge& vertices, Vector3* boxMin, Vector3* boxMax);
bool getPrimitiveBoxIntersection(Primitive& primitive, Vector3* boxCenter, Vector3* boxMin, Vector3* boxMax, Vector3* boxHalfSize, float offset = 0.0001f);

float getRayTriangleIntersection(Vector3& rayStart, Vector3& rayDirection, Tri* tri, float minParam, float maxParam);

// distances:
float getDistanceToATriangleSq(Vector3** vertices, Vector3* point);
float getDistanceToATriangleSq2(Tri* vertices, Vector3& point);
// it's hard to say which one of the previous two is faster
float getDistanceToAnEdgeSq(Edge* vertices, Vector3& point);
float getDistanceToAPrimitiveSq(Primitive& primitive, Vector3& point);

Vector3 getClosestPtOnATriangle(Tri* vertices, Vector3& point);
Vector3 getClosestPtOnAnEdge(Edge* vertices, Vector3& point);

#endif