#include "Geometry.h"
#include "Quaternion.h"

#define failedTriang 0

Geometry::Geometry(std::string name)
{
	this->name = name;
}

Geometry::~Geometry()
{
	clear();
}

Geometry::Geometry(const Geometry& other)
{
	name = std::string(other.name);
	normals = std::vector<float>(other.normals);
	uniqueVertices = std::vector<Vector3>(other.uniqueVertices);
	vertices = std::vector<float>(other.vertices);
	vertexIndices = std::vector<unsigned int>(other.vertexIndices);
	triangulations = std::vector<std::vector<unsigned int>>(other.triangulations);
}

bool Geometry::hasVertices()
{
	return this->vertices.size() > 0;
}

bool Geometry::hasUniqueVertices()
{
	return this->uniqueVertices.size() > 0;
}

bool Geometry::hasVertexIndices()
{
	return this->vertexIndices.size() > 0;
}

bool Geometry::hasNormals()
{
	return this->normals.size() > 0;
}

bool Geometry::hasTriangulations()
{
	return this->triangulations.size() > 0;
}

bool Geometry::hasScalarData()
{
	return scalarTables.size();
}

Box3 Geometry::getBoundingBox(Box3 bbox, Matrix4 matrix)
{
	Vector3 helperVector = Vector3();
	for (unsigned int i = 0; i < vertices.size(); i += 3) {
		helperVector.set(vertices[i], vertices[i + 1], vertices[i + 2]);
		helperVector.applyMatrix4(matrix);
		bbox.expandByPoint(helperVector);
	}
	return bbox;
}

void Geometry::computeNormals()
{
	StructGeom::Triangle faceVerts = { &Vector3(), &Vector3(), &Vector3() };

	Vector3 normal = Vector3();
	this->normals = std::vector<float>(this->vertices.size());
	for (unsigned int i = 0; i < this->vertices.size(); i += 9) {

		for (unsigned int j = 0; j < 3; j++) {
			faceVerts[j]->x = this->vertices[i + j * 3];
			faceVerts[j]->y = this->vertices[i + j * 3 + 1];
			faceVerts[j]->z = this->vertices[i + j * 3 + 2];
		}

		normal = getTriangleNormal(faceVerts, normal);

		for (unsigned int j = 0; j < 3; j++) {
			this->normals[i + j * 3] = normal.x;
			this->normals[i + j * 3 + 1] = normal.y;
			this->normals[i + j * 3 + 2] = normal.z;
		}
	}
}

void Geometry::computeTriangulations()
{
	if (this->hasVertexIndices() || this->hasVertices()) {
		unsigned int N = (this->hasVertexIndices() && !this->hasVertices()) ? this->vertexIndices.size() : this->vertices.size() / 3;
		for (unsigned int i = 0; i < N; i += 3) {
			this->triangulations.push_back({ i / 3 });
		}
	}
	// otherwise the geometry is just a point cloud and has to be triangulated separately
}

void Geometry::fillVerticesFromUniqueVertices()
{
	if (this->hasVertexIndices() && this->hasUniqueVertices()) {
		uint N = vertexIndices.size();
		this->vertices.clear();
		this->vertices = std::vector<float>(3 * N);
		for (uint i = 0; i < N; i++) {
			this->vertices[(size_t)3 * i] = this->uniqueVertices[this->vertexIndices[i]].x;
			this->vertices[(size_t)3 * i + 1] = this->uniqueVertices[this->vertexIndices[i]].y;
			this->vertices[(size_t)3 * i + 2] = this->uniqueVertices[this->vertexIndices[i]].z;
		}
	}
}

void Geometry::setScalarData(std::vector<float>* data, std::string name)
{
	if (!this->hasVertices() || this->uniqueVertices.size() != data->size()) return;

	this->scalarTables.push_back(VertexScalarData(data, name));
}

void Geometry::clearScalarData()
{
	scalarTables.clear();
}

std::vector<unsigned int> Geometry::getPolygonIndicesFromTriangulation(BufferGeom::Triangulation t)
{
	Geometry* g = this;
	std::vector<BufferGeom::Triangle> triangles = std::vector<BufferGeom::Triangle>();
	for (unsigned int k = 0; k < t.size(); k++) {
		triangles.push_back({ g->vertexIndices[3 * t[k]], g->vertexIndices[3 * t[k] + 1], g->vertexIndices[3 * t[k] + 2] });
	}
	std::vector<unsigned int> polygonIds = g->getPolygonIndicesFromTriangles(triangles);
	return polygonIds;
}

std::vector<unsigned int> Geometry::getPolygonIndicesFromTriangles(std::vector<BufferGeom::Triangle> triangles)
{
	if (triangles.size() > 1) {
		std::vector<std::vector<unsigned int>> edges = std::vector<std::vector<unsigned int>>();

		for (auto&& t : triangles) {
			edges.push_back({ t[0], t[1] });
			edges.push_back({ t[1], t[2] });
			edges.push_back({ t[2], t[0] });
		}
		for (unsigned int i = 0; i < edges.size() - 1; i++) {
			for (unsigned int j = i + 1; j < edges.size(); j++) {
				if ((edges[i][0] == edges[j][0] && edges[i][1] == edges[j][1]) || (edges[i][0] == edges[j][1] && edges[i][1] == edges[j][0])) {
					edges.erase(edges.begin() + j);
					edges.erase(edges.begin() + i);
					i--;
					break;
				}
			}
		}

		std::vector<unsigned int> faceVertIndices = std::vector<unsigned int>();
		faceVertIndices.push_back(edges[0][1]);

		std::vector<std::vector<unsigned int>> edgesSet = std::vector<std::vector<unsigned int>>(edges);
		edgesSet.erase(edgesSet.begin());

		for (unsigned int i = 1; i < edges.size(); i++) {
			for (unsigned int j = 0; j < edgesSet.size(); j++) {
				if (edgesSet[j][0] == faceVertIndices[faceVertIndices.size() - 1]) {
					faceVertIndices.push_back(edgesSet[j][1]);
					edgesSet.erase(edgesSet.begin() + j);
					break;
				}
			}
		}

		return faceVertIndices;
	}
	else if (triangles.size() == 1) {
		return triangles[0];
	}

	return std::vector<unsigned int>();
}

std::vector<Vector3> Geometry::getVertices()
{
	std::vector<Vector3> result = std::vector<Vector3>();

	for (unsigned int i = 0; i < this->vertices.size(); i += 3) {
		result.push_back(Vector3(this->vertices[i], this->vertices[i + 1], this->vertices[i + 2]));
	}
	return result;
}

std::vector<Vector3> Geometry::getUniqueVertices()
{
	return this->uniqueVertices;
}

std::vector<Primitive> Geometry::getPrimitives(PrimitiveType type)
{
	std::vector<Primitive> result = {};
	if (type == PrimitiveType::vert) {
		for (uint i = 0; i < this->uniqueVertices.size(); i++) {
			std::vector<Vector3*> v = { &this->uniqueVertices[i] };
			result.push_back(Primitive(v));
		}
	}
	else if (type == PrimitiveType::edge) {
		std::set<Edge> edgesSet = {};
		this->getEdgesSet(&edgesSet);
		std::set<Edge>::iterator it;
		for (it = edgesSet.begin(); it != edgesSet.end(); ++it) {
			Edge e = *it;
			result.push_back(Primitive({ e[0], e[1] }));
		}
	}
	else {	
		for (uint i = 0; i < this->vertexIndices.size(); i += 3) {
			Tri T = { 
				&this->uniqueVertices[this->vertexIndices[i]],
				&this->uniqueVertices[this->vertexIndices[i + 1]],
				&this->uniqueVertices[this->vertexIndices[i + 2]]
			};

			result.push_back(Primitive(T));
		}
	}
	return result;
}

void Geometry::getVertexToTriangleMap(std::multimap<Vector3, BufferGeom::TriWithMarkedVertex>* buffer)
{
	std::multimap<Vector3, BufferGeom::TriWithMarkedVertex>::iterator it;
	for (uint i = 0; i < this->vertexIndices.size(); i += 3) {
		BufferGeom::Triangle face = { this->vertexIndices[i] , this->vertexIndices[i + 1], this->vertexIndices[i + 2] };
		for (uint j = 0; j < 3; j++) {
			it = buffer->insert(
				std::pair<Vector3, BufferGeom::TriWithMarkedVertex>(
					this->uniqueVertices[this->vertexIndices[i + j]],
					BufferGeom::TriWithMarkedVertex(face, this->vertexIndices[i + j])
				)
			);
		}
	}
}

void Geometry::getEdgeToTriangleMap(std::multimap<Edge, BufferGeom::Triangle>* buffer)
{
	std::multimap<Edge, BufferGeom::Triangle>::iterator it;
	for (uint i = 0; i < this->vertexIndices.size(); i += 3) {
		BufferGeom::Triangle T = { this->vertexIndices[i],	this->vertexIndices[i + 1],	this->vertexIndices[i + 2] };
		for (uint j = 0; j < 3; j++) {
			Edge e = {
				&this->uniqueVertices[this->vertexIndices[i + j]],
				&this->uniqueVertices[this->vertexIndices[i + (j + 1) % 3]]
			};
			it = buffer->insert(std::pair<Edge, BufferGeom::Triangle>(e, T));
		}
	}
}

std::vector<Vector3> Geometry::getAngleWeightedVertexPseudoNormals()
{
	std::vector<Vector3> result = {};

	std::multimap<Vector3, BufferGeom::TriWithMarkedVertex> vertexToTriangles = {};
	this->getVertexToTriangleMap(&vertexToTriangles);
	std::multimap<Vector3, BufferGeom::TriWithMarkedVertex>::iterator it;
	Vector3* v; std::vector<Vector3*> verts = {};
	float alpha; Vector3 triNormal = Vector3();

	for (uint i = 0; i < this->uniqueVertices.size(); i++) {
		v = &this->uniqueVertices[i];
		it = vertexToTriangles.find(*v);
		Vector3 pseudoNormal = Vector3();
		while (it != vertexToTriangles.end()) {
			BufferGeom::TriWithMarkedVertex ti = it->second;
			Tri T = { &uniqueVertices[ti.first[0]], &uniqueVertices[ti.first[1]], &uniqueVertices[ti.first[2]] };
			// compute normal
			getTriangleNormal(T, triNormal);
			triNormal.normalize();
			for (auto&& vi : ti.first) {
				if (vi != ti.second) { // ignore the marked vertex because it's v
					verts.push_back(&this->uniqueVertices[vi]);
				}
			}
			// compute angle
			alpha = acos(dot(normalize(*verts[0] - *v), normalize(*verts[1] - *v)));

			pseudoNormal = pseudoNormal + alpha * triNormal;

			verts.clear();
			vertexToTriangles.erase(it);
			it = vertexToTriangles.find(*v);
		}
		result.push_back(normalize(pseudoNormal / (2 * M_PI)));
	}

	return result;
}

using Polygon = BufferGeom::Triangle;
using PolyWithMarkedVertex = BufferGeom::TriWithMarkedVertex;
void Geometry::getVertexToPolygonMap(std::multimap<Vector3, BufferGeom::TriWithMarkedVertex>* buffer)
{
	std::multimap<Vector3, BufferGeom::TriWithMarkedVertex>::iterator it;
	for (auto&& t : triangulations) {
		Polygon face = this->getPolygonIndicesFromTriangulation(t);
		for (uint j = 0; j < face.size(); j++) {
			it = buffer->insert(
				std::pair<Vector3, BufferGeom::TriWithMarkedVertex>(
					this->uniqueVertices[face[j]],
					BufferGeom::TriWithMarkedVertex(face, face[j])
				)
			);
		}
	}
}

bool sortPolysWithMarkedVertexByAdjacency(std::vector<PolyWithMarkedVertex>* polys) {
	std::list<PolyWithMarkedVertex> result = {};
	
	PolyWithMarkedVertex currPoly, leftPoly, rightPoly;
	// take first
	std::vector<PolyWithMarkedVertex>::iterator currIt = polys->begin();
	bool CCW = true; // Counter-Clock-Wise direction
	bool closedCycle = false; // closed cycle flag
	currPoly = *currIt;
	// find iterator to marked vertex of first poly
	auto firstMarkedIt = std::find(currPoly.first.begin(), currPoly.first.end(), currPoly.second);
	// select next vertex in the first poly as the admissable edge vertex
	auto firstEdgeIt = (std::next(firstMarkedIt) == currPoly.first.end() ? currPoly.first.begin() : std::next(firstMarkedIt));
	uint firstEdgeId = *firstEdgeIt;

	while (polys->size()) {
		// copy current poly to the result list
		currPoly = *currIt;
		if (CCW) result.push_back(currPoly);
		else result.push_front(currPoly);

		// find right edge ids of the selected polygon
		// iterator to marked vertex of current poly
		auto markedIt = std::find(currPoly.first.begin(), currPoly.first.end(), currPoly.second);
		// iterator to previous vertex of current poly's marked vertex
		auto prevIt = (markedIt == currPoly.first.begin() ? currPoly.first.end() - 1 : std::prev(markedIt));
		// iterator to next vertex of current poly's marked vertex (only used for CW rotation)
		auto nextIt = (std::next(markedIt) == currPoly.first.end() ? currPoly.first.begin() : std::next(markedIt));

		// erase current poly from the list
		polys->erase(currIt);

		if (CCW) {
			// search through all remaining polys for left neighbor
			for (auto it = polys->begin(); it != polys->end(); it++) {
				leftPoly = *it;
				// find left edge ids of the left poly candidate
				auto markedLeftIt = std::find(leftPoly.first.begin(), leftPoly.first.end(), leftPoly.second);
				auto nextLeftIt = (std::next(markedLeftIt) == leftPoly.first.end() ? leftPoly.first.begin() : std::next(markedLeftIt));
			
				if (*prevIt == *nextLeftIt) {
					// edges match
					currIt = it;
					closedCycle = (firstEdgeId == *nextLeftIt);
					break;
				}

				if (std::next(it) == polys->end()) {
					// start from beginning circling clockwise if no shared edges found
					currIt = polys->begin();
					CCW = false; // start adding nodes to the front
					break;
				}
			}
		} else if (polys->size()) {
			// search through all remaining polys for right neighbor
			for (auto it = polys->end() - 1; it != polys->end(); it--) {
				rightPoly = *it;
				// find right edge ids of the right poly candidate
				auto markedRightIt = std::find(rightPoly.first.begin(), rightPoly.first.end(), rightPoly.second);
				auto prevRightIt = (markedRightIt == rightPoly.first.begin() ? std::prev(rightPoly.first.end()) : std::prev(markedRightIt));

				if (*nextIt == *prevRightIt) {
					// edges match
					currIt = it;
					break;
				}
			}
		}
	}
	*polys = std::vector<BufferGeom::TriWithMarkedVertex>(result.begin(), result.end());
	return closedCycle; // closed cycle flag
}

void Geometry::getVertexFiniteVolumes(std::vector<std::vector<Vector3>>* vVolVerts, std::vector<std::vector<Polygon>>* adjacentPolyIds)
{
	// multimap for polygons adjacent to a vertex
	std::multimap<Vector3, PolyWithMarkedVertex> vertexToPolygons = {};
	this->getVertexToPolygonMap(&vertexToPolygons);
	std::multimap<Vector3, PolyWithMarkedVertex>::iterator it;
	Vector3* v;
	std::vector<Vector3> volRingVerts = {};

	for (uint i = 0; i < this->uniqueVertices.size(); i++) {
		v = &this->uniqueVertices[i];
		it = vertexToPolygons.find(*v);

		// extract adjacent polygons from multimap
		std::vector<PolyWithMarkedVertex> adjacentPolys = {};
		while (it != vertexToPolygons.end()) {
			adjacentPolys.push_back(it->second);
			vertexToPolygons.erase(it);
			it = vertexToPolygons.find(*v);
		}
		// sort by adjacency within poly ring
		bool closedCycle = sortPolysWithMarkedVertexByAdjacency(&adjacentPolys);
		adjacentPolyIds->push_back({ }); 

		Polygon P; Polygon::iterator markedIt;
		for (uint p = 0; p < adjacentPolys.size(); p++) {
			P = adjacentPolys[p].first;
			adjacentPolyIds->at(i).push_back({ i }); // all polys will start with the central vertex F_i

			markedIt = std::find(P.begin(), P.end(), adjacentPolys[p].second);
			Polygon::iterator nextIt = (std::next(markedIt) == P.end() ? P.begin() : std::next(markedIt));
			uint leftId = *nextIt;

			Vector3 centroid = uniqueVertices[*markedIt];
			while (nextIt != markedIt) {
				centroid = centroid + uniqueVertices[*nextIt];
				adjacentPolyIds->at(i).at(p).push_back(*nextIt);
				nextIt = (std::next(nextIt) == P.end() ? P.begin() : std::next(nextIt));
			}
			centroid = 1.0f / ((float)P.size()) * centroid;

			volRingVerts.push_back(0.5f * (*v + uniqueVertices[leftId]));
			volRingVerts.push_back(centroid);
		}
		if (!closedCycle) {
			// close ring if not closed
			P = adjacentPolys[adjacentPolys.size() - 1].first;
			markedIt = std::find(P.begin(), P.end(), adjacentPolys[adjacentPolys.size() - 1].second);
			Polygon::iterator prevIt = (markedIt == P.begin() ? P.end() - 1 : std::prev(markedIt));
			volRingVerts.push_back(0.5f * (*v + uniqueVertices[*prevIt]));
		}

		// the result should be a contour of the finite volume around v
		vVolVerts->push_back(volRingVerts);
		volRingVerts.clear();
	}
}

std::vector<Vector3> Geometry::getTriangleNormals()
{
	std::vector<Vector3> result = {};
	Vector3 triNorm = Vector3();

	for (uint i = 0; i < vertexIndices.size(); i += 3) {
		Tri T = { 
			&uniqueVertices[vertexIndices[i]],
			&uniqueVertices[vertexIndices[i + 1]],
			&uniqueVertices[vertexIndices[i + 2]]
		};
		getTriangleNormal(T, triNorm);
		Vector3 n = normalize(triNorm);
		result.push_back(n);
	}

	return result;
}


std::vector<Vector3> Geometry::getProjectionsAlongNormal(BufferGeom::Face& vertices)
{
	Vector3 normal = getNormal(vertices); // directions
	Vector3 referencePoint = *vertices[0];
	Vector3 up = Vector3(0, 0, 1);

	// ------ -----------
	// TODO: create Quaternion

	if (fabs(dot(normal, up)) >= 0.999) {
		if (fabs(normal.x) < fabs(normal.y) && fabs(normal.x) < fabs(normal.z)) {
			up = normalize(Vector3(1, 0, 0).cross(normal));
		}
		else if (fabs(normal.y) < fabs(normal.z)) {
			up = normalize(Vector3(0, 1, 0).cross(normal));
		}
		else {
			up = normalize(Vector3(0, 0, 1).cross(normal));
		}
	}

	Vector3 crossProduct = normalize(cross(up, normal));
	Vector3 pseudoUp = normalize(cross(normal, crossProduct));

	Matrix4 transform = Matrix4(
		normal.x, normal.y, normal.z, 0,
		crossProduct.x, crossProduct.y, crossProduct.z, 0,
		pseudoUp.x, pseudoUp.y, pseudoUp.z, 0,
		0, 0, 0, 1
	);
	transform.transpose();

	Vector3 projX = transform * Vector3(0, 1, 0);
	Vector3 projY = transform * Vector3(0, 0, 1);

	std::vector<Vector3> projections = std::vector<Vector3>();
	for (unsigned int k = 0; k < vertices.size(); k++) {
		Vector3 vec = *vertices[k] - referencePoint;
		float px = dot(vec, projX);
		float py = dot(vec, projY);
		projections.push_back(Vector3(px, py, 0.0f));
	}

	return projections;
}

std::vector<std::vector<unsigned int>> Geometry::getTriangulatedIndices(BufferGeom::Face& vertices)
{
	std::vector<std::vector<unsigned int>> faces = std::vector<std::vector<unsigned int>>();
	if (vertices.size() < 3) {
	}
	else if (vertices.size() == 3) {
		faces = { {0, 1, 2} };
	}
	else if (vertices.size() == 4) {
		faces = { {0, 1, 2}, {0, 2, 3} };
		Vector3 e2 = *vertices[2] - *vertices[0];
		Vector3 e1 = *vertices[1] - *vertices[0];
		Vector3 e3 = *vertices[3] - *vertices[0];
		Vector3 c21 = cross(e2, e1);
		Vector3 c23 = cross(e2, e3);
		if (dot(c21, c23) > 0.0) {
			faces = { {0, 1, 3}, {1, 2, 3} };
		}
	}

	else {

		std::vector<Vector3> projections = getProjectionsAlongNormal(vertices);

		p2t::CDT* cdt = NULL;

		// poly2tri doesn't work well with duplicate or collinear points. it is possible that it will fail on certain inputs
		// if it fails, we retry with increasingly higher jitter. it usually works on first retry
		// if it doesn't work with 42 retries, there is no hope
		std::vector<p2t::Point> points = {};

		unsigned int MAX_TRIES = 42;
		unsigned int tries_utilized = 0;
		srand((unsigned)time(NULL));
		for (tries_utilized = 0; tries_utilized < MAX_TRIES; tries_utilized++) {
			//// poly2tri triangulation, seems to work well with holes
			points.clear();
			for (unsigned int i = 0; i < projections.size(); i++) {
				p2t::Point pt = p2t::Point(
					projections[i].x + (tries_utilized * (rand() / RAND_MAX * 0.001 - 0.0005)),
					projections[i].y + (tries_utilized * (rand() / RAND_MAX * 0.001 - 0.0005))
				);
				pt._vIdx = i;
				points.push_back(pt);
			}
			
			try {
				std::vector<p2t::Point*> contour = {};
				for (unsigned int i = 0; i < points.size(); i++) {
					contour.push_back(&points[i]);
				}
				cdt = new p2t::CDT(contour);
				cdt->Triangulate();
				if (cdt == NULL) {
					throw(failedTriang);
				}
			}
			catch (int e) {
				continue;
			}
			break;// assuming it all goes well
		}
		if (cdt == NULL) {
			std::cout << "ERROR! Triangulation was not able to finish with" << MAX_TRIES << "tries" << std::endl;
			return {};
		}

		std::vector<p2t::Triangle*> triangles = cdt->GetTriangles();

		for (unsigned int i = 0; i < triangles.size(); i++) {
			faces.push_back({
				triangles[i]->GetPoint(0)->_vIdx,
				triangles[i]->GetPoint(1)->_vIdx,
				triangles[i]->GetPoint(2)->_vIdx
			});
		}
	}

	return faces;
}

template<class T>
void flipArray(std::vector<T> arr, unsigned int elementSize) {
	if (arr.size() && elementSize) {
		unsigned int iStep = 3 * elementSize;
		float tmp[3] = { 0.0f, 0.0f, 0.0f };
		for (unsigned int i = 0; i < arr.size(); i += iStep) {
			for (unsigned int j = 0; j < elementSize; j++) {
				const unsigned int i0 = i + j;
				const unsigned int i1 = i + j + 2 * elementSize;
				tmp[j] = arr[i0];
				arr[i0] = arr[i1];
				arr[i1] = tmp[j];
			}
		}
	}	
}

void Geometry::flipFaceOrientation()
{
	flipArray(vertices, 3);
	if (hasNormals()) {
		flipArray(normals, 3);
	}
	flipArray(vertexIndices, 1);
}

Matrix3 getSubMatrix3(Matrix4 m) {
	float* e = m.elements;
	float resultElems[9] = {
		e[0],	e[4],	e[8],
		e[1],	e[5],	e[9],
		e[2],	e[6],	e[10]
	};
	return Matrix3(resultElems);
}

void Geometry::applyMatrix(Matrix4 m)
{
	if (m.isIdentity()) {
		return;
	}

	const bool hasNormals = this->hasNormals();
	const bool changeOrientation = (m.determinant() < 0);
	Matrix4 mInverse = Matrix4().getInverse(m);
	Matrix4 tmp = transpose(mInverse);
	Matrix3 matrix3 = m.getSubMatrix3();
	Matrix3 matrix3InverseTransposed = getSubMatrix3(tmp);

	Vector3 helperVector = Vector3();

	for (int i = 0; i < vertices.size(); i += 3) {
		const unsigned int i0 = i, i1 = i + 1, i2 = i + 2;

		helperVector.set(vertices[i0], vertices[i1], vertices[i2]);
		helperVector.applyMatrix4(m);

		vertices[i0] = helperVector.x;
		vertices[i1] = helperVector.y;
		vertices[i2] = helperVector.z;

		this->uniqueVertices[this->vertexIndices[i / 3]] = helperVector;

		if (hasNormals) {
			helperVector.set(normals[i0], normals[i1], normals[i2]);
			helperVector.applyMatrix3(matrix3InverseTransposed);
			const float length = helperVector.length();

			normals[i0] = helperVector.x / length;
			normals[i1] = helperVector.y / length;
			normals[i2] = helperVector.z / length;
		}
	}

	if (changeOrientation) {
		flipFaceOrientation();
	}
}

Vector3 Geometry::getNormal(BufferGeom::Face f)
{
	Vector3 normal = Vector3();
	for (unsigned int i = 0; i < f.size(); i++) {
		Vector3 fromNext = *f[i] - *f[(i + 1) % f.size()];
		Vector3 plusNext = *f[i] + *f[(i + 1) % f.size()];
		normal.x = normal.x + fromNext.y * plusNext.z;
		normal.y = normal.y + fromNext.z * plusNext.x;
		normal.z = normal.z + fromNext.x * plusNext.y;
	}
	return normalize(normal);
}

void Geometry::clear()
{
	uniqueVertices.clear();
	vertices.clear();
	normals.clear();
	vertexIndices.clear();
	triangulations.clear();
	clearScalarData();
}

std::pair<std::vector<BufferGeom::Triangulation>, std::vector<size_t>> Geometry::getSortedPolygonTriangulationsAndSizes()
{
	std::vector<BufferGeom::Triangulation> T = triangulations;
	std::vector<size_t> sizes = std::vector<size_t>();
	std::sort(T.begin(), T.end(), [](const BufferGeom::Triangulation& a, const BufferGeom::Triangulation& b) { return a.size() < b.size(); });
	unsigned int tCount = 0; // triangulation count
	size_t currentSize = T.begin()->size();

	for (std::vector<BufferGeom::Triangulation>::iterator it = T.begin(); it < T.end(); ++it) {
		if (it->size() > currentSize) {
			currentSize = it->size();
			sizes.push_back(tCount);
			tCount = 0;
			continue;
		}
		tCount++;
	}
	sizes.push_back(tCount);

	std::pair<std::vector<BufferGeom::Triangulation>, std::vector<size_t>> result = { T, sizes };
	return result;
}

void Geometry::getTriangles(std::vector<Tri>* trianglesBuffer)
{
	for (unsigned int i = 0; i < this->vertexIndices.size(); i += 3) {
		Vector3* v0 = &this->uniqueVertices[this->vertexIndices[i]];
		Vector3* v1 = &this->uniqueVertices[this->vertexIndices[i + 1]];
		Vector3* v2 = &this->uniqueVertices[this->vertexIndices[i + 2]];

		Tri T = { v0, v1, v2 };

		trianglesBuffer->push_back(T);
	}
}

// TODO: "<" operator for an edge
void Geometry::getEdgesSet(std::set<Edge>* edgesSet)
{
	std::pair<std::set<Edge>::iterator, bool> it;
	for (uint i = 0; i < this->vertexIndices.size(); i += 3) {
		BufferGeom::Triangle T = { this->vertexIndices[i],	this->vertexIndices[i + 1],	this->vertexIndices[i + 2] };
		for (uint j = 0; j < 3; j++) {
			Edge e = {
				&this->uniqueVertices[this->vertexIndices[i + j]],
				&this->uniqueVertices[this->vertexIndices[i + (j + 1) % 3]]
			};
			it = edgesSet->insert(e);
		}
	}
}

Geometry mergeGeometries(std::vector<Geometry>& geometries)
{
	Geometry result = Geometry();
	unsigned int vertexCount = 0;
	unsigned int uniqueVertCount = 0;

	for (unsigned int i = 0; i < geometries.size(); i++) {
		if (!geometries[i].hasNormals()) {
			geometries[i].computeNormals();
		}
		if (!geometries[i].hasTriangulations()) {
			geometries[i].computeTriangulations();
		}
		vertexCount += geometries[i].vertices.size() / 3;
		uniqueVertCount += geometries[i].uniqueVertices.size();
	};

	result.uniqueVertices = {};
	result.vertices = {};
	result.normals = {};
	result.vertexIndices = std::vector<unsigned int>(vertexCount);
	result.triangulations = {};

	unsigned int idx = 0;
	unsigned int uIdx = 0;
	for (unsigned int i = 0; i < geometries.size(); i++) {
		const unsigned int currVertexCount = geometries[i].vertices.size() / 3;
		const unsigned int currUniqueCount = geometries[i].uniqueVertices.size();

		result.uniqueVertices.insert(result.uniqueVertices.end(), geometries[i].uniqueVertices.begin(), geometries[i].uniqueVertices.end());
		result.vertices.insert(result.vertices.end(), geometries[i].vertices.begin(), geometries[i].vertices.end());
		result.normals.insert(result.normals.end(), geometries[i].normals.begin(), geometries[i].normals.end());

		for (unsigned int j = 0; j < currVertexCount; j++) {
			const unsigned int vertexIdx = geometries[i].hasVertexIndices() ? geometries[i].vertexIndices[j] : j;
			result.vertexIndices[j + idx] = vertexIdx + uIdx;
		}

		if (geometries[i].hasTriangulations()) {
			for (unsigned int j = 0; j < geometries[i].triangulations.size(); j++) {
				BufferGeom::Triangulation triangulationCopy = geometries[i].triangulations[j];
				for (unsigned int k = 0; k < triangulationCopy.size(); k++) {
					triangulationCopy[k] += idx / 3;
				}
				result.triangulations.push_back(triangulationCopy);
			}

		}

		idx += currVertexCount;
		uIdx += currUniqueCount;
	}
	return result;
}

Vector3 getTriangleNormal(StructGeom::Triangle triangle, Vector3& resultNormal)
{
	Vector3 u = *triangle[1] - *triangle[0];
	Vector3 v = *triangle[2] - *triangle[0];

	resultNormal.set(
		(u.y * v.z) - (u.z * v.y),
		(u.z * v.x) - (u.x * v.z),
		(u.x * v.y) - (u.y * v.x)
	);
	return resultNormal;
}

// ====== Helper macros for vectors (to increase speed) ============

#define CROSS(dest, v1, v2)						\
          dest.x = v1.y * v2.z - v1.z * v2.y;   \
          dest.y = v1.z * v2.x - v1.x * v2.z;	\
          dest.z = v1.x * v2.y - v1.y * v2.x;

#define DOT(v1, v2) (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z)

#define SUB(dest, v1, v2)						\
          dest.x = v1.x - v2.x;					\
          dest.y = v1.y - v2.y;					\
          dest.z = v1.z - v2.z;

#define FINDMINMAX(x0, x1, x2, min, max)		\
		  min = max = x0;						\
		  if (x1 < min) min = x1;				\
		  if (x1 > max) max = x1;				\
		  if (x2 < min) min = x2;				\
		  if (x2 > max) max = x2;

bool getPlaneBoxIntersection(Vector3* normal, Vector3* vert, Vector3* boxMax) {
	int q;
	Vector3 vmin, vmax;
	float v;

	for (q = 0; q <= 2; q++) {
		v = vert->getCoordById(q);

		if (normal->getCoordById(q) > 0.0f) {
			vmin.setCoordById(-boxMax->getCoordById(q) - v, q);
			vmax.setCoordById(boxMax->getCoordById(q) - v, q);
		} else {
			vmin.setCoordById(boxMax->getCoordById(q) - v, q);
			vmax.setCoordById(-boxMax->getCoordById(q) - v, q);
		}
	}

	if (DOT((*normal), vmin) > 0.0f) return false;
	if (DOT((*normal), vmax) >= 0.0f) return true;
	return false;
}

// =================== Helper macros for axis tests ======================
// X-tests:

#define AXISTEST_X01(a, b, fa, fb)									   \
		p0 = a * v0.y - b * v0.z;			       				       \
		p2 = a * v2.y - b * v2.z;			       					   \
        if (p0 < p2) {min = p0; max = p2;} else {min = p2; max = p0;}  \
		rad = fa * boxHalfSize->y + fb * boxHalfSize->z;		       \
		if (min > rad || max < -rad) return false;

#define AXISTEST_X2(a, b, fa, fb)									   \
		p0 = a * v0.y - b * v0.z;									   \
		p1 = a * v1.y - b * v1.z;		       						   \
        if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;}  \
		rad = fa * boxHalfSize->y + fb * boxHalfSize->z;			   \
		if (min > rad || max < -rad) return false;

// Y-tests:

#define AXISTEST_Y02(a, b, fa, fb)									   \
		p0 = -a * v0.x + b * v0.z;				   					   \
		p2 = -a * v2.x + b * v2.z;				       				   \
		if (p0 < p2) {min = p0; max = p2;} else {min = p2; max = p0;}  \
		rad = fa * boxHalfSize->x + fb * boxHalfSize->z;		       \
		if (min > rad || max <- rad) return false;

#define AXISTEST_Y1(a, b, fa, fb)							    	   \
		p0 = -a * v0.x + b * v0.z;		      					       \
		p1 = -a * v1.x + b * v1.z;					                   \
		if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;}  \
		rad = fa * boxHalfSize->x + fb * boxHalfSize->z;			   \
		if (min > rad || max < -rad) return false;

// Z-tests:

#define AXISTEST_Z12(a, b, fa, fb)									   \
		p1 = a * v1.x - b * v1.y;			 			               \
		p2 = a * v2.x - b * v2.y;			       	                   \
		if (p2 < p1) {min = p2; max = p1;} else {min = p1; max = p2;}  \
		rad = fa * boxHalfSize->x + fb * boxHalfSize->y;			   \
		if (min > rad || max < -rad) return false;

#define AXISTEST_Z0(a, b, fa, fb)									   \
		p0 = a * v0.x - b * v0.y;									   \
		p1 = a * v1.x - b * v1.y;									   \
		if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;}  \
		rad = fa * boxHalfSize->x + fb * boxHalfSize->y;			   \
		if (min > rad || max < -rad) return false;

// =====================================================

bool getTriangleBoxIntersection(Vector3** T, Vector3* boxCenter, Vector3* boxHalfSize) {
	Vector3 v0, v1, v2;
	float min, max, p0, p1, p2, rad, fex, fey, fez;
	Vector3 normal, e0, e1, e2;

	SUB(v0, (*T[0]), (*boxCenter));
	SUB(v1, (*T[1]), (*boxCenter));
	SUB(v2, (*T[2]), (*boxCenter));

	// tri edges:
	SUB(e0, v1, v0);
	SUB(e1, v2, v1);
	SUB(e2, v0, v2);

	// 9 axis tests:
	fex = fabsf(e0.x);
	fey = fabsf(e0.y);
	fez = fabsf(e0.z);

	AXISTEST_X01(e0.z, e0.y, fez, fey);
	AXISTEST_Y02(e0.z, e0.x, fez, fex);
	AXISTEST_Z12(e0.y, e0.x, fey, fex);

	fex = fabsf(e1.x);
	fey = fabsf(e1.y);
	fez = fabsf(e1.z);

	AXISTEST_X01(e1.z, e1.y, fez, fey);
	AXISTEST_Y02(e1.z, e1.x, fez, fex);
	AXISTEST_Z0(e1.y, e1.x, fey, fex);

	fex = fabsf(e2.x);
	fey = fabsf(e2.y);
	fez = fabsf(e2.z);

	AXISTEST_X2(e2.z, e2.y, fez, fey);
	AXISTEST_Y1(e2.z, e2.x, fez, fex);
	AXISTEST_Z12(e2.y, e2.x, fey, fex);

	// test for AABB overlap in x, y, and z:
	// test in x:
	FINDMINMAX(v0.x, v1.x, v2.x, min, max);
	if (min > boxHalfSize->x || max < -boxHalfSize->x) return false;

	// test in y:
	FINDMINMAX(v0.y, v1.y, v2.y, min, max);
	if (min > boxHalfSize->y || max < -boxHalfSize->y) return false;

	// test in z:
	FINDMINMAX(v0.z, v1.z, v2.z, min, max);
	if (min > boxHalfSize->z || max < -boxHalfSize->z) return false;

	// test if the box intersects the triangle plane  dot(normal, x) + d = 0

	CROSS(normal, e0, e1);

	if (!getPlaneBoxIntersection(&normal, &v0, boxHalfSize)) return false;

	return true;
}

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}


bool getEdgeBoxIntersection(Edge& vertices, Vector3* boxMin, Vector3* boxMax)
{
	Vector3 edge_direction = *vertices[1] - *vertices[0];
	float edge_param = edge_direction.length();
	edge_direction = edge_direction / edge_param;

	float t_min = -FLT_MAX;
	float t_max = FLT_MAX;

	uint i; float t1, t2;

	for (i = 0; i < 3; i++) {
		if (edge_direction.getCoordById(i) > FLT_EPSILON) {
			t1 = (boxMin->getCoordById(i) - vertices[0]->getCoordById(i)) * (1.0f / edge_direction.getCoordById(i));
			t2 = (boxMax->getCoordById(i) - vertices[0]->getCoordById(i)) * (1.0f / edge_direction.getCoordById(i));

			t_min = std::max(t_min, std::min(t1, t2));
			t_max = std::min(t_max, std::max(t1, t2));
		} else if (
			vertices[0]->getCoordById(i) < boxMin->getCoordById(i) ||
			vertices[0]->getCoordById(i) > boxMax->getCoordById(i)
		) {
			return false;
		}
	}

	return (t_max >= t_min && t_max >= 0.0f && t_min <= edge_param);
}

bool getPrimitiveBoxIntersection(Primitive& primitive, Vector3* boxCenter, Vector3* boxMin, Vector3* boxMax, Vector3* boxHalfSize, float offset)
{
	if (primitive.vertices.size() == 1) {
		return (
			(primitive.vertices[0]->x >= boxMin->x && primitive.vertices[0]->x <= boxMax->x) &&
			(primitive.vertices[0]->y >= boxMin->y && primitive.vertices[0]->y <= boxMax->y) &&
			(primitive.vertices[0]->z >= boxMin->z && primitive.vertices[0]->z <= boxMax->z)
		);
	}
	else if (primitive.vertices.size() == 2) {
		return getEdgeBoxIntersection(primitive.vertices, boxMin, boxMax);
	}
	else if (primitive.vertices.size() == 3) {
		Vector3** t = new Vector3*[3];
		t[0] = primitive.vertices[0];
		t[1] = primitive.vertices[1];
		t[2] = primitive.vertices[2];
		return getTriangleBoxIntersection(t, boxCenter, boxHalfSize);
		delete[] t;
	}

	return false;
}

float getDistanceToATriangleSq(Vector3** vertices, Vector3* point)
{
	Vector3 diff = *point - *vertices[0];
	Vector3 edge0 = *vertices[1] - *vertices[0];
	Vector3 edge1 = *vertices[2] - *vertices[0];
	float a00 = DOT(edge0, edge0);
	float a01 = DOT(edge0, edge1);
	float a11 = DOT(edge1, edge1);
	float b0 = -DOT(diff, edge0);
	float b1 = -DOT(diff, edge1);
	float const zero = (float)0;
	float const one = (float)1;
	float det = a00 * a11 - a01 * a01;
	float t0 = a01 * b1 - a11 * b0;
	float t1 = a01 * b0 - a00 * b1;

	if (t0 + t1 <= det)	{
		if (t0 < zero) {
			if (t1 < zero) { // region 4			
				if (b0 < zero) {
					t1 = zero;
					if (-b0 >= a00) { // V1					
						t0 = one;
					} else { // E01					
						t0 = -b0 / a00;
					}
				} else {
					t0 = zero;
					if (b1 >= zero) { // V0					
						t1 = zero;
					} else if (-b1 >= a11) { // V2					
						t1 = one;
					} else { // E20					
						t1 = -b1 / a11;
					}
				}
			} else { // region 3			
				t0 = zero;
				if (b1 >= zero) { // V0				
					t1 = zero;
				} else if (-b1 >= a11) { // V2				
					t1 = one;
				} else { // E20				
					t1 = -b1 / a11;
				}
			}
		} else if (t1 < zero) { // region 5		
			t1 = zero;
			if (b0 >= zero) { // V0			
				t0 = zero;
			} else if (-b0 >= a00) { // V1			
				t0 = one;
			} else { // E01			
				t0 = -b0 / a00;
			}
		} else { // region 0, interior		
			float invDet = one / det;
			t0 *= invDet;
			t1 *= invDet;
		}
	} else {	
		float tmp0, tmp1, numer, denom;

		if (t0 < zero) { // region 2		
			tmp0 = a01 + b0;
			tmp1 = a11 + b1;
			if (tmp1 > tmp0) {			
				numer = tmp1 - tmp0;
				denom = a00 - ((float)2) * a01 + a11;
				if (numer >= denom) { // V1				
					t0 = one;
					t1 = zero;
				} else { // E12				
					t0 = numer / denom;
					t1 = one - t0;
				}
			} else {
				t0 = zero;
				if (tmp1 <= zero) { // V2				
					t1 = one;
				} else if (b1 >= zero) { // V0				
					t1 = zero;
				} else { // E20				
					t1 = -b1 / a11;
				}
			}
		} else if (t1 < zero) { // region 6		
			tmp0 = a01 + b1;
			tmp1 = a00 + b0;
			if (tmp1 > tmp0) {			
				numer = tmp1 - tmp0;
				denom = a00 - ((float)2) * a01 + a11;
				if (numer >= denom) { // V2				
					t1 = one;
					t0 = zero;
				} else { // E12				
					t1 = numer / denom;
					t0 = one - t1;
				}
			} else {
				t1 = zero;
				if (tmp1 <= zero) { // V1				
					t0 = one;
				} else if (b0 >= zero) { // V0				
					t0 = zero;
				} else { // E01				
					t0 = -b0 / a00;
				}
			}
		} else { // region 1		
			numer = a11 + b1 - a01 - b0;
			if (numer <= zero) { // V2			
				t0 = zero;
				t1 = one;
			} else {			
				denom = a00 - ((float)2) * a01 + a11;
				if (numer >= denom) { // V1				
					t0 = one;
					t1 = zero;
				} else { // 12				
					t0 = numer / denom;
					t1 = one - t0;
				}
			}
		}
	}

	Vector3 closest = *vertices[0] + t0 * edge0 + t1 * edge1;
	SUB(diff, (*point), closest);
	return DOT(diff, diff);
}

int const sign(float x)
{
	return (x > 0 ? 1 : -1);
}

float dot2(Vector3 v) { return dot(v, v); }

float getDistanceToATriangleSq2(Tri* vertices, Vector3& point)
{
	Vector3 ba = *vertices->at(1) - *vertices->at(0); Vector3 pa = point - *vertices->at(0);
	Vector3 cb = *vertices->at(2) - *vertices->at(1); Vector3 pb = point - *vertices->at(1);
	Vector3 ac = *vertices->at(0) - *vertices->at(2); Vector3 pc = point - *vertices->at(2);
	Vector3 nor = cross(ba, ac);

	return ((sign(dot(cross(ba, nor), pa)) +
			 sign(dot(cross(cb, nor), pb)) +
			 sign(dot(cross(ac, nor), pc)) < 2.0) ?		
		std::fminf(
			std::fminf(
				dot2(ba * clamp( dot(ba, pa) / dot2(ba), 0.0f, 1.0f) - pa),
				dot2(cb * clamp( dot(cb, pb) / dot2(cb), 0.0f, 1.0f) - pb)
			),
		dot2(ac * clamp( dot(ac, pc) / dot2(ac), 0.0f, 1.0f) - pc)
		) : dot(nor, pa) * dot(nor, pa) / dot2(nor));
}

float getDistanceToAnEdgeSq(Edge* vertices, Vector3& point)
{
	float h = (*vertices->at(1) - *vertices->at(0)).length();
	Vector3 p = point;
	p.y -= clamp(p.y, 0.0f, h);
	return p.length();
}

float getDistanceToAPrimitiveSq(Primitive& primitive, Vector3& point)
{
	if (primitive.vertices.size() == 1) {
		return (point - *primitive.vertices[0]).lengthSq();
	}
	else if (primitive.vertices.size() == 2) {
		return getDistanceToAnEdgeSq(&primitive.vertices, point);
	}
	else if (primitive.vertices.size() == 3) {
		Vector3** t = new Vector3 * [3];
		t[0] = primitive.vertices[0];
		t[1] = primitive.vertices[1];
		t[2] = primitive.vertices[2];
		return getDistanceToATriangleSq(t, &point);
		delete[] t;
	}

	std::cout << "invalid primitive vert count!" << std::endl;
	return -1.0f;
}

Vector3 getClosestPtOnATriangle(Tri* vertices, Vector3& point)
{
	Vector3 diff = point - *vertices->at(0);
	Vector3 edge0 = *vertices->at(1) - *vertices->at(0);
	Vector3 edge1 = *vertices->at(2) - *vertices->at(0);
	float a00 = dot(edge0, edge0);
	float a01 = dot(edge0, edge1);
	float a11 = dot(edge1, edge1);
	float b0 = -dot(diff, edge0);
	float b1 = -dot(diff, edge1);
	float const zero = (float)0;
	float const one = (float)1;
	float det = a00 * a11 - a01 * a01;
	float t0 = a01 * b1 - a11 * b0;
	float t1 = a01 * b0 - a00 * b1;

	if (t0 + t1 <= det)
	{
		if (t0 < zero)
		{
			if (t1 < zero)  // region 4
			{
				if (b0 < zero)
				{
					t1 = zero;
					if (-b0 >= a00)  // V1
					{
						t0 = one;
					}
					else  // E01
					{
						t0 = -b0 / a00;
					}
				}
				else
				{
					t0 = zero;
					if (b1 >= zero)  // V0
					{
						t1 = zero;
					}
					else if (-b1 >= a11)  // V2
					{
						t1 = one;
					}
					else  // E20
					{
						t1 = -b1 / a11;
					}
				}
			}
			else  // region 3
			{
				t0 = zero;
				if (b1 >= zero)  // V0
				{
					t1 = zero;
				}
				else if (-b1 >= a11)  // V2
				{
					t1 = one;
				}
				else  // E20
				{
					t1 = -b1 / a11;
				}
			}
		}
		else if (t1 < zero)  // region 5
		{
			t1 = zero;
			if (b0 >= zero)  // V0
			{
				t0 = zero;
			}
			else if (-b0 >= a00)  // V1
			{
				t0 = one;
			}
			else  // E01
			{
				t0 = -b0 / a00;
			}
		}
		else  // region 0, interior
		{
			float invDet = one / det;
			t0 *= invDet;
			t1 *= invDet;
		}
	}
	else
	{
		float tmp0, tmp1, numer, denom;

		if (t0 < zero)  // region 2
		{
			tmp0 = a01 + b0;
			tmp1 = a11 + b1;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a00 - ((float)2) * a01 + a11;
				if (numer >= denom)  // V1
				{
					t0 = one;
					t1 = zero;
				}
				else  // E12
				{
					t0 = numer / denom;
					t1 = one - t0;
				}
			}
			else
			{
				t0 = zero;
				if (tmp1 <= zero)  // V2
				{
					t1 = one;
				}
				else if (b1 >= zero)  // V0
				{
					t1 = zero;
				}
				else  // E20
				{
					t1 = -b1 / a11;
				}
			}
		}
		else if (t1 < zero)  // region 6
		{
			tmp0 = a01 + b1;
			tmp1 = a00 + b0;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a00 - ((float)2) * a01 + a11;
				if (numer >= denom)  // V2
				{
					t1 = one;
					t0 = zero;
				}
				else  // E12
				{
					t1 = numer / denom;
					t0 = one - t1;
				}
			}
			else
			{
				t1 = zero;
				if (tmp1 <= zero)  // V1
				{
					t0 = one;
				}
				else if (b0 >= zero)  // V0
				{
					t0 = zero;
				}
				else  // E01
				{
					t0 = -b0 / a00;
				}
			}
		}
		else  // region 1
		{
			numer = a11 + b1 - a01 - b0;
			if (numer <= zero)  // V2
			{
				t0 = zero;
				t1 = one;
			}
			else
			{
				denom = a00 - ((float)2) * a01 + a11;
				if (numer >= denom)  // V1
				{
					t0 = one;
					t1 = zero;
				}
				else  // 12
				{
					t0 = numer / denom;
					t1 = one - t0;
				}
			}
		}
	}

	return *vertices->at(0) + t0 * edge0 + t1 * edge1;
}

Vector3 getClosestPtOnAnEdge(Edge* vertices, Vector3& point)
{
	float h = (*vertices->at(1) - *vertices->at(0)).length();
	Vector3 p = point;
	p.y -= clamp(p.y, 0.0f, h);

	return p;
}


// https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
float getRayTriangleIntersection(Vector3& rayStart, Vector3& rayDirection, Tri* tri, float minParam, float maxParam)
{
	const float EPSILON = 0.0000001;
	Vector3 edge1, edge2, h, s, q;
	float a, f, u, v;
	edge1 = *tri->at(1) - *tri->at(0);
	edge2 = *tri->at(2) - *tri->at(0);
	h = cross(rayDirection, edge2);
	a = dot(edge1, h);
	if (a > -EPSILON && a < EPSILON) {
		return -1.0f;    // This ray is parallel to this triangle.
	}
		
	f = 1.0 / a;
	s = rayStart - *tri->at(0);
	u = f * dot(s, h);
	if (u < 0.0 || u > 1.0) {
		return -1.0f;
	}
	q = cross(s, edge1);
	v = f * dot(rayDirection, q);
	if (v < 0.0 || u + v > 1.0) {
		return -1.0f;
	}
		
	// At this stage we can compute t to find out where the intersection point is on the line.
	float t = f * dot(edge2, q);
	bool inOrNoRange = t < minParam;
	bool outOrNoRange = t > maxParam;
	if (t > EPSILON && t < 1 / EPSILON && !inOrNoRange && !outOrNoRange) {
		return t;
	}
	else {	// This means that there is a line intersection but not a ray intersection.
		return -1.0f;
	}
}

Primitive::Primitive(std::vector<Vector3*> verts)
{
	this->vertices = verts;
}

Primitive::~Primitive()
{
}

float Primitive::getMinById(uint id)
{
	if (vertices.size() == 1) {
		return vertices[0]->getCoordById(id);
	}
	else if (vertices.size() == 2) {
		return std::fminf(vertices[0]->getCoordById(id), vertices[1]->getCoordById(id));
	}
	else if (vertices.size() == 3) {
		return std::fminf(
			vertices[0]->getCoordById(id),
			std::fminf(
				vertices[1]->getCoordById(id),
				vertices[2]->getCoordById(id)
			)
		);
	}

	std::cout << "invalid primitive vertex count!" << std::endl;
	return 0.0f;
}

float Primitive::getMaxById(uint id)
{
	if (vertices.size() == 1) {
		return vertices[0]->getCoordById(id);
	}
	else if (vertices.size() == 2) {
		return std::fmaxf(vertices[0]->getCoordById(id), vertices[1]->getCoordById(id));
	}
	else if (vertices.size() == 3) {
		return std::fmaxf(
			vertices[0]->getCoordById(id),
			std::fmaxf(
				vertices[1]->getCoordById(id),
				vertices[2]->getCoordById(id)
			)
		);
	}

	std::cout << "invalid primitive vertex count!" << std::endl;
	return 0.0f;
}

VertexScalarData::VertexScalarData(std::vector<float>* data, std::string name)
{
	this->data = std::vector<float>(*data);
	this->name = name;
}

VertexScalarData::~VertexScalarData()
{
}

float& VertexScalarData::operator[](int i)
{
	return data[i];
}
