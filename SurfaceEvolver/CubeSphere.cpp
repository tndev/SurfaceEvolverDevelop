#include "CubeSphere.h"

CubeSphere::CubeSphere()
{
}

CubeSphere::CubeSphere(unsigned int detail, float radius)
{
	this->detail = detail; this->radius = radius;
	build();
}

CubeSphere::~CubeSphere()
{
}

// =======================================================
// ============ Pre-requisites for CubeSphere ============

struct Quad {
	unsigned int vertex[4];
};

using QuadList = std::vector<Quad>;
using VertexList = std::vector<Vector3>;

void CubeSphere::build()
{
	float a = 2 * radius / sqrt(3.);
	PrimitiveBox box = PrimitiveBox(a, a, a, detail, detail, detail);

	// translate to center
	Matrix4 T = Matrix4();
	T.makeTranslation(-a / 2., -a / 2., -a / 2.);
	box.applyMatrix(T);

	// spherify
	Deform def = Deform(&box);
	def.spherify(1.);

	Geometry result = def.result;
	this->uniqueVertices = result.uniqueVertices;
	this->vertices = result.vertices;
	this->normals = result.normals;
	this->vertexIndices = result.vertexIndices;
	this->triangulations = result.triangulations;
}
