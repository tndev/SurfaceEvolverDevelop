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

void CubeSphere::copy(CubeSphere other)
{
	Geometry::copy(other);
	detail = other.detail;
	radius = other.radius;
}

CubeSphere CubeSphere::clone()
{
	CubeSphere result = CubeSphere();
	result.copy(*this);
	return result;
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

}
