#ifndef CUBESPHERE_H_
#define CUBESPHERE_H_

#include "../Geometry/Geometry.h"
#include "PrimitiveBox.h"
#include "Deform.h"

class CubeSphere :
	public Geometry
{
public:
	unsigned int detail = 0;
	float radius = 1.;
	bool quad = true;

	CubeSphere();
	CubeSphere(unsigned int detail, float radius, bool quad = true, std::string name = "");
	~CubeSphere();

	void build();
};

#endif
