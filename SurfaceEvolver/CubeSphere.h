#ifndef CUBESPHERE_H_
#define CUBESPHERE_H_

#include "Geometry.h"
#include "PrimitiveBox.h"
#include "Deform.h"

class CubeSphere :
	public Geometry
{
public:
	unsigned int detail = 0;
	float radius = 1.;
	CubeSphere();
	CubeSphere(unsigned int detail, float radius);
	~CubeSphere();

	void build();
};

#endif
