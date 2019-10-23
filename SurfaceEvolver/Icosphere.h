#ifndef ICOSPHERE_H_
#define ICOSPHERE_H_

#include "Geometry.h"

class IcoSphere :
	public Geometry
{
public:
	unsigned int detail = 0;
	float radius = 50.0f;
	IcoSphere();
	IcoSphere(const IcoSphere& other);
	IcoSphere(unsigned int detail, float radius);
	~IcoSphere();

	void build();
};

#endif