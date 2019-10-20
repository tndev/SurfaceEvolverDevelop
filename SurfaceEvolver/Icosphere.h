#ifndef ICOSPHERE_H_
#define ICOSPHERE_H_

#include "Geometry.h"

class Icosphere :
	public Geometry
{
public:
	unsigned int detail = 0;
	float radius = 1.;
	Icosphere();
	Icosphere(unsigned int detail, float radius);
	~Icosphere();

	void copy(Icosphere other);
	Icosphere clone();

	void build();
};

#endif