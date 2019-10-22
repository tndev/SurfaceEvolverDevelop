#ifndef DEFORM_H_
#define DEFORM_H_

#include "Geometry.h"
#include "Box3.h"

class Deform
{
public:
	Geometry* geom = nullptr;
	Geometry result = Geometry();

	Deform();
	Deform(Geometry* geom);
	~Deform();

	void spherify(float param);
};

#endif