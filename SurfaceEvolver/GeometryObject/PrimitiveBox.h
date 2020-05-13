#ifndef PRIMITIVEBOX_H_
#define PRIMITIVEBOX_H_

#include "../Geometry/Geometry.h"

class PrimitiveBox :
	public Geometry
{
public:
	double dimensions[3] = { 50.0, 50.0, 50.0 };
	unsigned int segments[3] = { 2, 2, 2 };
	bool quad = true;
	bool lastWall = true;

	PrimitiveBox();
	PrimitiveBox(const PrimitiveBox& other);
	PrimitiveBox(double x, double y, double z, unsigned int sx, unsigned int sy, unsigned int sz, bool quad = true, std::string name = "", bool lastWall = true);
	~PrimitiveBox();

	void build();
};

#endif