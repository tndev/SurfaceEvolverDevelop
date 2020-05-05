#ifndef PRIMITIVEBOX_H_
#define PRIMITIVEBOX_H_

#include "../Geometry/Geometry.h"

class PrimitiveBox :
	public Geometry
{
public:
	float dimensions[3] = { 50.0f, 50.0f, 50.0f };
	unsigned int segments[3] = { 2, 2, 2 };
	bool quad = true;
	bool lastWall = true;

	PrimitiveBox();
	PrimitiveBox(const PrimitiveBox& other);
	PrimitiveBox(float x, float y, float z, unsigned int sx, unsigned int sy, unsigned int sz, bool quad = true, std::string name = "", bool lastWall = true);
	~PrimitiveBox();

	void build();
};

#endif