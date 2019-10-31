#ifndef BOX3_H_
#define BOX3_H_

#include "Vector3.h"

class Box3
{
public:
	Vector3 min = Vector3(INFINITY, INFINITY, INFINITY);
	Vector3 max = Vector3(-INFINITY, -INFINITY, -INFINITY);

	Box3();
	Box3(Vector3 min, Vector3 max);
	~Box3();
	Box3(const Box3& other);

	bool isEmpty();
	bool intersectsBox(Box3& other);

	void expandByPoint(Vector3 p);
	void expandByOffset(float offset);
	Vector3 getCenter();
	Vector3 getSize();
};

#endif
