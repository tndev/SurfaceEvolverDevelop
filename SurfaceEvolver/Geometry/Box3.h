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
	void expandByFactor(float factor);
	Vector3 getCenter();
	Vector3 getSize();
	void setToCenter(Vector3* target);
	void setToSize(Vector3* target);
	void setToHalfSize(Vector3* target);
	bool equals(Box3& other);

	Vector3* getBoundById(unsigned int id);

	bool isInside(Vector3& pt);
};

#endif
