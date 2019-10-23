#include "Box3.h"

Box3::Box3()
{
}

Box3::~Box3()
{
}

bool Box3::isEmpty()
{
	return (
		this->min.equals(Vector3(INFINITY, INFINITY, INFINITY)) &&
		this->max.equals(Vector3(-INFINITY, -INFINITY, -INFINITY))
	);
}

void Box3::expandByPoint(Vector3 p)
{
	this->min.min(p);
	this->max.max(p);
}

void Box3::expandByOffset(float offset)
{
	this->min.addScalar(-offset);
	this->max.addScalar(offset);
}

Vector3 Box3::getCenter()
{
	if (isEmpty()) {
		return Vector3();
	}

	Vector3 center = 0.5 * (min + max);
	return center;
}

Vector3 Box3::getSize()
{
	if (isEmpty()) {
		return Vector3();
	}

	Vector3 dims = (max - min);
	return dims;
}
