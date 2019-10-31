#include "Box3.h"

Box3::Box3()
{
}

Box3::Box3(Vector3 min, Vector3 max)
{
	this->min = min;
	this->max = max;
}

Box3::~Box3()
{
}

Box3::Box3(const Box3& other)
{
	this->min = other.min;
	this->max = other.max;
}

bool Box3::isEmpty()
{
	return (
		this->min.equals(Vector3(INFINITY, INFINITY, INFINITY)) &&
		this->max.equals(Vector3(-INFINITY, -INFINITY, -INFINITY))
	);
}

bool Box3::intersectsBox(Box3& other)
{
	return !(
		other.max.x < this->min.x || other.min.x > this->max.x ||
		other.max.y < this->min.y || other.min.y > this->max.y ||
		other.max.z < this->min.z || other.min.z > this->max.z
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
