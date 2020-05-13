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
		(this->min.x > FLT_MAX && this->min.y > FLT_MAX && this->min.z > FLT_MAX) && 
		(this->max.x < -FLT_MAX && this->max.y < -FLT_MAX && this->max.z < -FLT_MAX)
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

void Box3::expandByOffset(double offset)
{
	this->min.addScalar(-offset);
	this->max.addScalar(offset);
}

void Box3::expandByFactor(double factor)
{
	Vector3 scale = this->getSize();
	Vector3 scaled = factor * scale;
	Vector3 offset = 0.5 * (scaled - scale);
	this->min = this->min - offset;
	this->max = this->max + offset;
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

void Box3::setToCenter(Vector3* target)
{
	if (isEmpty()) {
		return;
	}

	target->set(
		0.5 * (this->min.x + this->max.x),
		0.5 * (this->min.y + this->max.y),
		0.5 * (this->min.z + this->max.z)
	);
}

void Box3::setToSize(Vector3* target)
{
	if (isEmpty()) {
		return;
	}

	target->set(
		(this->max.x - this->min.x),
		(this->max.y - this->min.y),
		(this->max.z - this->min.z)
	);
}

void Box3::setToHalfSize(Vector3* target)
{
	if (isEmpty()) {
		return;
	}

	target->set(
		0.5 * (this->max.x - this->min.x),
		0.5 * (this->max.y - this->min.y),
		0.5 * (this->max.z - this->min.z)
	);
}

bool Box3::equals(Box3& other)
{
	return (this->min.equals(other.min) && this->max.equals(other.max));
}

Vector3* Box3::getBoundById(unsigned int id)
{
	if (id == 0) {
		return &min;
	}
	return &max;
}

bool Box3::isInside(Vector3& pt)
{
	return (
		(pt.x >= this->min.x && pt.x <= this->max.x) &&
		(pt.y >= this->min.y && pt.y <= this->max.y) &&
		(pt.z >= this->min.z && pt.z <= this->max.z)
	);
}
