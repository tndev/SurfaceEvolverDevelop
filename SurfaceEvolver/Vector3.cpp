#include "Vector3.h"

Vector3::Vector3()
{
	this->x = 0.;
	this->y = 0.;
	this->z = 0.;
}

Vector3::Vector3(float x, float y, float z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

Vector3::~Vector3()
{
}

void Vector3::copy(Vector3 other)
{
	this->x = other.x;
	this->y = other.y;
	this->z = other.z;
}

Vector3 Vector3::clone()
{
	Vector3* result = new Vector3();
	result->copy(*this);
	return *result;
}

void Vector3::negate()
{
	this->x = -x;
	this->y = -y;
	this->z = -z;
}

Vector3 Vector3::operator+(Vector3 other)
{
	return Vector3(x + other.x, y + other.y, z + other.z);
}

Vector3 Vector3::operator-(Vector3 other)
{
	return Vector3(x - other.x, y - other.y, z - other.z);
}

Vector3 Vector3::operator*(float scalar)
{
	return Vector3(scalar * x, scalar * y, scalar * z);
}

Vector3 Vector3::operator/(float scalar)
{
	return Vector3(x / scalar, y / scalar, z / scalar);
}

Vector3 operator*(float scalar, Vector3 a)
{
	return Vector3(scalar * a.x, scalar * a.y, scalar * a.z);
}

Vector3 operator/(float scalar, Vector3 a)
{
	return Vector3(a.x / scalar, a.y / scalar, a.z / scalar);
}

std::ostream& operator<<(std::ostream& out, const Vector3& v)
{
	out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
	return out;
}
