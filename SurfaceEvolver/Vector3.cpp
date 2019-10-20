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

Vector3::Vector3(float v[3])
{
	this->x = v[0]; this->y = v[1]; this->z = v[2];
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

float Vector3::lengthSq()
{
	return x * x + y * y + z * z;
}

float Vector3::length()
{
	return sqrt(lengthSq());
}

float* Vector3::toArray()
{
	float result[3] = {x, y, z};
	return result;
}

void Vector3::normalize()
{
	float len = length();
	try {
		if (len < 10 * DBL_MIN && len > -10 * DBL_MIN) {
			throw(0);
		}
		x /= len;
		y /= len;
		z /= len;
	}
	catch (int e) {
		std::cout << "Attempting to normalize a zero-length vector: " << this << std::endl;
	}
}

Vector3 normalize(Vector3 target)
{
	target.normalize();
	return target;
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
