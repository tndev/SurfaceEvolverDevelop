#include "Vector3.h"
#include "Matrix4.h"

#define zeroVect 0

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

float Vector3::dot(Vector3 other)
{
	return x * other.x + y * other.y + z * other.z;
}

float Vector3::lengthSq()
{
	return x * x + y * y + z * z;
}

float Vector3::length()
{
	return sqrt(lengthSq());
}

void Vector3::toArray(float* a)
{
	a[0] = x;
	a[1] = y;
	a[2] = z;
}

void Vector3::applyMatrix(Matrix4 m)
{
	float* e = m.elements;
	Vector3 a = this->clone();
	float w = 1.f / (e[3] * a.x + e[7] * a.y + e[11] * a.z + e[15]);

	x = (e[0] * a.x + e[4] * a.y + e[8] * a.z + e[12]) * w;
	y = (e[1] * a.x + e[5] * a.y + e[9] * a.z + e[13]) * w;
	z = (e[2] * a.x + e[6] * a.y + e[10] * a.z + e[14]) * w;
}

void Vector3::normalize()
{
	float len = length();
	int e = zeroVect;
	try {
		if (len < 10 * DBL_MIN && len > -10 * DBL_MIN) {
			throw(zeroVect);
		}
		x /= len;
		y /= len;
		z /= len;
	}
	catch (int e) {
		std::cout << "Attempting to normalize a zero-length vector: length(" << this << ") = " << e << std::endl;
	}
}

Vector3 normalize(Vector3 target)
{
	target.normalize();
	return target;
}

float dot(Vector3 a, Vector3 b)
{
	return a.clone().dot(b);
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

Vector3 operator*(Matrix4 m, Vector3 a)
{
	Vector3 result = a.clone();
	result.applyMatrix(m);
	return result;
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
