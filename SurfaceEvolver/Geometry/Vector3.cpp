#include "Vector3.h"
#include "Matrix4.h"
#include "Quaternion.h"

#define zeroVect 0

Vector3::Vector3()
{
	this->x = 0.;
	this->y = 0.;
	this->z = 0.;
}

Vector3::Vector3(const Vector3& other)
{
	this->x = other.x;
	this->y = other.y;
	this->z = other.z;
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

void Vector3::set(float x, float y, float z)
{
	this->x = x; this->y = y; this->z = z;
}

Vector3 Vector3::setAndReturn(float x, float y, float z)
{
	this->set(x, y, z);
	return *this;
}

void Vector3::setCoordById(float val, unsigned int id)
{
	if (id == 0) {
		x = val;
	}
	else if (id == 1) {
		y = val;
	}
	else {
		z = val;
	}
}

void Vector3::min(Vector3 other)
{
	this->x = std::fminf(this->x, other.x);
	this->y = std::fminf(this->y, other.y);
	this->z = std::fminf(this->z, other.z);
}

void Vector3::max(Vector3 other)
{
	this->x = std::fmaxf(this->x, other.x);
	this->y = std::fmaxf(this->y, other.y);
	this->z = std::fmaxf(this->z, other.z);
}

float Vector3::getCoordById(unsigned int id)
{
	if (id == 0) {
		return this->x;
	}
	else if (id == 1) {
		return this->y;
	}
	return this->z;
}

bool Vector3::equals(Vector3 other)
{
	return (
		fabs(this->x - other.x) < 2.0f * FLT_EPSILON &&
		fabs(this->y - other.y) < 2.0f * FLT_EPSILON &&
		fabs(this->z - other.z) < 2.0f * FLT_EPSILON
	);
}

bool Vector3::equalsWithEpsilon(Vector3 other, float epsilon)
{
	return (
		fabs(this->x - other.x) < epsilon &&
		fabs(this->y - other.y) < epsilon &&
		fabs(this->z - other.z) < epsilon
		);
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

void Vector3::applyQuaternion(Quaternion& q)
{
	// calculate quat * vector

	float ix = q.w * x + q.y * z - q.z * y;
	float iy = q.w * y + q.z * x - q.x * z;
	float iz = q.w * z + q.x * y - q.y * x;
	float iw = -q.x * x - q.y * y - q.z * z;

	// calculate result * inverse quat

	this->x = ix * q.w + iw * -q.x + iy * -q.z - iz * -q.y;
	this->y = iy * q.w + iw * -q.y + iz * -q.x - ix * -q.z;
	this->z = iz * q.w + iw * -q.z + ix * -q.y - iy * -q.x;
}

void Vector3::applyAxisAngle(Vector3& axis, float angle)
{
	Quaternion q = Quaternion().setFromAxisAngleAndReturn(&axis, angle);
	applyQuaternion(q);
}

void Vector3::applyMatrix4(Matrix4& m)
{
	Vector3 a = *this;
	float* e = m.elements;
	float w = 1.0f / (e[12] * a.x + e[13] * a.y + e[14] * a.z + e[15]);

	this->x = (e[0] * a.x + e[1] * a.y + e[2] * a.z + e[3]) * w;
	this->y = (e[4] * a.x + e[5] * a.y + e[6] * a.z + e[7]) * w;
	this->z = (e[8] * a.x + e[9] * a.y + e[10] * a.z + e[11]) * w;
}

void Vector3::normalize()
{
	float len = length();
	int e = zeroVect;
	try {
		if (fabs(len) < FLT_MIN) {
			throw(zeroVect);
		}
		x /= len;
		y /= len;
		z /= len;
	}
	catch (int e) {
		std::exception("Attempting to normalize a zero-length vector");
	}
}

void Vector3::lerp(Vector3 other, float param)
{
	this->x += (other.x - x) * param;
	this->y += (other.y - y) * param;
	this->z += (other.z - z) * param;
}

Vector3 Vector3::cross(Vector3 other)
{
	float ax = this->x, ay = this->y, az = this->z;
	float bx = other.x, by = other.y, bz = other.z;

	this->x = ay * bz - az * by;
	this->y = az * bx - ax * bz;
	this->z = ax * by - ay * bx;

	return *this;
}

Vector3 normalize(Vector3 target)
{
	target.normalize();
	return target;
}

float dot(Vector3 a, Vector3 b)
{
	return Vector3(a).dot(b);
}

Vector3 cross(Vector3 a, Vector3 b)
{
	Vector3 result = a;
	return a.cross(b);
}

Vector3 lerp(Vector3 v1, Vector3 v2, float param)
{
	Vector3 result = v1;
	result.lerp(v2, param);
	return result;
}

Vector3 multiply(Vector3 a, Vector3 b)
{
	a.multiply(b);
	return a;
}

bool equal(Vector3& a, Vector3& b)
{
	return a.equals(b);
}

bool notEqual(Vector3& a, Vector3& b)
{
	return !a.equals(b);
}

void Vector3::applyMatrix3(Matrix3& m)
{
	Vector3 a = *this;
	float* e = m.elements;
	this->x = e[0] * x + e[3] * y + e[6] * z;
	this->y = e[1] * x + e[4] * y + e[7] * z;
	this->z = e[2] * x + e[5] * y + e[8] * z;
}

void Vector3::addScalar(float scalar)
{
	this->x += scalar;
	this->y += scalar;
	this->z += scalar;
}

void Vector3::multiply(Vector3& other)
{
	this->x *= other.x;
	this->y *= other.y;
	this->z *= other.z;
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

Vector3& Vector3::operator+=(const Vector3& other)
{
	this->x += other.x;	this->y += other.y;	this->z += other.z;
	return *this;
}

Vector3& Vector3::operator-=(const Vector3& other)
{
	this->x -= other.x;	this->y -= other.y;	this->z -= other.z;
	return *this;
}

Vector3& Vector3::operator*=(const float& scalar)
{
	this->x *= scalar; this->y *= scalar; this->z *= scalar;
	return *this;
}

Vector3& Vector3::operator/=(const float& scalar)
{
	this->x /= scalar; this->y /= scalar; this->z /= scalar;
	return *this;
}

Vector3 operator*(Matrix3 m, Vector3 a)
{
	Vector3 result = a;
	result.applyMatrix3(m);
	return result;
}

Vector3 operator*(Matrix4 m, Vector3 a)
{
	Vector3 result = a;
	result.applyMatrix4(m);
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

bool operator<(const Vector3& left, const Vector3& right)
{
	return (
		left.x < right.x ||
		(left.x == right.x && left.y < right.y) ||
		left.x == right.x && left.y == right.y && left.z < right.z
	);
}
