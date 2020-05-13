#include "Quaternion.h"
#include "Vector3.h"
#include "Matrix4.h"

Quaternion::Quaternion()
{
	this->x = 0.0; this->y = 0.0; this->z = 0.0; this->w = 1.0;
}

Quaternion::Quaternion(const Quaternion& other)
{
	set(other.x, other.y, other.z, other.w);
}

Quaternion::Quaternion(double x, double y, double z, double w)
{
	set(x, y, z, w);
}

Quaternion::~Quaternion()
{
}

void Quaternion::set(double x, double y, double z, double w)
{
	this->x = x; this->y = y; this->z = z; this->w = w;
}

double Quaternion::lengthSq()
{
	return this->x * this->x + this->y * this->y + this->z * this->z + this->w * this->w;
}

double Quaternion::length()
{
	return sqrt(this->x * this->x + this->y * this->y + this->z * this->z + this->w * this->w);
}

double Quaternion::dot(Quaternion& q)
{
	return this->x * q.x + this->y * q.y + this->z * q.z + this->w * q.w;
}

double clamp(double v, double lo, double hi)
{
	assert(!(hi < lo));
	return (v < lo) ? lo : (hi < v) ? hi : v;
}

double Quaternion::angleTo(Quaternion& q)
{
	return 2.0 * acos(fabs(clamp(this->dot(q), -1.0, 1.0)));
}

void Quaternion::conjugate()
{
	this->x *= -1.0;
	this->y *= -1.0;
	this->z *= -1.0;
}

void Quaternion::normalize()
{
	double len = length();
	if (len == 0.0f) {
		this->x = 0.0;
		this->y = 0.0;
		this->z = 0.0;
		this->w = 1.0;
	}
	else {
		len = 1.0f / len;

		this->x *= len;
		this->y *= len;
		this->z *= len;
		this->w *= len;
	}
}

void Quaternion::setFromAxisAngle(Vector3* axis, double angle)
{
	// http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm
	// assumes axis is normalized

	double halfAngle = angle / 2.0f, s = sin(halfAngle);

	this->x = axis->x * s;
	this->y = axis->y * s;
	this->z = axis->z * s;
	this->w = cos(halfAngle);
}

void Quaternion::setFromRotationMatrix(Matrix4* m)
{
	// http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
	// assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

	double* e = m->elements;

	double m11 = e[0], m12 = e[4], m13 = e[8],
		  m21 = e[1], m22 = e[5], m23 = e[9],
		  m31 = e[2], m32 = e[6], m33 = e[10];

	double trace = m11 + m22 + m33, s;

	if (trace > 0) {
		s = 0.5 / sqrt(trace + 1.0);

		this->w = 0.25 / s;
		this->x = (m32 - m23) * s;
		this->y = (m13 - m31) * s;
		this->z = (m21 - m12) * s;
	}
	else if (m11 > m22 && m11 > m33) {
		s = 2.0 * sqrt(1.0 + m11 - m22 - m33);

		this->w = (m32 - m23) / s;
		this->x = 0.25 * s;
		this->y = (m12 + m21) / s;
		this->z = (m13 + m31) / s;
	}
	else if (m22 > m33) {
		s = 2.0 * sqrt(1.0 + m22 - m11 - m33);

		this->w = (m13 - m31) / s;
		this->x = (m12 + m21) / s;
		this->y = 0.25 * s;
		this->z = (m23 + m32) / s;
	}
	else {
		s = 2.0 * sqrt(1.0 + m33 - m11 - m22);

		this->w = (m21 - m12) / s;
		this->x = (m13 + m31) / s;
		this->y = (m23 + m32) / s;
		this->z = 0.25 * s;
	}
}

Quaternion Quaternion::setFromAxisAngleAndReturn(Vector3* axis, double angle)
{
	// assumes axis is normalized

	double halfAngle = angle / 2.0f, s = sin(halfAngle);

	this->x = axis->x * s;
	this->y = axis->y * s;
	this->z = axis->z * s;
	this->w = cos(halfAngle);

	return *this;
}

Quaternion Quaternion::setFromRotationMatrixAndReturn(Matrix4* m)
{
	// assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

	double* e = m->elements;

	double m11 = e[0], m12 = e[4], m13 = e[8];
	double m21 = e[1], m22 = e[5], m23 = e[9];
	double m31 = e[2], m32 = e[6], m33 = e[10];

	double trace = m11 + m22 + m33, s;

	if (trace > 0) {
		s = 0.5 / sqrt(trace + 1.0);

		this->w = 0.25 / s;
		this->x = (m32 - m23) * s;
		this->y = (m13 - m31) * s;
		this->z = (m21 - m12) * s;
	}
	else if (m11 > m22&& m11 > m33) {
		s = 2.0 * sqrt(1.0 + m11 - m22 - m33);

		this->w = (m32 - m23) / s;
		this->x = 0.25 * s;
		this->y = (m12 + m21) / s;
		this->z = (m13 + m31) / s;
	}
	else if (m22 > m33) {
		s = 2.0 * sqrt(1.0 + m22 - m11 - m33);

		this->w = (m13 - m31) / s;
		this->x = (m12 + m21) / s;
		this->y = 0.25 * s;
		this->z = (m23 + m32) / s;
	}
	else {
		s = 2.0 * sqrt(1.0 + m33 - m11 - m22);

		this->w = (m21 - m12) / s;
		this->x = (m13 + m31) / s;
		this->y = (m23 + m32) / s;
		this->z = 0.25 * s;
	}

	return *this;
}

void Quaternion::multiplyQuaternions(Quaternion& a, Quaternion& b)
{
	// from http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/code/index.htm

	double qax = a.x, qay = a.y, qaz = a.z, qaw = a.w;
	double qbx = b.x, qby = b.y, qbz = b.z, qbw = b.w;

	this->x = qax * qbw + qaw * qbx + qay * qbz - qaz * qby;
	this->y = qay * qbw + qaw * qby + qaz * qbx - qax * qbz;
	this->z = qaz * qbw + qaw * qbz + qax * qby - qay * qbx;
	this->w = qaw * qbw - qax * qbx - qay * qby - qaz * qbz;
}

void Quaternion::multiply(Quaternion& q)
{
	return this->multiplyQuaternions(*this, q);
}

void Quaternion::premultiply(Quaternion& q)
{
	return this->multiplyQuaternions(q, *this);
}

Quaternion conjugate(Quaternion& q)
{
	return Quaternion(-q.x, -q.y, -q.z, q.w);
}

Quaternion multiply(Quaternion& a, Quaternion& b)
{
	Quaternion result;
	result.multiplyQuaternions(a, b);
	return result;
}

Quaternion normalize(Quaternion& q)
{
	q.normalize();
	return q;
}
