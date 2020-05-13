#ifndef QUATERNION_H_
#define QUATERNION_H_

// to fix circular dependency of methods setFromAxisAngle and setFromRotationMatrix
class Matrix4;
class Vector3;

#include <cassert>
#include <cmath>

class Quaternion
{
public:
	double x, y, z, w;

	Quaternion();
	Quaternion(const Quaternion& other);
	Quaternion(double x, double y, double z, double w);
	~Quaternion();

	void set(double x, double y, double z, double w);

	double lengthSq();
	double length();
	double dot(Quaternion& q);
	double angleTo(Quaternion& q);

	void conjugate();
	void normalize();

	void setFromAxisAngle(Vector3* axis, double angle);
	void setFromRotationMatrix(Matrix4* m);
	Quaternion setFromAxisAngleAndReturn(Vector3* axis, double angle);
	Quaternion setFromRotationMatrixAndReturn(Matrix4* m);
	void multiplyQuaternions(Quaternion& a, Quaternion& b);
	void multiply(Quaternion& q);
	void premultiply(Quaternion& q);
};

Quaternion conjugate(Quaternion& q);
Quaternion multiply(Quaternion& a, Quaternion& b);
Quaternion normalize(Quaternion& q);
double clamp(double v, double lo, double hi);

#endif