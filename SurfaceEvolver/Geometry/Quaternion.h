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
	float x, y, z, w;

	Quaternion();
	Quaternion(const Quaternion& other);
	Quaternion(float x, float y, float z, float w);
	~Quaternion();

	void set(float x, float y, float z, float w);

	float lengthSq();
	float length();
	float dot(Quaternion& q);
	float angleTo(Quaternion& q);

	void conjugate();
	void normalize();

	void setFromAxisAngle(Vector3* axis, float angle);
	void setFromRotationMatrix(Matrix4* m);
	Quaternion setFromAxisAngleAndReturn(Vector3* axis, float angle);
	Quaternion setFromRotationMatrixAndReturn(Matrix4* m);
	void multiplyQuaternions(Quaternion& a, Quaternion& b);
	void multiply(Quaternion& q);
	void premultiply(Quaternion& q);
};

Quaternion conjugate(Quaternion& q);
Quaternion multiply(Quaternion& a, Quaternion& b);
Quaternion normalize(Quaternion& q);
float clamp(float v, float lo, float hi);

#endif