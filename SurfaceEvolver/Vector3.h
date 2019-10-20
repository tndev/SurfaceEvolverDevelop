#ifndef VECTOR3_H_
#define VECTOR3_H_

#include <iostream>
#include "Matrix4.h"

class Vector3
{
public:
	float x, y, z;

	Vector3();
	Vector3(float x, float, float z);
	Vector3(float v[3]);
	~Vector3();
	void copy(Vector3 other);
	Vector3 clone();

	void negate();
	float dot(Vector3 other);
	void normalize();

	float lengthSq();
	float length();

	void toArray(float* a);

	void applyMatrix(Matrix4 m);

	Vector3 operator+ (Vector3 other);
	Vector3 operator- (Vector3 other);
	Vector3 operator* (float scalar);
	Vector3 operator/ (float scalar);
	friend Vector3 operator* (Matrix4 m, Vector3 a);
	friend Vector3 operator* (float scalar, Vector3 a);
	friend Vector3 operator/ (float scalar, Vector3 a);
	friend std::ostream& operator<< (std::ostream& out, const Vector3& v);
};

Vector3 normalize(Vector3 target);
float dot(Vector3 a, Vector3 b);

#endif

