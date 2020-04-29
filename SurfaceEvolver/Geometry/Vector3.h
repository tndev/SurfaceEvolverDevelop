#ifndef VECTOR3_H_
#define VECTOR3_H_

// to fix circular dependency of methods applyQuaternion and applyMatrix4
class Quaternion;
class Matrix4;

#include <iostream>
#include "Matrix3.h"

class Vector3
{
public:
	float x, y, z;

	Vector3();
	Vector3(const Vector3& other);
	Vector3(float x, float, float z);
	Vector3(float v[3]);
	~Vector3();

	void set(float x, float y, float z);
	Vector3 setAndReturn(float x, float y, float z);
	void setCoordById(float val, unsigned int id);
	void min(Vector3 other);
	void max(Vector3 other);
	float getCoordById(unsigned int id);

	bool equals(Vector3 other);
	bool equalsWithEpsilon(Vector3 other, float epsilon);

	void negate();
	float dot(Vector3 other);
	void normalize();
	void lerp(Vector3 other, float param);
	Vector3 cross(Vector3 other);

	float lengthSq();
	float length();

	void toArray(float* a);

	void applyQuaternion(Quaternion& q);
	void applyAxisAngle(Vector3& axis, float angle);
	void applyMatrix4(Matrix4& m);
	void applyMatrix3(Matrix3& m);
	void addScalar(float scalar);
	void multiply(Vector3& other);

	Vector3 operator+ (Vector3 other);
	Vector3 operator- (Vector3 other);
	Vector3 operator* (float scalar);
	Vector3 operator/ (float scalar);
	friend Vector3 operator* (Matrix3 m, Vector3 a);
	friend Vector3 operator* (Matrix4 m, Vector3 a);
	friend Vector3 operator* (float scalar, Vector3 a);
	friend Vector3 operator/ (float scalar, Vector3 a);
	Vector3& operator+=(const Vector3& other);
	Vector3& operator-=(const Vector3& other);
	Vector3& operator*=(const float& scalar);
	Vector3& operator/=(const float& scalar);
	friend std::ostream& operator<< (std::ostream& out, const Vector3& v);

	friend bool operator< (const Vector3& left, const Vector3& right);
};

Vector3 normalize(Vector3 target);
float dot(Vector3 a, Vector3 b);
Vector3 cross(Vector3 a, Vector3 b);
Vector3 lerp(Vector3 v1, Vector3 v2, float param);
Vector3 multiply(Vector3 a, Vector3 b);
bool equal(Vector3& a, Vector3& b);
bool notEqual(Vector3& a, Vector3& b);

#endif

