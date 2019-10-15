#ifndef VECTOR3_H_
#define VECTOR3_H_
#include <iostream>

class Vector3
{
public:
	float x, y, z;

	Vector3();
	Vector3(float x, float, float z);
	~Vector3();
	void copy(Vector3 other);
	Vector3 clone();

	void negate();

	Vector3 operator+ (Vector3 other);
	Vector3 operator- (Vector3 other);
	Vector3 operator* (float scalar);
	Vector3 operator/ (float scalar);
	friend Vector3 operator*(float scalar, Vector3 a);
	friend Vector3 operator/(float scalar, Vector3 a);
	friend std::ostream& operator<< (std::ostream& out, const Vector3& v);
};

#endif

