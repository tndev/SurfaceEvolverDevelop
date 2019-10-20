#ifndef MATRIX3_H_
#define MATRIX3_H_

#include <iostream>

class Matrix3
{
public:
	float elements[9] = {
		1.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 1.0f
	};

	Matrix3();
	Matrix3(float* elems);
	~Matrix3();

	void copy(Matrix3 other);
	Matrix3 clone();

	bool isIdentity();
	void setToIdentity();
	void transpose();
	void multiplyScalar(float scalar);
	float determinant();
	Matrix3 getInverse(Matrix3& from);

	Matrix3 operator+ (Matrix3 other);
	Matrix3 operator- (Matrix3 other);
	Matrix3 operator* (Matrix3 other);
	friend Matrix3 operator*(float scalar, Matrix3 m);
};

Matrix3 inverse(Matrix3& m);
Matrix3 transpose(Matrix3& m);

#endif