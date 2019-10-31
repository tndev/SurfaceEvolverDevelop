#ifndef MATRIX4_H_
#define MATRIX4_H_

#include <iostream>
#include "Matrix3.h"

class Matrix4
{
public:
	float elements[16] = {
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};

	Matrix4();
	Matrix4(const Matrix4& other);
	Matrix4(float* elems);
	Matrix4(
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33
	);
	~Matrix4();

	void set(
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33
	);

	Matrix3 getSubMatrix3();

	bool isIdentity();
	void setToIdentity();

	// the following methods do not return a value
	// use external methods instead
	void transpose();
	void multiplyScalar(float scalar);

	Matrix4 setToScale(float sx, float sy, float sz);
	Matrix4 makeRotationAxis(float ax, float ay, float az, float angle);
	Matrix4 makeTranslation(float tx, float ty, float tz);

	float determinant();
	Matrix4 getInverse(Matrix4& from);

	Matrix4 operator+ (Matrix4 other);
	Matrix4 operator- (Matrix4 other);
	Matrix4 operator* (Matrix4 other);
	friend Matrix4 operator*(float scalar, Matrix4 m);
};

Matrix4 inverse(Matrix4& m);
Matrix4 transpose(Matrix4& m);

#endif