#ifndef MATRIX4_H_
#define MATRIX4_H_

// there will be circular dependencies upon classes Quaternion and Vector3 
// thanks to compose and decompose methods
#include <iostream>
#include "Vector3.h"
#include "Matrix3.h"
#include "Quaternion.h"

class Matrix4
{
public:
	double elements[16] = {
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0
	};

	Matrix4();
	Matrix4(const Matrix4& other);
	Matrix4(double* elems);
	Matrix4(
		double m00, double m01, double m02, double m03,
		double m10, double m11, double m12, double m13,
		double m20, double m21, double m22, double m23,
		double m30, double m31, double m32, double m33
	);
	~Matrix4();

	void set(
		double m00, double m01, double m02, double m03,
		double m10, double m11, double m12, double m13,
		double m20, double m21, double m22, double m23,
		double m30, double m31, double m32, double m33
	);

	Matrix3 getSubMatrix3();

	bool isIdentity();
	void setToIdentity();
	void compose(Vector3* position, Quaternion* quaternion, Vector3* scale);

	// the following methods do not return a value
	// use external methods instead
	void transpose();
	void multiplyScalar(double scalar);
	void multiplyMatrices(Matrix4& a, Matrix4& b);
	void multiply(Matrix4& m);
	Matrix4 multiply(Matrix4 m);
	void premultiply(Matrix4& m);
	Matrix4 premultiply(Matrix4 m);

	Matrix4 setToScale(double sx, double sy, double sz);
	Matrix4 makeRotationAxis(double ax, double ay, double az, double angle);
	Matrix4 makeTranslation(double tx, double ty, double tz);

	void decompose(Vector3* position, Quaternion* quaternion, Vector3* scale);

	double determinant();
	Matrix4 getInverse(Matrix4& from);

	Matrix4 operator+ (Matrix4 other);
	Matrix4 operator- (Matrix4 other);
	Matrix4 operator* (Matrix4 other);
	friend Matrix4 operator*(double scalar, Matrix4 m);
};

Matrix4 inverse(Matrix4& m);
Matrix4 transpose(Matrix4& m);

#endif