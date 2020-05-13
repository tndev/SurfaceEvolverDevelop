#ifndef MATRIX3_H_
#define MATRIX3_H_

#include <iostream>

class Matrix3
{
public:
	double elements[9] = {
		1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0
	};

	Matrix3();
	Matrix3(const Matrix3& other);
	Matrix3(double* elems);
	~Matrix3();

	bool isIdentity();
	void setToIdentity();
	void set(
		double m00, double m01, double m02,
		double m10, double m11, double m12,
		double m20, double m21, double m22
	);

	void transpose();
	void multiplyScalar(double scalar);

	double determinant();
	Matrix3 getInverse(Matrix3& from);

	Matrix3 operator+ (Matrix3 other);
	Matrix3 operator- (Matrix3 other);
	Matrix3 operator* (Matrix3 other);
	friend Matrix3 operator*(double scalar, Matrix3 m);
};

Matrix3 inverse(Matrix3& m);
Matrix3 transpose(Matrix3& m);

#endif