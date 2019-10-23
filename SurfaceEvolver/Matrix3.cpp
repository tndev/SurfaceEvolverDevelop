#include "Matrix3.h"

#define id(i, j) 3 * i + j
#define singularity 0

Matrix3::Matrix3()
{
}

Matrix3::Matrix3(const Matrix3& other)
{
	const float* e = other.elements;
	this->set(
		e[0], e[1], e[2],
		e[3], e[4], e[5],
		e[6], e[7], e[8]
	);
}

Matrix3::Matrix3(float* elems)
{
	float* e = elements;
	e[0] = elems[0];	e[1] = elems[1];	e[2] = elems[2];
	e[3] = elems[3];	e[4] = elems[4];	e[5] = elems[5];
	e[6] = elems[6];	e[7] = elems[7];	e[8] = elems[8];
}

Matrix3::~Matrix3()
{
}

bool Matrix3::isIdentity()
{
	float* e = elements;
	return (
		e[0] == 1.0f && e[1] == 0.0f && e[2] == 0.0f &&
		e[3] == 0.0f &&	e[4] == 1.0f && e[5] == 0.0f && 
		e[6] == 0.0f && e[7] == 0.0f &&	e[8] == 1.0f
	);
}

void Matrix3::setToIdentity()
{
	this->set(
		1.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 1.0f
	);
}

void Matrix3::set(float m00, float m01, float m02, float m10, float m11, float m12, float m20, float m21, float m22)
{
	float* e = elements;
	e[0] = m00;		e[1] = m01;		e[2] = m02;
	e[3] = m10;		e[4] = m11;		e[5] = m12;
	e[6] = m20;		e[7] = m21;		e[8] = m22;
}

void Matrix3::transpose()
{

	float* e = this->elements;
	this->set(
		e[0], e[6], e[9],
		e[1], e[4], e[7],
		e[2], e[5], e[8]
	);
}

void Matrix3::multiplyScalar(float scalar)
{
	float* e = this->elements;
	for (int i = 0; i < 9; i++) {
		e[i] *= scalar;
	}
}

float Matrix3::determinant()
{
	float* m = elements;
	float value = 
		m[id(0, 0)] * m[id(1, 1)] * m[id(2, 2)] - m[id(0, 2)] * m[id(1, 1)] * m[id(2, 0)] +
		m[id(0, 1)] * m[id(1, 2)] * m[id(2, 0)] - m[id(0, 1)] * m[id(1, 0)] * m[id(2, 2)] +
		m[id(0, 2)] * m[id(1, 0)] * m[id(2, 1)] - m[id(0, 0)] * m[id(1, 2)] * m[id(2, 1)];
	return value;
}

Matrix3 Matrix3::getInverse(Matrix3& from)
{
	float* m = from.elements;
	float* e = elements;
	float det = from.determinant();
	Matrix3 result = Matrix3();
	int ex = singularity;
	try {
		if (fabs(det) < FLT_MIN) {
			throw singularity;
		}
		else {
			float invDet = 1.0f / from.determinant();

			e[id(0, 0)] = m[id(1, 1)] * m[id(2, 2)] - m[id(2, 1)] * m[id(1, 2)];
			e[id(0, 1)] = m[id(0, 2)] * m[id(2, 1)] - m[id(0, 1)] * m[id(2, 2)];
			e[id(0, 2)] = m[id(0, 1)] * m[id(1, 2)] - m[id(0, 2)] * m[id(1, 1)];
			e[id(1, 0)] = m[id(1, 2)] * m[id(2, 0)] - m[id(1, 0)] * m[id(2, 2)];
			e[id(1, 1)] = m[id(0, 0)] * m[id(2, 2)] - m[id(0, 2)] * m[id(2, 0)];
			e[id(1, 2)] = m[id(1, 0)] * m[id(0, 2)] - m[id(0, 0)] * m[id(1, 2)];
			e[id(2, 0)] = m[id(1, 0)] * m[id(2, 1)] - m[id(2, 0)] * m[id(1, 1)];
			e[id(2, 1)] = m[id(2, 0)] * m[id(0, 1)] - m[id(0, 0)] * m[id(2, 1)];
			e[id(2, 2)] = m[id(0, 0)] * m[id(1, 1)] - m[id(1, 0)] * m[id(0, 1)];
			result = *this;
			return invDet * result;
		}
	}
	catch (int ex) {
		std::cout << "det(A) = 0.0! Unable to invert a singular matrix!" << std::endl;
	}

	return result;
}

Matrix3 Matrix3::operator+(Matrix3 other)
{
	Matrix3 result = Matrix3();
	for (int i = 0; i < 9; i++) {
		result.elements[i] = elements[i] + other.elements[i];
	}
	return result;
}

Matrix3 Matrix3::operator-(Matrix3 other)
{
	Matrix3 result = Matrix3();
	for (int i = 0; i < 9; i++) {
		result.elements[i] = elements[i] - other.elements[i];
	}
	return result;
}

Matrix3 Matrix3::operator*(Matrix3 other)
{
	Matrix3 result = Matrix3();
	float* re = result.elements;
	float* ae = elements;
	float* be = other.elements;

	float a11 = ae[0], a12 = ae[3], a13 = ae[6];
	float a21 = ae[1], a22 = ae[4], a23 = ae[7];
	float a31 = ae[2], a32 = ae[5], a33 = ae[8];

	float b11 = be[0], b12 = be[3], b13 = be[6];
	float b21 = be[1], b22 = be[4], b23 = be[7];
	float b31 = be[2], b32 = be[5], b33 = be[8];

	re[0] = a11 * b11 + a12 * b21 + a13 * b31;
	re[3] = a11 * b12 + a12 * b22 + a13 * b32;
	re[6] = a11 * b13 + a12 * b23 + a13 * b33;

	re[1] = a21 * b11 + a22 * b21 + a23 * b31;
	re[4] = a21 * b12 + a22 * b22 + a23 * b32;
	re[7] = a21 * b13 + a22 * b23 + a23 * b33;

	re[2] = a31 * b11 + a32 * b21 + a33 * b31;
	re[5] = a31 * b12 + a32 * b22 + a33 * b32;
	re[8] = a31 * b13 + a32 * b23 + a33 * b33;

	return result;
}

Matrix3 operator*(float scalar, Matrix3 m)
{
	Matrix3 result = m;
	result.multiplyScalar(scalar);
	return result;
}

Matrix3 inverse(Matrix3& m)
{
	Matrix3 result = Matrix3().getInverse(m);
	return result;
}

Matrix3 transpose(Matrix3& m)
{
	Matrix3 result = m;
	result.transpose();
	return result;
}
