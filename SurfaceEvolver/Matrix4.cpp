#include "Matrix4.h"

#define id(i, j) 4 * i + j

float Matrix4::determinant()
{
	float* m = elements;
	float value =
		m[id(0, 3)] * m[id(1, 2)] * m[id(2, 1)] * m[id(3, 0)] - m[id(0, 2)] * m[id(1, 3)] * m[id(2, 1)] * m[id(3, 0)] - m[id(0, 3)] * m[id(1, 1)] * m[id(2, 2)] * m[id(3, 0)] + m[id(0, 1)] * m[id(1, 3)] * m[id(2, 2)] * m[id(3, 0)] +
		m[id(0, 2)] * m[id(1, 1)] * m[id(2, 3)] * m[id(3, 0)] - m[id(0, 1)] * m[id(1, 2)] * m[id(2, 3)] * m[id(3, 0)] - m[id(0, 3)] * m[id(1, 2)] * m[id(2, 0)] * m[id(3, 1)] + m[id(0, 2)] * m[id(1, 3)] * m[id(2, 0)] * m[id(3, 1)] +
		m[id(0, 3)] * m[id(1, 0)] * m[id(2, 2)] * m[id(3, 1)] - m[id(0, 0)] * m[id(1, 3)] * m[id(2, 2)] * m[id(3, 1)] - m[id(0, 2)] * m[id(1, 0)] * m[id(2, 3)] * m[id(3, 1)] + m[id(0, 0)] * m[id(1, 2)] * m[id(2, 3)] * m[id(3, 1)] +
		m[id(0, 3)] * m[id(1, 1)] * m[id(2, 0)] * m[id(3, 2)] - m[id(0, 1)] * m[id(1, 3)] * m[id(2, 0)] * m[id(3, 2)] - m[id(0, 3)] * m[id(1, 0)] * m[id(2, 1)] * m[id(3, 2)] + m[id(0, 0)] * m[id(1, 3)] * m[id(2, 1)] * m[id(3, 2)] +
		m[id(0, 1)] * m[id(1, 0)] * m[id(2, 3)] * m[id(3, 2)] - m[id(0, 0)] * m[id(1, 1)] * m[id(2, 3)] * m[id(3, 2)] - m[id(0, 2)] * m[id(1, 1)] * m[id(2, 0)] * m[id(3, 3)] + m[id(0, 1)] * m[id(1, 2)] * m[id(2, 0)] * m[id(3, 3)] +
		m[id(0, 2)] * m[id(1, 0)] * m[id(2, 1)] * m[id(3, 3)] - m[id(0, 0)] * m[id(1, 2)] * m[id(2, 1)] * m[id(3, 3)] - m[id(0, 1)] * m[id(1, 0)] * m[id(2, 2)] * m[id(3, 3)] + m[id(0, 0)] * m[id(1, 1)] * m[id(2, 2)] * m[id(3, 3)];
	return 0.0f;
}

Matrix4 Matrix4::getInverse(Matrix4& from)
{
	float* m = from.elements;
	float* e = elements;
	float invDet = 1.0f / from.determinant();

	e[id(0, 0)] = m[id(1, 2)] * m[id(2, 3)] * m[id(3, 1)] - m[id(1, 3)] * m[id(2, 2)] * m[id(3, 1)] + m[id(1, 3)] * m[id(2, 1)] * m[id(3, 2)] - m[id(1, 1)] * m[id(2, 3)] * m[id(3, 2)] - m[id(1, 2)] * m[id(2, 1)] * m[id(3, 3)] + m[id(1, 1)] * m[id(2, 2)] * m[id(3, 3)];
	e[id(0, 1)] = m[id(0, 3)] * m[id(2, 2)] * m[id(3, 1)] - m[id(0, 2)] * m[id(2, 3)] * m[id(3, 1)] - m[id(0, 3)] * m[id(2, 1)] * m[id(3, 2)] + m[id(0, 1)] * m[id(2, 3)] * m[id(3, 2)] + m[id(0, 2)] * m[id(2, 1)] * m[id(3, 3)] - m[id(0, 1)] * m[id(2, 2)] * m[id(3, 3)];
	e[id(0, 2)] = m[id(0, 2)] * m[id(1, 3)] * m[id(3, 1)] - m[id(0, 3)] * m[id(1, 2)] * m[id(3, 1)] + m[id(0, 3)] * m[id(1, 1)] * m[id(3, 2)] - m[id(0, 1)] * m[id(1, 3)] * m[id(3, 2)] - m[id(0, 2)] * m[id(1, 1)] * m[id(3, 3)] + m[id(0, 1)] * m[id(1, 2)] * m[id(3, 3)];
	e[id(0, 3)] = m[id(0, 3)] * m[id(1, 2)] * m[id(2, 1)] - m[id(0, 2)] * m[id(1, 3)] * m[id(2, 1)] - m[id(0, 3)] * m[id(1, 1)] * m[id(2, 2)] + m[id(0, 1)] * m[id(1, 3)] * m[id(2, 2)] + m[id(0, 2)] * m[id(1, 1)] * m[id(2, 3)] - m[id(0, 1)] * m[id(1, 2)] * m[id(2, 3)];
	e[id(1, 0)] = m[id(1, 3)] * m[id(2, 2)] * m[id(3, 0)] - m[id(1, 2)] * m[id(2, 3)] * m[id(3, 0)] - m[id(1, 3)] * m[id(2, 0)] * m[id(3, 2)] + m[id(1, 0)] * m[id(2, 3)] * m[id(3, 2)] + m[id(1, 2)] * m[id(2, 0)] * m[id(3, 3)] - m[id(1, 0)] * m[id(2, 2)] * m[id(3, 3)];
	e[id(1, 1)] = m[id(0, 2)] * m[id(2, 3)] * m[id(3, 0)] - m[id(0, 3)] * m[id(2, 2)] * m[id(3, 0)] + m[id(0, 3)] * m[id(2, 0)] * m[id(3, 2)] - m[id(0, 0)] * m[id(2, 3)] * m[id(3, 2)] - m[id(0, 2)] * m[id(2, 0)] * m[id(3, 3)] + m[id(0, 0)] * m[id(2, 2)] * m[id(3, 3)];
	e[id(1, 2)] = m[id(0, 3)] * m[id(1, 2)] * m[id(3, 0)] - m[id(0, 2)] * m[id(1, 3)] * m[id(3, 0)] - m[id(0, 3)] * m[id(1, 0)] * m[id(3, 2)] + m[id(0, 0)] * m[id(1, 3)] * m[id(3, 2)] + m[id(0, 2)] * m[id(1, 0)] * m[id(3, 3)] - m[id(0, 0)] * m[id(1, 2)] * m[id(3, 3)];
	e[id(1, 3)] = m[id(0, 2)] * m[id(1, 3)] * m[id(2, 0)] - m[id(0, 3)] * m[id(1, 2)] * m[id(2, 0)] + m[id(0, 3)] * m[id(1, 0)] * m[id(2, 2)] - m[id(0, 0)] * m[id(1, 3)] * m[id(2, 2)] - m[id(0, 2)] * m[id(1, 0)] * m[id(2, 3)] + m[id(0, 0)] * m[id(1, 2)] * m[id(2, 3)];
	e[id(2, 0)] = m[id(1, 1)] * m[id(2, 3)] * m[id(3, 0)] - m[id(1, 3)] * m[id(2, 1)] * m[id(3, 0)] + m[id(1, 3)] * m[id(2, 0)] * m[id(3, 1)] - m[id(1, 0)] * m[id(2, 3)] * m[id(3, 1)] - m[id(1, 1)] * m[id(2, 0)] * m[id(3, 3)] + m[id(1, 0)] * m[id(2, 1)] * m[id(3, 3)];
	e[id(2, 1)] = m[id(0, 3)] * m[id(2, 1)] * m[id(3, 0)] - m[id(0, 1)] * m[id(2, 3)] * m[id(3, 0)] - m[id(0, 3)] * m[id(2, 0)] * m[id(3, 1)] + m[id(0, 0)] * m[id(2, 3)] * m[id(3, 1)] + m[id(0, 1)] * m[id(2, 0)] * m[id(3, 3)] - m[id(0, 0)] * m[id(2, 1)] * m[id(3, 3)];
	e[id(2, 2)] = m[id(0, 1)] * m[id(1, 3)] * m[id(3, 0)] - m[id(0, 3)] * m[id(1, 1)] * m[id(3, 0)] + m[id(0, 3)] * m[id(1, 0)] * m[id(3, 1)] - m[id(0, 0)] * m[id(1, 3)] * m[id(3, 1)] - m[id(0, 1)] * m[id(1, 0)] * m[id(3, 3)] + m[id(0, 0)] * m[id(1, 1)] * m[id(3, 3)];
	e[id(2, 3)] = m[id(0, 3)] * m[id(1, 1)] * m[id(2, 0)] - m[id(0, 1)] * m[id(1, 3)] * m[id(2, 0)] - m[id(0, 3)] * m[id(1, 0)] * m[id(2, 1)] + m[id(0, 0)] * m[id(1, 3)] * m[id(2, 1)] + m[id(0, 1)] * m[id(1, 0)] * m[id(2, 3)] - m[id(0, 0)] * m[id(1, 1)] * m[id(2, 3)];
	e[id(3, 0)] = m[id(1, 2)] * m[id(2, 1)] * m[id(3, 0)] - m[id(1, 1)] * m[id(2, 2)] * m[id(3, 0)] - m[id(1, 2)] * m[id(2, 0)] * m[id(3, 1)] + m[id(1, 0)] * m[id(2, 2)] * m[id(3, 1)] + m[id(1, 1)] * m[id(2, 0)] * m[id(3, 2)] - m[id(1, 0)] * m[id(2, 1)] * m[id(3, 2)];
	e[id(3, 1)] = m[id(0, 1)] * m[id(2, 2)] * m[id(3, 0)] - m[id(0, 2)] * m[id(2, 1)] * m[id(3, 0)] + m[id(0, 2)] * m[id(2, 0)] * m[id(3, 1)] - m[id(0, 0)] * m[id(2, 2)] * m[id(3, 1)] - m[id(0, 1)] * m[id(2, 0)] * m[id(3, 2)] + m[id(0, 0)] * m[id(2, 1)] * m[id(3, 2)];
	e[id(3, 2)] = m[id(0, 2)] * m[id(1, 1)] * m[id(3, 0)] - m[id(0, 1)] * m[id(1, 2)] * m[id(3, 0)] - m[id(0, 2)] * m[id(1, 0)] * m[id(3, 1)] + m[id(0, 0)] * m[id(1, 2)] * m[id(3, 1)] + m[id(0, 1)] * m[id(1, 0)] * m[id(3, 2)] - m[id(0, 0)] * m[id(1, 1)] * m[id(3, 2)];
	e[id(3, 3)] = m[id(0, 1)] * m[id(1, 2)] * m[id(2, 0)] - m[id(0, 2)] * m[id(1, 1)] * m[id(2, 0)] + m[id(0, 2)] * m[id(1, 0)] * m[id(2, 1)] - m[id(0, 0)] * m[id(1, 2)] * m[id(2, 1)] - m[id(0, 1)] * m[id(1, 0)] * m[id(2, 2)] + m[id(0, 0)] * m[id(1, 1)] * m[id(2, 2)];
	return invDet * (*this);
}



Matrix4 Matrix4::operator+(Matrix4 other)
{
	Matrix4 result = Matrix4();
	for (int i = 0; i < 16; i++) {
		result.elements[i] = elements[i] + other.elements[i];
	}
	return result;
}

Matrix4 Matrix4::operator-(Matrix4 other)
{
	Matrix4 result = Matrix4();
	for (int i = 0; i < 16; i++) {
		result.elements[i] = elements[i] - other.elements[i];
	}
	return result;
}

Matrix4 Matrix4::operator*(Matrix4 other)
{
	Matrix4 result = Matrix4();
	float* re = result.elements;
	float* ae = elements;
	float* be = other.elements;

	float a11 = ae[0], a12 = ae[4], a13 = ae[8], a14 = ae[12];
	float a21 = ae[1], a22 = ae[5], a23 = ae[9], a24 = ae[13];
	float a31 = ae[2], a32 = ae[6], a33 = ae[10], a34 = ae[14];
	float a41 = ae[3], a42 = ae[7], a43 = ae[11], a44 = ae[15];

	float b11 = be[0], b12 = be[4], b13 = be[8], b14 = be[12];
	float b21 = be[1], b22 = be[5], b23 = be[9], b24 = be[13];
	float b31 = be[2], b32 = be[6], b33 = be[10], b34 = be[14];
	float b41 = be[3], b42 = be[7], b43 = be[11], b44 = be[15];

	re[0] = a11 * b11 + a12 * b21 + a13 * b31 + a14 * b41;
	re[4] = a11 * b12 + a12 * b22 + a13 * b32 + a14 * b42;
	re[8] = a11 * b13 + a12 * b23 + a13 * b33 + a14 * b43;
	re[12] = a11 * b14 + a12 * b24 + a13 * b34 + a14 * b44;

	re[1] = a21 * b11 + a22 * b21 + a23 * b31 + a24 * b41;
	re[5] = a21 * b12 + a22 * b22 + a23 * b32 + a24 * b42;
	re[9] = a21 * b13 + a22 * b23 + a23 * b33 + a24 * b43;
	re[13] = a21 * b14 + a22 * b24 + a23 * b34 + a24 * b44;

	re[2] = a31 * b11 + a32 * b21 + a33 * b31 + a34 * b41;
	re[6] = a31 * b12 + a32 * b22 + a33 * b32 + a34 * b42;
	re[10] = a31 * b13 + a32 * b23 + a33 * b33 + a34 * b43;
	re[14] = a31 * b14 + a32 * b24 + a33 * b34 + a34 * b44;

	re[3] = a41 * b11 + a42 * b21 + a43 * b31 + a44 * b41;
	re[7] = a41 * b12 + a42 * b22 + a43 * b32 + a44 * b42;
	re[11] = a41 * b13 + a42 * b23 + a43 * b33 + a44 * b43;
	re[15] = a41 * b14 + a42 * b24 + a43 * b34 + a44 * b44;

	return result;
}

Matrix4 operator*(float scalar, Matrix4 m)
{
	for (int i = 0; i < 16; i++) {
		m.elements[i] *= scalar;
	}
	return m;
}

Matrix4 inverse(Matrix4& m)
{
	Matrix4 result = Matrix4().getInverse(m);
	return result;
}
