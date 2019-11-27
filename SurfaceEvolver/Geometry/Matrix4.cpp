#include "Matrix4.h"

#define id(i, j) 4 * i + j
#define singularity 0

Matrix4::Matrix4()
{
}

Matrix4::Matrix4(const Matrix4& other)
{
	const float* e = other.elements;
	this->set(
		e[0], e[1], e[2], e[3],
		e[4], e[5], e[6], e[7],
		e[8], e[9], e[10], e[11],
		e[12], e[13], e[14], e[15]
	);
}

Matrix4::Matrix4(float* elems)
{
	float* e = elements;
	e[0] = elems[0];	e[1] = elems[1];	e[2] = elems[2];	e[3] = elems[3];
	e[4] = elems[4];	e[5] = elems[5];	e[6] = elems[6];	e[7] = elems[7];
	e[8] = elems[8];	e[9] = elems[9];	e[10] = elems[10];	e[11] = elems[11];
	e[12] = elems[12];	e[13] = elems[13];	e[14] = elems[14];	e[15] = elems[15];
}

Matrix4::Matrix4(float m00, float m01, float m02, float m03, float m10, float m11, float m12, float m13, float m20, float m21, float m22, float m23, float m30, float m31, float m32, float m33)
{
	this->set(
		m00, m01, m02, m03,
		m10, m11, m12, m13,
		m20, m21, m22, m23,
		m30, m31, m32, m33
	);
}

Matrix4::~Matrix4()
{
}

void Matrix4::set(float m00, float m01, float m02, float m03, float m10, float m11, float m12, float m13, float m20, float m21, float m22, float m23, float m30, float m31, float m32, float m33)
{
	float* e = elements;
	e[0] = m00;		e[1] = m01;		e[2] = m02;		e[3] = m03;
	e[4] = m10;		e[5] = m11;		e[6] = m12;		e[7] = m13;
	e[8] = m20;		e[9] = m21;		e[10] = m22;	e[11] = m23;
	e[12] = m30;	e[13] = m31;	e[14] = m32;	e[15] = m33;
}

Matrix3 Matrix4::getSubMatrix3()
{
	float* e = elements;
	float resultElems[9] = {
		e[0],	e[4],	e[8],
		e[1],	e[5],	e[9],
		e[2],	e[6],	e[10]
	};
	return Matrix3(resultElems);
}

bool Matrix4::isIdentity()
{
	float* e = elements;
	return (
		e[0] == 1.0f  &&   e[1] == 0.0f  &&   e[2] == 0.0f  &&   e[3] == 0.0f &&
		e[4] == 0.0f  &&   e[5] == 1.0f  &&   e[6] == 0.0f  &&   e[7] == 0.0f &&
		e[8] == 0.0f  &&   e[9] == 0.0f  &&  e[10] == 1.0f  &&  e[11] == 0.0f &&
		e[12] == 0.0f &&  e[13] == 0.0f  &&  e[14] == 0.0f  &&  e[15] == 1.0f
	);
}

void Matrix4::setToIdentity()
{
	this->set(
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);
}

void Matrix4::compose(Vector3* position, Quaternion* quaternion, Vector3* scale)
{
	float* te = this->elements;

	float x = quaternion->x, y = quaternion->y, z = quaternion->z, w = quaternion->w;
	float x2 = x + x, y2 = y + y, z2 = z + z;
	float xx = x * x2, xy = x * y2, xz = x * z2;
	float yy = y * y2, yz = y * z2, zz = z * z2;
	float wx = w * x2, wy = w * y2, wz = w * z2;

	te[0] = (1 - (yy + zz)) * scale->x;
	te[1] = (xy + wz) * scale->x;
	te[2] = (xz - wy) * scale->x;
	te[3] = 0;

	te[4] = (xy - wz) * scale->y;
	te[5] = (1 - (xx + zz)) * scale->y;
	te[6] = (yz + wx) * scale->y;
	te[7] = 0;

	te[8] = (xz + wy) * scale->z;
	te[9] = (yz - wx) * scale->z;
	te[10] = (1 - (xx + yy)) * scale->z;
	te[11] = 0;

	te[12] = position->x;
	te[13] = position->y;
	te[14] = position->z;
	te[15] = 1;
}

void Matrix4::transpose()
{
	float* e = elements;
	this->set(
		e[0], e[4], e[8], e[12],
		e[1], e[5], e[9], e[13],
		e[2], e[6], e[10], e[14],
		e[3], e[7], e[11], e[15]
	);
}

void Matrix4::multiplyScalar(float scalar)
{
	float* e = this->elements;
	for (int i = 0; i < 16; i++) {
		e[i] *= scalar;
	}
}

void Matrix4::multiplyMatrices(Matrix4& a, Matrix4& b)
{
	float* ae = a.elements;
	float* be = b.elements;
	float* te = this->elements;

	float a11 = ae[0], a12 = ae[4], a13 = ae[8], a14 = ae[12];
	float a21 = ae[1], a22 = ae[5], a23 = ae[9], a24 = ae[13];
	float a31 = ae[2], a32 = ae[6], a33 = ae[10], a34 = ae[14];
	float a41 = ae[3], a42 = ae[7], a43 = ae[11], a44 = ae[15];

	float b11 = be[0], b12 = be[4], b13 = be[8], b14 = be[12];
	float b21 = be[1], b22 = be[5], b23 = be[9], b24 = be[13];
	float b31 = be[2], b32 = be[6], b33 = be[10], b34 = be[14];
	float b41 = be[3], b42 = be[7], b43 = be[11], b44 = be[15];

	te[0] = a11 * b11 + a12 * b21 + a13 * b31 + a14 * b41;
	te[4] = a11 * b12 + a12 * b22 + a13 * b32 + a14 * b42;
	te[8] = a11 * b13 + a12 * b23 + a13 * b33 + a14 * b43;
	te[12] = a11 * b14 + a12 * b24 + a13 * b34 + a14 * b44;

	te[1] = a21 * b11 + a22 * b21 + a23 * b31 + a24 * b41;
	te[5] = a21 * b12 + a22 * b22 + a23 * b32 + a24 * b42;
	te[9] = a21 * b13 + a22 * b23 + a23 * b33 + a24 * b43;
	te[13] = a21 * b14 + a22 * b24 + a23 * b34 + a24 * b44;

	te[2] = a31 * b11 + a32 * b21 + a33 * b31 + a34 * b41;
	te[6] = a31 * b12 + a32 * b22 + a33 * b32 + a34 * b42;
	te[10] = a31 * b13 + a32 * b23 + a33 * b33 + a34 * b43;
	te[14] = a31 * b14 + a32 * b24 + a33 * b34 + a34 * b44;

	te[3] = a41 * b11 + a42 * b21 + a43 * b31 + a44 * b41;
	te[7] = a41 * b12 + a42 * b22 + a43 * b32 + a44 * b42;
	te[11] = a41 * b13 + a42 * b23 + a43 * b33 + a44 * b43;
	te[15] = a41 * b14 + a42 * b24 + a43 * b34 + a44 * b44;
}

void Matrix4::multiply(Matrix4& m)
{
	this->multiplyMatrices(*this, m);
}

Matrix4 Matrix4::multiply(Matrix4 m)
{
	Matrix4 result = Matrix4();
	result.multiplyMatrices(*this, m);
	return result;
}

void Matrix4::premultiply(Matrix4& m)
{
	this->multiplyMatrices(m, *this);
}

Matrix4 Matrix4::premultiply(Matrix4 m)
{
	Matrix4 result = Matrix4();
	result.multiplyMatrices(m, *this);
	return result;
}

Matrix4 Matrix4::setToScale(float sx, float sy, float sz)
{
	this->set(
		sx,   0.0f, 0.0f, 0.0f,
		0.0f, sy,   0.0f, 0.0f,
		0.0f, 0.0f, sz,   0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);

	return *this;
}

Matrix4 Matrix4::makeRotationAxis(float ax, float ay, float az, float angle)
{
	float c = cos(angle), s = sin(angle);
	float t = 1 - c;
	float tx = t * ax, ty = t * ay;

	this->set(
		tx * ax + c, tx * ay - s * az, tx * az + s * ay, 0.0f,
		tx * ay + s * az, ty * ay + c, ty * az - s * ax, 0.0f,
		tx * az - s * ay, ty * az + s * ax, t * az * az + c, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);

	return *this;
}

Matrix4 Matrix4::makeTranslation(float tx, float ty, float tz)
{
	this->set(
		1.0f, 0.0f, 0.0f, tx,
		0.0f, 1.0f, 0.0f, ty,
		0.0f, 0.0f, 1.0f, tz,
		0.0f, 0.0f, 0.0f, 1.0f
	);

	return *this;
}

void Matrix4::decompose(Vector3* position, Quaternion* quaternion, Vector3* scale)
{
	float* te = this->elements;

	Vector3 _v1 = Vector3();
	float sx = _v1.setAndReturn(te[0], te[1], te[2]).length();
	float sy = _v1.setAndReturn(te[4], te[5], te[6]).length();
	float sz = _v1.setAndReturn(te[8], te[9], te[10]).length();

	// if determine is negative, we need to invert one scale
	float det = this->determinant();
	if (det < 0) sx = -sx;

	position->x = te[12];
	position->y = te[13];
	position->z = te[14];

	// scale the rotation part
	Matrix4 _m1 = *this;

	float invSX = 1.0f / sx;
	float invSY = 1.0f / sy;
	float invSZ = 1.0f / sz;

	_m1.elements[0] *= invSX;
	_m1.elements[1] *= invSX;
	_m1.elements[2] *= invSX;

	_m1.elements[4] *= invSY;
	_m1.elements[5] *= invSY;
	_m1.elements[6] *= invSY;

	_m1.elements[8] *= invSZ;
	_m1.elements[9] *= invSZ;
	_m1.elements[10] *= invSZ;

	quaternion->setFromRotationMatrix(&_m1);

	scale->x = sx;
	scale->y = sy;
	scale->z = sz;
}

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
	return value;
}

Matrix4 Matrix4::getInverse(Matrix4& from)
{
	float* m = from.elements;
	float* e = elements;
	float det = from.determinant();
	Matrix4 result = Matrix4();
	int ex = singularity;
	try {
		if (fabs(det) < FLT_MIN) {
			throw singularity;
		} else {
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
			result = *this;
			return invDet * result;
		}
	}
	catch (int ex) {
		std::cout << "det(A) = 0.0! Unable to invert a singular matrix!" << std::endl;
	}

	return result;
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
	Matrix4 result = m;
	result.multiplyScalar(scalar);
	return result;
}

Matrix4 inverse(Matrix4& m)
{
	Matrix4 result = Matrix4().getInverse(m);
	return result;
}

Matrix4 transpose(Matrix4& m)
{
	Matrix4 result = m;
	result.transpose();
	return result;
}
