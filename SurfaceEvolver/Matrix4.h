#ifndef MATRIX4_H_
#define MATRIX4_H_

#include <iostream>

class Matrix4
{
public:
	float elements[16] = {
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};

	float determinant();
	Matrix4 getInverse(Matrix4& from);

	Matrix4 operator+ (Matrix4 other);
	Matrix4 operator- (Matrix4 other);
	Matrix4 operator* (Matrix4 other);
	friend Matrix4 operator*(float scalar, Matrix4 m);
};

Matrix4 inverse(Matrix4& m);

#endif;