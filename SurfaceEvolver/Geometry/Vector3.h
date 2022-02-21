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
	double x, y, z;

	Vector3();
	Vector3(const Vector3& other);
	Vector3(double x, double, double z);
	Vector3(double v[3]);
	~Vector3();

	Vector3 clone();

	void set(double x, double y, double z);
	Vector3 setAndReturn(double x, double y, double z);
	void setCoordById(double val, unsigned int id);
	void min(Vector3 other);
	void max(Vector3 other);
	double getCoordById(unsigned int id);

	bool equals(Vector3 other);
	bool equalsWithEpsilon(Vector3 other, double epsilon);

	void negate();
	double dot(const Vector3& other) const;
	void normalize();
	void lerp(Vector3 other, double param);
	Vector3 cross(Vector3 other);

	double lengthSq() const;
	double length() const;

	void toArray(double* a);

	Vector3& applyQuaternion(Quaternion& q);
	Vector3& applyAxisAngle(Vector3& axis, double angle);
	Vector3& applyMatrix4(Matrix4& m);
	Vector3& applyMatrix3(Matrix3& m);
	Vector3& addScalar(double scalar);
	Vector3& subScalar(double scalar);
	Vector3& multiply(Vector3& other);

	double angleTo(const Vector3& other) const;

	Vector3 operator+ (Vector3& other);
	Vector3 operator- (Vector3& other);
	Vector3 operator* (double scalar);
	Vector3 operator/ (double scalar);
	Vector3 operator+ (const Vector3& other) const;
	Vector3 operator- (const Vector3& other) const;
	friend Vector3 operator* (Matrix3 m, Vector3 a);
	friend Vector3 operator* (Matrix4 m, Vector3 a);
	friend Vector3 operator* (double scalar, Vector3 a);
	friend Vector3 operator/ (double scalar, Vector3 a);
	Vector3& operator+=(const Vector3& other);
	Vector3& operator-=(const Vector3& other);
	Vector3& operator*=(const double& scalar);
	Vector3& operator/=(const double& scalar);
	friend std::ostream& operator<< (std::ostream& out, const Vector3& v);

	friend bool operator< (const Vector3& left, const Vector3& right);
};

Vector3 normalize(Vector3 target);
double dot(Vector3 a, Vector3 b);
Vector3 cross(Vector3 a, Vector3 b);
Vector3 lerp(Vector3 v1, Vector3 v2, double param);
Vector3 multiply(Vector3 a, Vector3 b);
bool equal(Vector3& a, Vector3& b);
bool notEqual(Vector3& a, Vector3& b);

#endif

