
#include <iostream>
#include "Geometry.h"
#include "Matrix4.h"
#include "Matrix3.h"
#include "Vector3.h"
#include "Icosphere.h"

int main()
{
	Vector3 a(1., 1., 2.);
	Vector3 b(-1., 3., -1.);
	Vector3 c = a + b;
	Vector3 d = a - 2 * b;
	float melems[16] = {
		2.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 2.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 2.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	Matrix4 m = Matrix4(melems);
	m.multiplyScalar(3.0f);
	Vector3 e = d.clone();
	e.normalize();
	e.applyMatrix4(m);
	std::cout << "M_00 = " << m.elements[0] << std::endl;
	std::cout << "normalize(d) = " << normalize(d) << std::endl;
    std::cout << "M * normalize(d) = " << e << std::endl;

	Icosphere ico = Icosphere(0, 2);
}
