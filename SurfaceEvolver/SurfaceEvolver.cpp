
#include <iostream>
#include "Geometry.h"
#include "Matrix4.h"
#include "Vector3.h"
#include "Icosphere.h"

int main()
{
	Vector3 a(1., 1., 2.);
	Vector3 b(-1., 3., -1.);
	Vector3 c = a + b;
	Vector3 d = a - 2 * b;
    std::cout << "d = " << normalize(d);

	Matrix4 m = Matrix4();
	Icosphere ico = Icosphere(0, 2);
}
