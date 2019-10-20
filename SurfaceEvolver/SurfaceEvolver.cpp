
#include "Geometry.h"
#include "Vector3.h"
#include "Icosphere.h"
#include <iostream>

int main()
{
	Vector3 a(1., 1., 2.);
	Vector3 b(-1., 3., -1.);
	Vector3 c = a + b;
	Vector3 d = a - 2 * b;
    std::cout << "d = " << normalize(d);

	Icosphere ico = Icosphere(0, 2);
}
