
#include "Geometry.h"
#include "Vector3.h"
#include <iostream>

int main()
{
	// Geometry geom();
	Vector3 a(1., 1., 2.);
	Vector3 b(-1., 3., -1.);
	Vector3 c = a + b;
	Vector3 d = a - 2 * b;
    std::cout << "d = " << d;
}
