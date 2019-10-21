
#include <iostream>
#include "Geometry.h"
#include "Matrix4.h"
#include "Matrix3.h"
#include "Vector3.h"
#include "Icosphere.h"
#include "VTKExporter.h"

int main()
{
	Icosphere ico = Icosphere(6, 50.);

	VTKExporter e = VTKExporter();
	e.initExport(ico, "icosphere");

	return 1;
}
