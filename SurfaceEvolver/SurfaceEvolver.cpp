
#include <iostream>
#include "Geometry.h"
#include "Matrix4.h"
#include "Matrix3.h"
#include "Vector3.h"
#include "Icosphere.h"
#include "PrimitiveBox.h"
#include "CubeSphere.h"
#include "VTKExporter.h"

int main()
{
	IcoSphere ico = IcoSphere(1, 50.);
	PrimitiveBox box = PrimitiveBox(53., 52., 51., 4, 3, 3);

	VTKExporter e = VTKExporter();
	e.initExport(ico, "icosphere");
	e.initExport(box, "box");

	return 1;
}
