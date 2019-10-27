
#include <iostream>
#include "Geometry.h"
#include "Matrix4.h"
#include "Matrix3.h"
#include "Vector3.h"
#include "Icosphere.h"
#include "PrimitiveBox.h"
#include "CubeSphere.h"
#include "VTKExporter.h"
#include "OBJImporter.h"
#include "AABBTree.h"

//   TODO:
// - Add an AABBTree structure
// - Make a fast cell intersection query
// - Set intersected cell values to 0 and INFINITY everywhere else
// - Apply Fast Sweeping Method
// - Alternatively: Make a fast distance query (CUDA?)
// - Interior/Exterior sign
// - SDF

int main()
{
	float r = 50.0f;
	unsigned int d = 3;
	IcoSphere ico = IcoSphere(d, r);
	float a = 2 * r / sqrt(3.);
	unsigned int ns = 10;
	PrimitiveBox box = PrimitiveBox(a, a, a, ns, ns, ns);
	CubeSphere cs = CubeSphere(ns, r);

	VTKExporter e = VTKExporter();
	e.initExport(ico, "icosphere");
	e.initExport(box, "box");

	box.applyMatrix(Matrix4().makeTranslation(-a / 2., -a / 2., -a / 2.));
	e.initExport(box, "boxTranslated");
	e.initExport(cs, "cubesphere");

	/*
	OBJImporter im = OBJImporter();
	Geometry cube = im.importOBJGeometry("Cube.obj");
	e.initExport(cube, "cube");
	Geometry polySphere = im.importOBJGeometry("PolygonalSphere.obj");
	e.initExport(polySphere, "polySphere");

	PrimitiveBox cubeCorrect = PrimitiveBox(100., 100., 100., 1, 1, 1);
	cubeCorrect.applyMatrix(Matrix4().makeTranslation(-50., -50., 0.));

	e.initExport(cubeCorrect, "cubeCorrect");

	Geometry unstructCube = im.importOBJGeometry("unstructCube.obj");
	e.initExport(unstructCube, "unstructCube");
	*/

	std::vector<StructGeom::Triangle> t = cs.getTriangles();
	// AABBTree<StructGeom::Triangle> tree = AABBTree<StructGeom::Triangle>(t, 100);

	return 1;
}
