
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
#include "Octree.h"

//   TODO:
// - Add an AABBTree structure (done)
// - Add an Octree for activating intersected grid cells
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
	unsigned int ns = 2;
	PrimitiveBox box = PrimitiveBox(a, a, a, ns, ns, ns);
	CubeSphere cs = CubeSphere(ns, r);

	VTKExporter e = VTKExporter();
	e.initExport(ico, "icosphere");
	e.initExport(box, "box");

	box.applyMatrix(Matrix4().makeTranslation(-a / 2., -a / 2., -a / 2.));
	e.initExport(box, "boxTranslated");
	e.initExport(cs, "cubesphere");

	std::vector<Tri> triangs = cs.getTriangles();
	AABBTree T = AABBTree(triangs, cs.getBoundingBox());
	Octree O = Octree(&T, T.bbox, 10.0f);
	std::cout << "octree construction finished" << std::endl;

	std::vector<Geometry> lboxes = O.getLeafBoxGeoms();
	Geometry resultLboxGeom = mergeGeometries(lboxes);
	e.initExport(resultLboxGeom, "leafBoxesOctree");

	unsigned int maxDepth = depth(&T);

	for (unsigned int d = 0; d < maxDepth; d++) {
		std::vector<Geometry> boxes = T.getAABBGeomsOfDepth(d);
		std::vector<Geometry> triangles = T.getAABBTrianglesOfDepth(d);
		Geometry resultGeom = mergeGeometries(boxes);
		Geometry resultTriGeom = mergeGeometries(triangles);
		e.initExport(resultGeom, "boxes" + std::to_string(d) + "AABB");
		std::cout << "boxes" << d << "AABB saved" << std::endl;
		e.initExport(resultTriGeom, "triangles" + std::to_string(d) + "AABB");
		std::cout << "triangles" << d << "AABB saved" << std::endl;
	}

	std::vector<Geometry> leafBoxes = T.getAABBLeafGeoms();
	Geometry leafBoxesGeom = mergeGeometries(leafBoxes);
	e.initExport(leafBoxesGeom, "leafBoxesAABB");

	return 1;
}

