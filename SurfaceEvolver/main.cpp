
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
	unsigned int ns = 10;
	PrimitiveBox box = PrimitiveBox(a, a, a, ns, ns, ns);
	CubeSphere cs = CubeSphere(ns, r);

	VTKExporter e = VTKExporter();
	e.initExport(ico, "icosphere");
	e.initExport(box, "box");

	box.applyMatrix(Matrix4().makeTranslation(-a / 2., -a / 2., -a / 2.));
	e.initExport(box, "boxTranslated");
	e.initExport(cs, "cubesphere");

	OBJImporter obj = OBJImporter();
	Geometry bunny = obj.importOBJGeometry("bunny.obj");
	e.initExport(bunny, "sfBunny");

	std::vector<Tri> triangs = bunny.getTriangles();
	AABBTree T = AABBTree(triangs, bunny.getBoundingBox());
	Octree O = Octree(&T, T.bbox, 2.5f);
	std::cout << "octree construction finished" << std::endl;

	/*
	std::vector<Tri> triangs = cs.getTriangles();
	AABBTree T = AABBTree(triangs, cs.getBoundingBox());
	Octree O = Octree(&T, T.bbox, 2.5f);
	std::cout << "octree construction finished" << std::endl;*/

	/*
	Box3 csBBox = cs.getBoundingBox();
	Vector3 csBBoxSize = csBBox.getSize();
	std::vector<Geometry> lboxes = {};
	uint N = 20;
	for (uint i = 0; i < N; i++) {
		for (uint j = 0; j < N; j++) {
			for (uint k = 0; k < N; k++) {
				Vector3 offset0 = multiply(Vector3((float)i / (float)N, (float)j / (float)N, (float)k / (float)N), csBBoxSize);
				Vector3 offset1 = multiply(Vector3((float)(i + 1) / (float)N, (float)(j + 1) / (float)N, (float)(k + 1) / (float)N), csBBoxSize);
				Box3 b = Box3(csBBox.min + offset0, csBBox.min + offset1);
				std::vector<Tri> tris = T.getTrianglesInABox(b);
				if (tris.size()) {
					float dimX = b.max.x - b.min.x;
					float dimY = b.max.y - b.min.y;
					float dimZ = b.max.z - b.min.z;
					PrimitiveBox box = PrimitiveBox(dimX, dimY, dimZ, 1, 1, 1);
					Vector3 t = b.min;
					box.applyMatrix(Matrix4().makeTranslation(t.x, t.y, t.z));
					lboxes.push_back(box);
				}
			}
		}
	}
	Geometry resultLboxGeom = mergeGeometries(lboxes);
	e.initExport(resultLboxGeom, "leafBoxesOctree");*/

	std::vector<Geometry> otlBoxGeoms = O.getLeafBoxGeoms();
	Geometry leafBoxGeom = mergeGeometries(otlBoxGeoms);
	e.initExport(leafBoxGeom, "leafBoxesOctree");

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

