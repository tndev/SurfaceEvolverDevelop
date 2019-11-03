
#include <iostream>
#include <chrono>
#include "../Geometry/Geometry.h"
#include "../Geometry/Matrix4.h"
#include "../Geometry/Matrix3.h"
#include "../Geometry/Vector3.h"
#include "../GeometryObject/Icosphere.h"
#include "../GeometryObject/PrimitiveBox.h"
#include "../GeometryObject/CubeSphere.h"
#include "../ExportImport/VTKExporter.h"
#include "../ExportImport/OBJImporter.h"
#include "../BVH/AABBTree.h"
#include "../BVH/Octree.h"
#include "../SDF/Grid.h"
#include "../SDF/FastSweep3D.h"

//   DONE:
//
// - Add an AABBTree structure
// - Add an Octree for activating intersected grid cells
// - Make a "fast" cell intersection query
//

//   WIP:
//
// - Set intersected cell values to 0 and INFINITY everywhere else (WIP)

//   TODO:
//
// - Apply Fast Sweeping Method
// - Alternatively: Make a fast distance query (CUDA?)
// - Interior/Exterior sign
// - SDF
// - adaptive resampling for AABB Tree construction
//

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


	auto startObjLoad = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	OBJImporter obj = OBJImporter();
	Geometry bunny = obj.importOBJGeometry("bunny.obj");
	// bunny.applyMatrix(Matrix4().setToScale(100.0f, 100.0f, 100.0f));
	e.initExport(bunny, "sfBunny");
	// === Timed code ============
	auto endObjLoad = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedObj = (endObjLoad - startObjLoad);
	std::cout << "Model loaded after " << elapsedObj.count() << " seconds" << std::endl;


	auto startAABBtree = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	std::cout << "fetching triangles..." << std::endl;
	std::vector<Tri> triangs = bunny.getTriangles();
	std::cout << "initializing AABB construction for " << triangs.size() << " triangles..." << std::endl;
	AABBTree T = AABBTree(triangs, bunny.getBoundingBox());
	// === Timed code ============
	auto endAABBtree = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedAABB = (endAABBtree - startAABBtree);
	std::cout << "AABBTree construction finished after " << elapsedAABB.count() << " seconds" << std::endl;


	auto startOctree = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	uint res = 60; // octree resolution
	std::cout << "initializing Octree construction for " << triangs.size() << " triangles with resolution " << res << std::endl;
	Octree O = Octree(&T, T.bbox, res);
	// === Timed code ============
	auto endOctree = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedOctree = (endOctree - startOctree);
	std::cout << "Octree construction finished after " << elapsedOctree.count() << " seconds" << std::endl;

	bool OctreeLeafBoxes = false;

	if (OctreeLeafBoxes) {
		auto startOctreeBoxes = std::chrono::high_resolution_clock::now();
		// === Timed code ============
		std::cout << "Exporting octree leaf voxels..." << std::endl;
		std::vector<Geometry> otlBoxGeoms = {};
		O.getLeafBoxGeoms(&otlBoxGeoms);
		Geometry leafBoxGeom = mergeGeometries(otlBoxGeoms);
		e.initExport(leafBoxGeom, "leafBoxesOctree"); // this is a major bottleneck
		// === Timed code ============
		auto endOctreeBoxes = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedOctreeBox = (endOctreeBoxes - startOctreeBoxes);
		std::cout << "Octree voxel export finished after " << elapsedOctreeBox.count() << " seconds" << std::endl;
	}

	bool OctreeLeafGrid = true;

	if (OctreeLeafGrid) {
		auto startOctreeGrid = std::chrono::high_resolution_clock::now();
		// === Timed code ============
		std::cout << "Exporting octree leaf voxels into grid..." << std::endl;
		Grid voxField = Grid(res, res, res, O.cubeBox);
		O.setLeafValueToScalarGrid(&voxField, 0.0f, true);
		voxField.exportToVTI("voxField"); // save initial cond
		Grid g;
		for (uint s = 0; s < 9; s++) {
			g = voxField;
			FastSweep3D fs = FastSweep3D(&g, s, s == 8); // computes distance field
			g.exportToVTI("voxField" + std::to_string(s));
		}
		voxField = g;
		voxField.exportToVTI("voxFieldSDF"); // save final SDF
		// === Timed code ============
		auto endOctreeGrid = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedOctreeField = (endOctreeGrid - startOctreeGrid);
		std::cout << "Octree voxel export finished after " << elapsedOctreeField.count() << " seconds" << std::endl;
	}

	auto startAABBdepth = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	std::cout << "calculating AABB depth..." << std::endl;
	unsigned int maxDepth = depth(&T);
	// === Timed code ============
	auto endAABBdepth = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedAABBdepth = (endAABBdepth - startAABBdepth);
	std::cout << "AABB depth calculation finished after " << elapsedAABBdepth.count() << " seconds" << std::endl;

	bool AABBleafExport = false;

	if (AABBleafExport) {
		auto startAABBleafExport = std::chrono::high_resolution_clock::now();
		// === Timed code ============
		std::cout << "exporting AABB leaf geoms of depth = 0, ... " << maxDepth << "." << std::endl;
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
		// === Timed code ============
		auto endAABBleafExport = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedLeafExport = (endAABBleafExport - startAABBleafExport);
		std::cout << "AABB voxel export finished after " << elapsedLeafExport.count() << " seconds" << std::endl;
	}

	return 1;
}

