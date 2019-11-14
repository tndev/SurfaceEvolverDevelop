
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
// - Set intersected cell values to 0 and INFINITY everywhere else
// - Apply Fast Sweeping Method
// - SDF
// - Optimize initial condition for FastSweep3D by actual distance computation for Octree leaves
// - Get exact octree centroid to closest triangle dist

//   DONE, BUT MIGHT BE IMPROVED:
//
// - AABB and Octree have to take as little space as possible

//   WIP:
// 


//   TODO (weekend Nov.15th - Nov.17th):
//
// - generalize AABB for points, edges, and triangles as generic sufrace features
// - Unite AABB, Octree and FastSweep3D into a single class
// - flat AABB
// - debug and optimize FastSweep3D
// - perform simple DF tests for geom primitives like sphere, icosphere, cubesphere
// - cleanup main & prep for VTK window form
// - Interior/Exterior mesh Sign

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

	bool iterateCubeSphereTest = true;

	if (iterateCubeSphereTest) {
		size_t min_Res = 40, max_Res = 50;
		size_t min_Ns = 3, max_Ns = 4;
		std::fstream timing("timing.txt", std::fstream::out);

		for (unsigned int res = min_Res; res < max_Res; res += 20) {
			for (unsigned int n = min_Ns; n < max_Ns; n++) {
				std::cout << "cubeSphere(" << n << "), grid_res = " << res << std::endl;
				timing << "cubeSphere(" << n << "), grid_res = " << res << std::endl;
				CubeSphere c = CubeSphere(n, r);
				e.initExport(c, "cubeSphere" + std::to_string(res) + "-" + std::to_string(n));
				auto startSDF = std::chrono::high_resolution_clock::now();
				// === Timed code ============

				auto startSDF_AABB = std::chrono::high_resolution_clock::now();
				AABBTree cT = AABBTree(&c);
				auto endSDF_AABB = std::chrono::high_resolution_clock::now();
				std::chrono::duration<float> elapsedSDF_AABB = (endSDF_AABB - startSDF_AABB);

				auto startSDF_Octree = std::chrono::high_resolution_clock::now();
				Octree O = Octree(&cT, cT.bbox, res);
				auto endSDF_Octree = std::chrono::high_resolution_clock::now();
				std::chrono::duration<float> elapsedSDF_Octree = (endSDF_Octree - startSDF_Octree);

				auto startSDF_FS = std::chrono::high_resolution_clock::now();
				Grid voxField_SDF = Grid(res, res, res, O.cubeBox);
				O.setLeafValueToScalarGrid(&voxField_SDF);
				FastSweep3D fs = FastSweep3D(&voxField_SDF, 8); // computes distance field
				auto endSDF_FS = std::chrono::high_resolution_clock::now();
				std::chrono::duration<float> elapsedSDF_FS = (endSDF_FS - startSDF_FS);
				// auto startSDF_Sign = std::chrono::high_resolution_clock::now();
				// voxField_SDF.computeSignField(&cT);
				// auto endSDF_Sign = std::chrono::high_resolution_clock::now();
				// std::chrono::duration<float> elapsedSDF_Sign = (endSDF_Sign - startSDF_Sign);

				auto endSDF = std::chrono::high_resolution_clock::now();
				std::chrono::duration<float> elapsedSDF = (endSDF - startSDF);

				std::cout << "computation times:  AABB: " << elapsedSDF_AABB.count() <<
					" s , Octree: " << elapsedSDF_Octree.count() << " s, FastSweep3D: " << elapsedSDF_FS.count() <<
					// ", Sign: " << elapsedSDF_Sign.count() <<
					", TOTAL: " << elapsedSDF.count() << " s" << std::endl;
				timing << "computation times:  AABB: " << elapsedSDF_AABB.count() <<
					" s , Octree: " << elapsedSDF_Octree.count() << " s, FastSweep3D: " << elapsedSDF_FS.count() <<
					// ", Sign: " << elapsedSDF_Sign.count() <<
					", TOTAL: " << elapsedSDF.count() << " s" << std::endl;
			
				// === Timed code ============
				voxField_SDF.exportToVTI("voxFieldSDF" + std::to_string(res) + "-" + std::to_string(n)); // save to vti		
			}
		}
		timing.close();
	}


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
	std::cout << "initializing AABB construction ..." << std::endl;
	AABBTree T = AABBTree(&bunny);
	// === Timed code ============
	auto endAABBtree = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedAABB = (endAABBtree - startAABBtree);
	std::cout << "AABBTree construction finished after " << elapsedAABB.count() << " seconds" << std::endl;


	auto startOctree = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	uint res = 30; // octree resolution
	std::cout << "initializing Octree construction for " << T.triangles.size() << " triangles with resolution " << res << std::endl;

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
		O.setLeafValueToScalarGrid(&voxField);
		voxField.exportToVTI("voxField"); // save initial cond

		FastSweep3D fs = FastSweep3D(&voxField, 8); // computes distance field

		bool signCompute = false;
		if (signCompute) {
			std::cout << "computing signs for " << voxField.field.size() << " grid pts ..." << std::endl;
			auto startSDF_Sign = std::chrono::high_resolution_clock::now();
			voxField.computeSignField(&T);
			auto endSDF_Sign = std::chrono::high_resolution_clock::now();
			std::chrono::duration<float> elapsedSDF_Sign = (endSDF_Sign - startSDF_Sign);
			std::cout << "sign computation finished after " << elapsedSDF_Sign.count() << " seconds" << std::endl;
		}

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

