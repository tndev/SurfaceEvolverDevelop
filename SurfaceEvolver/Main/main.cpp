#define _USE_MATH_DEFINES

#include <iostream>
#include <chrono>
#include <math.h>
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
#include "../SDF/SDF.h"

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
// - Unite AABB, Octree and FastSweep3D into a single class
// - debug and optimize FastSweep3D

//  POSTPONED:
//
// - Interior/Exterior mesh Sign (Needs to build on top of FastSweep3D)


//   DONE, BUT MIGHT BE IMPROVED:
//
// - AABB and Octree have to take as little space as possible

//   WIP:
// 
// - cleanup main & prep for VTK window form


//   TODO (weekend Nov.15th - Nov.17th):
//
// - flat AABB
// - perform simple DF tests for geom primitives like sphere, icosphere, cubesphere

int main()
{
	float r = 50.0f;
	unsigned int d = 3;
	float a = 2 * r / sqrt(3.);
	unsigned int ns = 10;

	VTKExporter e = VTKExporter();

	bool iterateCubeSphereTest = true;

	if (iterateCubeSphereTest) {
		size_t min_Res = 30, max_Res = 40;
		size_t min_Ns = 0, max_Ns = 1;
		std::fstream timing("timing.txt", std::fstream::out);

		for (unsigned int res = min_Res; res < max_Res; res += 20) {
			for (unsigned int n = min_Ns; n < max_Ns; n++) {
				std::cout << "icosphere(" << n << "), grid_res = " << res << std::endl;
				IcoSphere g = IcoSphere(n, r);
				Vector3 axis = normalize(Vector3(1, 1, 1));
				g.applyMatrix(Matrix4().makeRotationAxis(axis.x, axis.y, axis.z, M_PI / 6.));
				e.initExport(g, "icosphere" + std::to_string(res) + "-" + std::to_string(n));
				
				std::cout << "init SDF..." << std::endl;

				SDF sdf_FS = SDF(&g, res);

				std::cout << sdf_FS.getComputationProperties();
				timing << sdf_FS.getComputationProperties();
				
				sdf_FS.exportGrid(&e); // save to vti	

				/*
				SDF sdf_AABB = SDF(&g, res, SDF_Method::aabb_dist);

				std::cout << sdf_AABB.getComputationProperties();
				timing << sdf_AABB.getComputationProperties();

				sdf_AABB.exportGrid(&e);

				SDF sdf_Brute = SDF(&g, res, SDF_Method::brute_force);

				std::cout << sdf_Brute.getComputationProperties();
				timing << sdf_Brute.getComputationProperties();

				sdf_Brute.exportGrid(&e);

				Grid FSerror = absGrid(subGrids(*sdf_FS.grid, *sdf_Brute.grid));
				Grid AABBerror = absGrid(subGrids(*sdf_AABB.grid, *sdf_Brute.grid));

				FSerror.exportToVTI("voxField_" + std::to_string(res) + "FS_ERROR");
				AABBerror.exportToVTI("voxField_" + std::to_string(res) + "AABB_ERROR");
				*/
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


	uint res = 30; // octree resolution	
	SDF bunny_sdf = SDF(&bunny, res);

	std::cout << bunny_sdf.getComputationProperties();

	bunny_sdf.exportGrid(&e, "bunnySDF");

	return 1;
}

