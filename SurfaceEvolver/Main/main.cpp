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
// - perform simple DF tests for geom primitives like sphere, icosphere, cubesphere


//   TODO (weekend Nov.15th - Nov.17th):
//
// - cleanup main & prep for VTK window form
// - flat AABB
// - Inverse transform grid upon transforming mesh
// - interpolate distance field for higher resolutions
// - debug AABB closest primitive lookup

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

				// Fast sweeping DF, resized from 20 and interpolated
				SDF sdf_FS_r = SDF(&g, res, true, SDF_Method::fast_sweeping);

				std::cout << sdf_FS_r.getComputationProperties();
				timing << sdf_FS_r.getComputationProperties();
				
				sdf_FS_r.exportGrid(&e); // save to vti

				// Fast sweeping DF
				SDF sdf_FS = SDF(&g, res);

				std::cout << sdf_FS.getComputationProperties();
				timing << sdf_FS.getComputationProperties();

				sdf_FS.exportGrid(&e);

				// AABB DF
				SDF sdf_AABB = SDF(&g, res, false, SDF_Method::aabb_dist);

				std::cout << sdf_AABB.getComputationProperties();
				timing << sdf_AABB.getComputationProperties();

				sdf_AABB.exportGrid(&e);

				// Brute force DF
				SDF sdf_Brute = SDF(&g, res, false, SDF_Method::brute_force);

				std::cout << sdf_Brute.getComputationProperties();
				timing << sdf_Brute.getComputationProperties();

				sdf_Brute.exportGrid(&e);

				Grid FSerror_r = absGrid(subGrids(*sdf_FS_r.grid, *sdf_Brute.grid));
				Grid FSerror = absGrid(subGrids(*sdf_FS.grid, *sdf_Brute.grid));
				Grid AABBerror = absGrid(subGrids(*sdf_AABB.grid, *sdf_Brute.grid));

				float error = FSerror_r.getL2Norm();
				std::cout << "FS_ERROR_resized L2 Norm: " << error << std::endl;
				timing << "FS_ERROR_resized L2 Norm: " << error << std::endl;

				error = FSerror.getL2Norm();
				std::cout << "FS_ERROR L2 Norm: " << error << std::endl;
				timing << "FS_ERROR L2 Norm: " << error << std::endl;

				error = AABBerror.getL2Norm();
				std::cout << "AABB_ERROR L2 Norm: " << error << std::endl;
				timing << "AABB_ERROR L2 Norm: " << error << std::endl;

				FSerror_r.exportToVTI("voxField_" + std::to_string(res) + "FS_ERROR_resized");
				FSerror.exportToVTI("voxField_" + std::to_string(res) + "FS_ERROR");
				AABBerror.exportToVTI("voxField_" + std::to_string(res) + "AABB_ERROR");

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

	SDF bunny_sdf_r = SDF(&bunny, res, true);

	std::cout << bunny_sdf_r.getComputationProperties();

	bunny_sdf_r.exportGrid(&e, "bunnySDF_r");

	return 1;
}

