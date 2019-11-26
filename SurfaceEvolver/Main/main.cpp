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
// - interpolate distance field for higher resolutions (not much improvement)
// - perform simple DF tests for geom primitives like sphere, icosphere, cubesphere
// - optimize Box-Triangle intersection
// - add Quaternion class and TRS decomposition of Matrix4
// - adaptive resampling of split cost function (done for 2 * 4 samples - 265-bit registers)
// - minimize split cost function using a piecewise-quadratic interpolation to find the minimum (20% slower than simple cost(x) < minCost comparison)

//  POSTPONED:
//
// - implement sort order function of a 256-bit AVX vector
// - AABB update for transformations + Timing test
// - Inverse transform grid upon transforming mesh
// - implement adaptive resampling for 512-bit registers - 2 * 8 sampling positions (if possible)
// - compare results with CGAL distance query implementation


//   DONE, BUT MIGHT BE IMPROVED:
//
// - AABB and Octree have to take as little space as possible

//   WIP:
// 
// - flat AABB and Octree


//   TODO:
//
// - sign computation \w (arbitrary)ray-mesh intersection (even # of intersections = 1, odd # of intersections = -1)
// - implement a method/class to get CPU instruction set, mainly whether it supports AVX, an alternate resampling method has to be implemented for CPU's that do not support AVX
// - implement cutoff offset for the bounding cube to compute the field on minimum necessary subset (box)
// - compose a linear system for evolution from CubeSphere to PrimitiveBox of the same subdivision level
// - implement a VTK window form using a working example for mesh rendering and SDF volume rendering

void performTest(uint res, Geometry& g, std::fstream& timing, VTKExporter& e) {
	std::cout << "init SDF..." << std::endl;

	// Fast sweeping DF, resized from 20 and interpolated
	SDF sdf_FS_r = SDF(&g, res, false, true, SDF_Method::fast_sweeping);

	std::cout << sdf_FS_r.getComputationProperties();
	timing << sdf_FS_r.getComputationProperties();

	sdf_FS_r.exportGrid(&e); // save to vti

	// Fast sweeping DF
	SDF sdf_FS = SDF(&g, res);

	std::cout << sdf_FS.getComputationProperties();
	timing << sdf_FS.getComputationProperties();

	sdf_FS.exportGrid(&e);

	// AABB DF
	SDF sdf_AABB = SDF(&g, res, false, false, SDF_Method::aabb_dist);

	std::cout << sdf_AABB.getComputationProperties();
	timing << sdf_AABB.getComputationProperties();

	sdf_AABB.exportGrid(&e);

	// Brute force DF
	SDF sdf_Brute = SDF(&g, res, false, false, SDF_Method::brute_force);

	std::cout << sdf_Brute.getComputationProperties();
	timing << sdf_Brute.getComputationProperties();

	sdf_Brute.exportGrid(&e, "voxField_" + g.name + std::to_string(res));

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
	std::cout << "AABB_ERROR L2 Norm: " << error << std::endl << std::endl;
	timing << "AABB_ERROR L2 Norm: " << error << std::endl;

	FSerror_r.exportToVTI("voxField_" + g.name + std::to_string(res) + "FS_ERROR_resized");
	FSerror.exportToVTI("voxField_" + g.name + std::to_string(res) + "FS_ERROR");
	AABBerror.exportToVTI("voxField_" + g.name + std::to_string(res) + "AABB_ERROR");
}


int main()
{
	float r = 50.0f;
	unsigned int d = 3;
	float a = 2.0f * r / sqrt(3.0f);
	unsigned int ns = 10;

	VTKExporter e = VTKExporter();

	bool iterateCubeSphereTest = false;

	if (iterateCubeSphereTest) {
		size_t min_Res = 20, max_Res = 60;
		size_t min_Ns = 0, max_Ns = 5;
		std::fstream timing_cubes("timing_cubes.txt", std::fstream::out);

		Vector3 axis = normalize(Vector3(1, 1, 1));
        
		for (uint n = min_Ns; n < max_Ns; n++) {
			for (uint i = 0; i <= 2; i++) {
				uint res = min_Res * pow(2, i);

				std::cout << "cube(" << n + 1 << "), grid_res = " << res << std::endl;
				PrimitiveBox g1 = PrimitiveBox(a, a, a, n + 1, n + 1, n + 1);
				g1.applyMatrix(Matrix4().makeRotationAxis(axis.x, axis.y, axis.z, M_PI / 6.));
				e.initExport(g1, "cube" + std::to_string(res) + "-" + std::to_string(n));

				performTest(res, g1, timing_cubes, e);

				std::cout << "cubesphere(" << n + 1 << "), grid_res = " << res << std::endl;
				CubeSphere g2 = CubeSphere(n + 1, r);
				g2.applyMatrix(Matrix4().makeRotationAxis(axis.x, axis.y, axis.z, M_PI / 6.));
				e.initExport(g2, "cubesphere" + std::to_string(res) + "-" + std::to_string(n));

				performTest(res, g2, timing_cubes, e);
			}
		}
		timing_cubes.close();
	}


	bool iterateIcoSphereTest = false;

	if (iterateIcoSphereTest) {
		size_t min_Res = 20, max_Res = 60;
		size_t min_Ns = 0, max_Ns = 2;
		std::fstream timing_ico("timing_ico.txt", std::fstream::out);

		Vector3 axis = normalize(Vector3(1, 1, 1));

		for (uint n = min_Ns; n < max_Ns; n++) {
			for (uint i = 0; i <= 2; i++) {
				uint res = min_Res * pow(2, i);

				std::cout << "icosphere(" << n << "), grid_res = " << res << std::endl;
				IcoSphere g0 = IcoSphere(n, r);
				g0.applyMatrix(Matrix4().makeRotationAxis(axis.x, axis.y, axis.z, M_PI / 6.));
				e.initExport(g0, "icosphere" + std::to_string(res) + "-" + std::to_string(n));

				performTest(res, g0, timing_ico, e);
			}
		}

		timing_ico.close();
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

	Vector3 axis = normalize(Vector3(1, 1, 1));
	Matrix4 sdfTransform = Matrix4().makeTranslation(0.5, 0.5, 0.5);
	sdfTransform.setToScale(2.0f, 2.0f, 2.0f);
	//.makeRotationAxis(axis.x, axis.y, axis.z, M_PI / 6.);
	bunny_sdf.applyMatrix(sdfTransform);

	std::cout << bunny_sdf.last_transform;

	bunny_sdf.exportGrid(&e, "bunnySDF_scaled");
	e.initExport(*bunny_sdf.geom, "sfBunny_scaled");

	// tree visualisation
	/* bunny_sdf.tri_aabb->GenerateFullTreeBoxVisualisation(e);
	bunny_sdf.tri_aabb->GenerateFullLeafBoxVisualisation(e);
	bunny_sdf.tri_aabb->GenerateStepwiseLeafBoxVisualisation(e);
	bunny_sdf.octree->GenerateFullOctreeBoxVisualisation(e);
	bunny_sdf.octree->GenerateLeafCellVisualisation(e); */

	/* Interpolated bunny DF:
	SDF bunny_sdf_r = SDF(&bunny, res, false, true);

	std::cout << bunny_sdf_r.getComputationProperties();

	bunny_sdf_r.exportGrid(&e, "bunnySDF_r");
	*/

	/* The brute force DF of the bunny model will take ~27 min ! 
	
	SDF bunny_sdf_b = SDF(&bunny, res, false, SDF_Method::brute_force);

	std::cout << bunny_sdf_b.getComputationProperties();

	bunny_sdf_b.exportGrid(&e, "bunnySDF_b");

	Grid bunnySDF_r_Error = absGrid(subGrids(*bunny_sdf_r.grid, *bunny_sdf_b.grid));
	bunnySDF_r_Error.exportToVTI("bunnySDF_" + std::to_string(res) + "FS_ERROR_resized");

	Grid bunnySDF_Error = absGrid(subGrids(*bunny_sdf.grid, *bunny_sdf_b.grid));
	bunnySDF_Error.exportToVTI("bunnySDF_" + std::to_string(res) + "FS_ERROR");

	float error = bunnySDF_r_Error.getL2Norm();
	std::cout << "FS_ERROR_resized L2 Norm: " << error << std::endl;

	error = bunnySDF_Error.getL2Norm();
	std::cout << "FS_ERROR L2 Norm: " << error << std::endl; */

	return 1;
}

