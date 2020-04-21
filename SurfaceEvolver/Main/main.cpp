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
#include "../SDF/SDF.h"
#include "../Utils/CPUInfo.h"
#include "../EvolutionCore/Evolver.h"
#include "../EvolutionCore/Parameters.h"
#include "../EvolutionCore/EvolutionRemesher.h"

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
// - matrix multiplication for Matrix4
// - flood fill for sign computation of SDF
// - Test mesh angle weighted pseudonormals
// - test if grid gradient is computed correctly by exporting to a vtk vector file.
// - finite volume normal derivatives (Laplace-Beltrami)
// - compose a linear system for evolution from CubeSphere to PrimitiveBox of the same subdivision level
// - mean curvature flow for sphere test (cotan scheme)
// - fix lagging numerical solution for sphere test
// - test evolution without tangential redistribution on different objects
// - compute mean area for finite volumes in each time step
// - test a non-convex model (e.g.: bunny)
// - make separate outputs for mean co-volume measure
// - implement cutoff offset for the bounding cube to compute the field on minimum necessary subset (box)
// - visualize angle-weighted pseudo-normals with interpolated -grad(SDF) vectors
// - add scalar data (fvAreas, distances, curvatures) to mesh vertices
// - Special types: SEvolverParams, SDFParams,...
// - SurfaceEvolutionSolver -> Evolver, LinearSolver
// - refactor and separate console and log outputs for specific situations
// - catch all NaNs as exceptions (breaks)

//  POSTPONED:
//
// - implement a method/class to get CPU instruction set, mainly whether it supports AVX, an alternate resampling method has to be implemented for CPU's that do not support AVX
// - implement sort order function of a 256-bit AVX vector (needs a proper lookup hash)
// - AABB update for transformations + Timing test
// - Inverse transform grid upon transforming mesh
// - implement adaptive resampling for 512-bit registers - 2 * 8 sampling positions (if possible)
// - compare results with CGAL distance query implementation
// - sign computation \w (arbitrary)ray-mesh intersection (even # of intersections = 1, odd # of intersections = -1)
//   Important notes:
//		- split position matters, yet sometimes a less precise estimate yields better result than quad min. Trying 8 sample positions could help
//		- still no idea why the near/far cycle's written this way. Sometimes it leads the traversal to only one intersection, ignoring the rest, because they're "far"
// - flat AABB and Octree (we can try, but this would require dynamic arrays (i.e.: std::vectors) which would be slower than allocation)
// - fix apparent interpolation bug for SDF values and gradients (is there one?)

//   DONE, BUT MIGHT BE IMPROVED:
//
// - AABB and Octree have to take as little space as possible
// - implement generic triangle/quad centroid co-volume procedure (centroid on triangulated polygon)

//   WIP:
// 
// - test evolution for extremal cases: MCF dominant (eta = 0.01, eps = 1.0) and SDF dominant (eta = 1, eps = 0.01)


//   TODO:
//
// - co-volume measure-driven time step: dt ~ m(V)
//
// - fix SDF coordinates (use global grid indexing)
// - implement global grid and cellSize-based Octree & SDF (just like in Vctr Engine Meta Object)
// - finish class EvolutionRemesher with all params
//
// - quad co-volume scheme
// - mean curvature flow for sphere test (quad scheme)
// - mean curvature flow for sphere test (tri interp scheme)

void performSDFTest(uint res, Geometry& g, std::fstream& timing, VTKExporter& e) {
	std::cout << "init SDF..." << std::endl;

	// Fast sweeping DF, resized from 20 and interpolated
	SDF sdf_FS_r = SDF(&g, res, false, false, false, true, SDF_Method::fast_sweeping);

	std::cout << sdf_FS_r.getComputationProperties();
	timing << sdf_FS_r.getComputationProperties();

	sdf_FS_r.exportGrid(&e); // save to vti

	// Fast sweeping DF
	SDF sdf_FS = SDF(&g, res);

	std::cout << sdf_FS.getComputationProperties();
	timing << sdf_FS.getComputationProperties();

	sdf_FS.exportGrid(&e);

	// AABB DF
	SDF sdf_AABB = SDF(&g, res, false, false, false, false, SDF_Method::aabb_dist);

	std::cout << sdf_AABB.getComputationProperties();
	timing << sdf_AABB.getComputationProperties();

	sdf_AABB.exportGrid(&e);

	// Brute force DF
	SDF sdf_Brute = SDF(&g, res, false, false, false, false, SDF_Method::brute_force);

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


void performUnitSphereTest() {
	EvolutionParams eParams;
	eParams.tStop = 0.06f; eParams.elType = ElementType::tri;
	eParams.name = "testSphere"; eParams.saveStates = true;

	SphereTestParams stParams;

	std::fstream errLog("testSphere_errorLog.txt", std::fstream::out);
	errLog << "================================\n";
	errLog << ">>> Evolution error log ........\n";
	float errPrev, err, EOC;
	for (uint i = 0; i < 4; i++) {
		if (i > 0) errPrev = err;

		float dt = 0.01f / pow(4, i);
		eParams.dt = dt; eParams.NSteps = i + 1;
		stParams.testId = i;

		Evolver sphereTest(eParams, stParams);
		err = sphereTest.sphereTestL2Error;

		errLog << "dt = " << dt << ", Nsteps = " << sphereTest.NSteps << ", Nverts = " << sphereTest.N << std::endl;
		errLog << "L2Error = " << err;
		if (i > 0) {
			EOC = log2(errPrev / err);
			std::cout << "EOC = " << EOC << std::endl;
			errLog << ", EOC = " << EOC;
		}
		errLog << std::endl;
	}
	errLog.close();
}


int main()
{
	float r = 50.0f;
	unsigned int d = 3;
	float a = 2.0f * r / sqrt(3.0f);
	unsigned int ns = 10;

	VTKExporter e = VTKExporter();

	/*
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
				e.initExport(&g1, "cube" + std::to_string(res) + "-" + std::to_string(n));

				performSDFTest(res, g1, timing_cubes, e);

				std::cout << "cubesphere(" << n + 1 << "), grid_res = " << res << std::endl;
				CubeSphere g2 = CubeSphere(n + 1, r);
				g2.applyMatrix(Matrix4().makeRotationAxis(axis.x, axis.y, axis.z, M_PI / 6.));
				e.initExport(&g2, "cubesphere" + std::to_string(res) + "-" + std::to_string(n));

				performSDFTest(res, g2, timing_cubes, e);
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
				e.initExport(&g0, "icosphere" + std::to_string(res) + "-" + std::to_string(n));

				performSDFTest(res, g0, timing_ico, e);
			}
		}

		timing_ico.close();
	}*/

	// ===== BUNNY SDF tests =============================

	/*
	uint res = 40; // octree resolution

	// Vector3 axis = normalize(Vector3(1, 1, 1));
	
	auto startObjLoad = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	OBJImporter obj = OBJImporter();
	Geometry bunny = obj.importOBJGeometry("bunny_no_holes.obj");
	bunny.applyMatrix(Matrix4().setToScale(0.2f, 0.2f, 0.2f));
	e.initExport(&bunny, "sfBunny");
	// === Timed code ============
	auto endObjLoad = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedObj = (endObjLoad - startObjLoad);
	std::cout << "Model loaded after " << elapsedObj.count() << " seconds" << std::endl;

	
	SDF bunny_sdf = SDF(&bunny, res, true);

	std::cout << bunny_sdf.getComputationProperties();

	//bunny_sdf.exportGrid(&e, "bunnySDF");
	//bunny_sdf.exportGradientField(&e, "bunnySDF_grad");

	// ====== BUNNY Evolution =============================
	EvolutionParams eParams;
	eParams.name = "Bunny";
	eParams.dt = 0.03f; eParams.NSteps = 100; eParams.subdiv = (uint)4; eParams.elType = ElementType::tri;
	eParams.saveStates = true; eParams.printStepOutput = true; eParams.writeTimeLog = true;
	MeanCurvatureParams mcfParams;
	mcfParams.saveAreaStates = true; mcfParams.writeMeanAreaLog = true;
	GradDistanceParams sdfParams;
	sdfParams.targetGeom = &bunny; sdfParams.sdfGrid = bunny_sdf.grid;
	sdfParams.saveDistanceStates = true;
	sdfParams.saveGradientStates = true;

	Evolver evolver(eParams, mcfParams, sdfParams);*/

	// cube with holes
	/*
	OBJImporter obj = OBJImporter();
	Geometry cwh = obj.importOBJGeometry("cubeWithHoles.obj");
	cwh.applyMatrix(Matrix4().setToScale(0.02f, 0.02f, 0.02f));
	std::string name = "evolvingCubeWithHoles";
	e.initExport(&cwh, "cubeWithHoles");

	uint res = 40; // octree resolution
	SDF cwh_sdf = SDF(&cwh, res);
	// cwh_sdf.exportGrid(&e, "cubeWithHoles_SDF");
	// cwh_sdf.exportGradientField(&e, "cubeWithHoles_SDF_grad");

	std::cout << cwh_sdf.getComputationProperties();

	EvolutionParams eParams;
	eParams.name = name;
	eParams.dt = 0.03f; eParams.NSteps = 150; eParams.subdiv = (uint)4; eParams.elType = ElementType::tri;
	eParams.saveStates = true; eParams.printStepOutput = true; eParams.writeTimeLog = true;
	MeanCurvatureParams mcfParams;
	mcfParams.saveAreaStates = true; mcfParams.writeMeanAreaLog = true;
	GradDistanceParams sdfParams;
	sdfParams.targetGeom = &cwh; sdfParams.sdfGrid = cwh_sdf.grid;
	sdfParams.saveDistanceStates = true;
	sdfParams.saveGradientStates = true;
	*/

	/*
	// MCF dominant config (eta = 0.01, eps = 1.0)
	sdfParams.C = -0.01f;
	sdfParams.constant = true;

	mcfParams.C1 = 1.0f;
	mcfParams.constant = true;*/

	/*
	// SDF dominant config (eta = 1, eps = 0.01)
	sdfParams.C = -1.0f;
	sdfParams.constant = false;

	mcfParams.C1 = 0.01f;
	mcfParams.constant = true;*/

	// Evolver evolver(eParams, mcfParams, sdfParams);
	

	// arc
	/*
	OBJImporter obj = OBJImporter();
	Geometry arc = obj.importOBJGeometry("arc.obj");
	arc.applyMatrix(Matrix4().setToScale(0.02f, 0.02f, 0.02f));
	std::string name = "evolvingArc";
	e.initExport(&arc, "arc");

	uint res = 40; // octree resolution
	SDF arc_sdf = SDF(&arc, res);
	// cwh_sdf.exportGrid(&e, "cubeWithHoles_SDF");
	// cwh_sdf.exportGradientField(&e, "cubeWithHoles_SDF_grad");

	std::cout << arc_sdf.getComputationProperties();

	EvolutionParams eParams;
	eParams.name = name;
	eParams.dt = 0.03f; eParams.NSteps = 150; eParams.subdiv = (uint)4; eParams.elType = ElementType::tri;
	eParams.saveStates = true; eParams.printStepOutput = true; eParams.writeTimeLog = true;
	MeanCurvatureParams mcfParams;
	mcfParams.saveAreaStates = true; mcfParams.writeMeanAreaLog = true;
	GradDistanceParams sdfParams;
	sdfParams.targetGeom = &arc; sdfParams.sdfGrid = arc_sdf.grid;
	sdfParams.saveDistanceStates = true;
	sdfParams.saveGradientStates = true;

	Evolver evolver(eParams, mcfParams, sdfParams);
	*/
	/*
	e.exportGeometryVertexNormals(&bunny, "bunnyNormals");
	e.exportGeometryFiniteVolumeGrid(&bunny, "bunnyFVs");*/

	// performUnitSphereTest();

	/*
	std::string name = "testBox";
	PrimitiveBox b = PrimitiveBox(1, 1, 1, 3, 3, 3, true, name);
	Vector3 axis = Vector3(1, 1, 1);
	axis.normalize();
	Matrix4 M = Matrix4().makeRotationAxis(axis.x, axis.y, axis.z, M_PI / 6);
	b.applyMatrix(M);
	e.initExport(&b, name);
	uint res = 40; // octree resolution
	SDF boxSDF = SDF(&b, res, true, true);
	//boxSDF.exportGrid(&e, "boxSDF");
	//boxSDF.exportGradientField(&e, "boxSDF_Grad");

	EvolutionParams eParams;
	eParams.name = name;
	eParams.dt = 0.03f; eParams.NSteps = 150; eParams.subdiv = (uint)3; eParams.elType = ElementType::tri;
	eParams.saveStates = true; 
	eParams.printStepOutput = true; eParams.writeTimeLog = true;
	MeanCurvatureParams mcfParams;
	mcfParams.saveAreaStates = true; 
	mcfParams.writeMeanAreaLog = true;
	GradDistanceParams sdfParams;
	sdfParams.targetGeom = &b; sdfParams.sdfGrid = boxSDF.grid;
	sdfParams.saveDistanceStates = true;
	//sdfParams.saveGradientStates = true;
	
	Evolver evolver(eParams, mcfParams, sdfParams);*/

	/**/
	std::string name = "testEllipsoid";
	// CubeSphere cs = CubeSphere(10, 50.0f, true, name);
	IcoSphere is = IcoSphere(0, 1.0f, name);
	Matrix4 M = Matrix4().setToScale(1.5f, 1.0f, 1.0f);
	is.applyMatrix(M);
	e.initExport(&is, name);
	uint res = 40; // octree resolution
	SDF isSDF = SDF(&is, res, true, true);
	//isSDF.exportGrid(&e, name + "SDF");
	//isSDF.exportGradientField(&e, name + "SDF_Grad");

	EvolutionParams eParams;
	eParams.name = name;
	eParams.dt = 0.03f; eParams.NSteps = 150; eParams.subdiv = (uint)3; eParams.elType = ElementType::tri;
	eParams.saveStates = true;
	eParams.printStepOutput = true; eParams.writeTimeLog = true;
	MeanCurvatureParams mcfParams;
	mcfParams.saveAreaStates = true; 
	mcfParams.writeMeanAreaLog = true;
	GradDistanceParams sdfParams;
	sdfParams.targetGeom = &is; sdfParams.sdfGrid = isSDF.grid;
	sdfParams.saveDistanceStates = true;
	//sdfParams.saveGradientStates = true;

	Evolver evolver(eParams, mcfParams, sdfParams);
	
	
	/*
	std::string name = "testCubeSphere";
	CubeSphere cs = CubeSphere(4, 1.0f, true, name);
	Matrix4 M = Matrix4().setToScale(1.5f, 1.0f, 1.0f);
	cs.applyMatrix(M);
	e.initExport(&cs, name);
	uint res = 40; // octree resolution
	SDF csSDF = SDF(&cs, res, true, true);
	//isSDF.exportGrid(&e, name + "SDF");
	//isSDF.exportGradientField(&e, name + "SDF_Grad");

	EvolutionParams eParams;
	eParams.name = name;
	eParams.dt = 0.03f; eParams.NSteps = 150; eParams.subdiv = (uint)3; eParams.elType = ElementType::tri;
	eParams.saveStates = true;
	eParams.printStepOutput = true; eParams.writeTimeLog = true;
	MeanCurvatureParams mcfParams;
	mcfParams.saveAreaStates = true; 
	mcfParams.writeMeanAreaLog = true;
	GradDistanceParams sdfParams;
	sdfParams.targetGeom = &cs; sdfParams.sdfGrid = csSDF.grid;
	//sdfParams.saveDistanceStates = true;
	//sdfParams.saveGradientStates = true;

	Evolver evolver(eParams, mcfParams, sdfParams);*/

	/*
	IcoSphere is(1, 1.0f);
	e.initExport(&is, "icoSphere");
	std::vector<std::vector<Vector3>> fvVerts = {};
	std::vector<std::vector<std::vector<uint>>> adjacentPolys = {};
	is.getVertexFiniteVolumes(&fvVerts, &adjacentPolys);
	uint geomId = 0;
	for (int i = 0; i < is.uniqueVertices.size(); i++) {
		for (int j = 0; j < fvVerts[i].size(); j++) {
			e.exportGeometryFiniteVolumeGrid(&is, fvVerts, adjacentPolys, "icoSphereFV_" + std::to_string(geomId), i, j);
			if (i > 0) {
				e.exportGeometryFiniteVolumeGrid(&is, fvVerts, adjacentPolys, "icoSphereFVs_" + std::to_string(geomId), i - 1, -1, true);
			}
			else {
				Geometry empty = Geometry();
				e.initExport(&empty, "icoSphereFVs_" + std::to_string(geomId));
			}
			std::cout << "triangle " << geomId << " exported" << std::endl;
			geomId++;
		}
	}*/		

	/*
	uint i = 0;
	float dt = 0.01f / pow(4, i);
	Evolver sphereTest(dt, 0.06f, i, ElementType::tri, "testSphere", false, false, true);*/

	/*
	Matrix4 sdfTransform = Matrix4().makeTranslation(0.5, 0.5, 0.5).multiply(Matrix4().setToScale(2.0f, 2.0f, 2.0f));
	//.makeRotationAxis(axis.x, axis.y, axis.z, M_PI / 6.);
	bunny_sdf.applyMatrix(sdfTransform);

	std::cout << bunny_sdf.last_transform;

	bunny_sdf.exportGrid(&e, "bunnySDF_scaled");
	e.initExport(*bunny_sdf.geom, "sfBunny_scaled");*/

	/*
	PrimitiveBox cube = PrimitiveBox(a, a, a, 4, 4, 4);
	cube.applyMatrix(Matrix4().makeRotationAxis(axis.x, axis.y, axis.z, M_PI / 6.));
	e.initExport(cube, "cube");

	SDF cubeSDF = SDF(&cube, res, true); // signed dist to cube
	std::cout << std::endl << cubeSDF.getComputationProperties();
	cubeSDF.exportGrid(&e, "cubeSDF");*/
	

	// tree visualisation
	/*
	bunny_sdf.tri_aabb->GenerateFullTreeBoxVisualisation(e);
	bunny_sdf.tri_aabb->GenerateFullLeafBoxVisualisation(e);
	bunny_sdf.tri_aabb->GenerateStepwiseLeafBoxVisualisation(e);
	bunny_sdf.octree->GenerateFullOctreeBoxVisualisation(e);
	bunny_sdf.octree->GenerateLeafCellVisualisation(e); 
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

