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
// - test evolution for extremal cases: MCF dominant (eta = 0.01, eps = 1.0) and SDF dominant (eta = 1, eps = 0.01)
// - previous step mean curvature values H => H N = h (mean curvature vector)
// - previous step normal velocity term: v_N = eps(d) * h + eta(d) * N
// - tangential redist. lin. system: Laplace-Beltrami(psi) = dot( v_N, h ) - mean(dot( v_N , h) + omega * (A/G - 1)
// - Evolution from an input geometry
// - Sphere test with tangential redistribution TR (exceptionally slow for the last icosphere subdiv)
// - MCF-TR sphere test with an inhomogeneous distribution of vertices
// - MCF-TR(& w\o TR) evolution test on an ellipsoid
// - Angle-Based Tangential redistribution
// - polygon adjacency for boundary vertices

//  POSTPONED:
//
// - co-volume measure-driven time step: dt ~ m(V)
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
// - fix cotan scheme for boundary vertices


//   TODO:
//
// - fix finite volume computation for obtuse angles according to Desbrun
// - impose Dirichlet boundary conditions on Evolver
// - Debug volume-based TR - psi rhs has to have a zero sum
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

	// Fast sweeping D, resized from 20 and interpolated
	SDF sdf_FS_r = SDF(g, res, "", false, false, false, true, SDF_Method::fast_sweeping);

	std::cout << sdf_FS_r.getComputationProperties();
	timing << sdf_FS_r.getComputationProperties();

	sdf_FS_r.exportGrid(&e); // save to vti

	// Fast sweeping DF
	SDF sdf_FS = SDF(g, res, "");

	std::cout << sdf_FS.getComputationProperties();
	timing << sdf_FS.getComputationProperties();

	sdf_FS.exportGrid(&e);

	// AABB DF
	SDF sdf_AABB = SDF(g, res, "", false, false, false, false, SDF_Method::aabb_dist);

	std::cout << sdf_AABB.getComputationProperties();
	timing << sdf_AABB.getComputationProperties();

	sdf_AABB.exportGrid(&e);

	// Brute force DF
	SDF sdf_Brute = SDF(g, res, "", false, false, false, false, SDF_Method::brute_force);

	std::cout << sdf_Brute.getComputationProperties();
	timing << sdf_Brute.getComputationProperties();

	sdf_Brute.exportGrid(&e, "voxField_" + g.name + std::to_string(res));

	Grid FSerror_r = absGrid(subGrids(*sdf_FS_r.grid, *sdf_Brute.grid));
	Grid FSerror = absGrid(subGrids(*sdf_FS.grid, *sdf_Brute.grid));
	Grid AABBerror = absGrid(subGrids(*sdf_AABB.grid, *sdf_Brute.grid));

	double error = FSerror_r.getL2Norm();
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


void performSimpleSDFTest(Geometry& g, const uint octreeResolution, std::fstream& timing, VTKExporter& e) {
	std::cout << "init SDF for " << g.name << std::endl;

	// Fast sweeping DF
	SDF sdf_FS = SDF(g, octreeResolution, "");

	std::cout << sdf_FS.getComputationProperties();
	timing << sdf_FS.getComputationProperties();

	sdf_FS.exportGrid(&e);
}

SDFTimeLog performSDFTestWithOutput(Geometry& g, const uint octreeResolution, VTKExporter& e)
{
    std::cout << "init SDF for " << g.name << std::endl;

    const uint nAveragedCount = 10;
    SDFTimeLog averagedTimeLog{0, 0, 0, 0, 0, 0};

    std::cout << "Averaging time output for " << nAveragedCount << " runs:\n";

    uint gridExtent = 0;
	for (uint i = 0; i < nAveragedCount; i++)
	{
        std::cout << g.name << ", run " << (i + 1) << "...\n";
		auto sdf_FS = SDF(g, octreeResolution, e.pathPrefix, false);
        averagedTimeLog += sdf_FS.timeLog;
		std::cout << sdf_FS.getComputationProperties();

        if (i == 0) // save state of first simulation
        {
            gridExtent = sdf_FS.grid->gridExtent;
            sdf_FS.exportGrid(&e);
        }
	}

    averagedTimeLog /= nAveragedCount;
    averagedTimeLog.GridRes = gridExtent;
	
    return averagedTimeLog;
}

void performLagrangianEvolutionTest(
	Geometry& g, const uint octreeResolution, std::fstream& timing, VTKExporter& e, 
	EvolutionParams& evolParams, MeanCurvatureParams& mcfParams, 
	TangentialRedistParams& tRedistParams)
{
	if (evolParams.scale != 1.0)
	{
		g.applyMatrix(Matrix4().setToScale(evolParams.scale, evolParams.scale, evolParams.scale));		
	}

    e.initExport(&g, g.name);

	std::cout << "init SDF for " << g.name << std::endl;

	// Fast sweeping DF
	SDF sdf_FS = SDF(g, octreeResolution, e.pathPrefix, true);

	std::cout << sdf_FS.getComputationProperties();
	timing << sdf_FS.getComputationProperties();
	sdf_FS.exportGrid(&e);

	GradDistanceParams sdfParams;
	sdfParams.targetGeom = &g;
	sdfParams.sdfGrid = sdf_FS.grid.get();
	sdfParams.saveDistanceStates = true;

	Evolver evolver(evolParams, mcfParams, sdfParams, &tRedistParams);
}

void PerformFastSweepSdfTestForObjModel(const std::string& fileName, const uint& targetGridResolution)
{
	std::fstream timing(fileName + "_" + std::to_string(targetGridResolution) + ".txt", std::fstream::out);

	OBJImporter importer;
	Geometry geom = importer.importOBJGeometry(fileName);

	const uint octreeResolution = targetGridResolution;

	const bool computeSign = true;
	const bool computeGradient = false;
	const bool saveGridStates = false;
	const bool scaleAndInterpolate = false;

	SDF sdf_FS_r = SDF(geom, octreeResolution, "./CESCG_TestResults/",
		computeSign, computeGradient, saveGridStates, scaleAndInterpolate, 
		SDF_Method::fast_sweeping);

	std::cout << sdf_FS_r.getComputationProperties();
	timing << sdf_FS_r.getComputationProperties();

	VTKExporter e;
	std::cout << "Exporting to VTI: " << std::endl;
	sdf_FS_r.exportGrid(&e, fileName + "_" + std::to_string(targetGridResolution)); // save to vti
	std::cout << "... done" << std::endl;

	timing.close();
}

void performUnitSphereTest(bool tan_redistribute = false) {
	EvolutionParams eParams;
	eParams.tStop = 0.06; eParams.elType = ElementType::tri;
	eParams.name = "testSphere"; eParams.saveStates = false;

	SphereTestParams stParams;
	TangentialRedistParams* tanRedistParams = (tan_redistribute ? new TangentialRedistParams() : nullptr);
	tanRedistParams->type = 1;

	std::fstream errLog("testSphere_errorLog.txt", std::fstream::out);
	errLog << "================================\n";
	errLog << ">>> Evolution error log ........\n";
	double errPrev, err, EOC;
	for (uint i = 0; i < 4; i++) {
		if (i > 0) errPrev = err;

		double dt = 0.01 / pow(4, i);
		eParams.dt = dt;
		eParams.subdiv = i + 1;
		stParams.testId = i;

		Evolver sphereTest(eParams, stParams, tanRedistParams);
		err = sphereTest.testL2Error();

		errLog << "dt = " << dt << ", Nsteps = " << sphereTest.nSteps() << ", Nverts = " << sphereTest.nVerts() << std::endl;
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
    const std::vector<std::string> importedFilenames{
        //"armadillo.obj",
        //                      "blub.obj",
        "bunny.obj",
        "max-planck.obj",
        "nefertiti.obj",
        "ogre.obj",
        "spot.obj"
    };

    const std::string sourcePath = "./CESCG_TestData/";
    const std::string targetPath = "./CESCG_TestResults/";

    // ========= !!!!!!!!!!!!!!!!!!!!!!! ===================
	// >>>>>>> CHANGE THIS WHEN YOU RUN ON A DIFFERENT MACHINE
	// ========= !!!!!!!!!!!!!!!!!!!!!!! ===================
    const std::string cpuName = "IntelI7";

    VTKExporter vtkExp;
    vtkExp.pathPrefix = targetPath;
    OBJImporter objImp;
    objImp.pathPrefix = sourcePath;

    const bool performMeshSDFTests = true;
    const bool performEvolutionTests = false;
    const bool performUnitSphereTest = false;

    // ================== S D F    T E S T S ========================================
    if (performMeshSDFTests)
    {
        const std::vector<uint> octreeResolutions = {
            20, 30, 40, 50, 60, 70, 80, 90
        };

        /**/
        for (const auto& meshFileName : importedFilenames)
        {
            auto geom = objImp.importOBJGeometry(meshFileName);
            auto timing_file_aabb = std::fstream(targetPath + "timing_" + geom.name + "_aabb.txt", std::fstream::out);
            auto timing_file_octree = std::fstream(targetPath + "timing_" + geom.name + "_octree.txt", std::fstream::out);
            auto timing_file_fs = std::fstream(targetPath + "timing_" + geom.name + "_fs.txt", std::fstream::out);
            auto timing_file_flood = std::fstream(targetPath + "timing_" + geom.name + "_flood.txt", std::fstream::out);
            auto timing_file_total = std::fstream(targetPath + "timing_" + geom.name + "_total.txt", std::fstream::out);

            auto geomCapitalFirstLetter = std::string(1, toupper(geom.name[0]));
            const size_t geomNameSize = geom.name.size();
            timing_file_aabb << "timingAABB" + geomCapitalFirstLetter + geom.name.substr(1, geomNameSize) + cpuName + " = {";
            timing_file_octree << "timingOctree" + geomCapitalFirstLetter + geom.name.substr(1, geomNameSize) + cpuName + " = {";
            timing_file_fs << "timingFS" + geomCapitalFirstLetter + geom.name.substr(1, geomNameSize) + cpuName + " = {";
            timing_file_flood << "timingFlood" + geomCapitalFirstLetter + geom.name.substr(1, geomNameSize) + cpuName + " = {";
            timing_file_total << "timingTotal" + geomCapitalFirstLetter + geom.name.substr(1, geomNameSize) + cpuName + " = {";

            for (const auto& res : octreeResolutions)
            {
                const auto tLog = performSDFTestWithOutput(geom, res, vtkExp);
                timing_file_aabb << tLog.AABBTreeTime << (res != *(octreeResolutions.end() - 1) ? ", " : "");
                timing_file_octree << tLog.OctreeTime << (res != *(octreeResolutions.end() - 1) ? ", " : "");
                timing_file_fs << tLog.FastSweepingTime << (res != *(octreeResolutions.end() - 1) ? ", " : "");
                timing_file_flood << tLog.SignFloodFillTime << (res != *(octreeResolutions.end() - 1) ? ", " : "");
                timing_file_total << tLog.TotalTime << (res != *(octreeResolutions.end() - 1) ? ", " : "");
            }

            timing_file_aabb << "};";
            timing_file_octree << "};";
            timing_file_fs << "};";
            timing_file_flood << "};";
            timing_file_total << "};";
        	
            timing_file_aabb.close();
            timing_file_octree.close();
            timing_file_fs.close();
            timing_file_flood.close();
            timing_file_total.close();
        }
    }

    // ========= E V O L V E R    T E S T S ========================================

    if (performEvolutionTests)
    {
        const uint baselineOctreeRes = 40;

        // open Evolver timings:    
        auto timing_armadillo2 = std::fstream(targetPath + "timing_armadillo2.txt", std::fstream::out);
        auto timing_blub2 = std::fstream(targetPath + "timing_blub2.txt", std::fstream::out);
        auto timing_bunny2 = std::fstream(targetPath + "timing_bunny2.txt", std::fstream::out);
        auto timing_max_planck2 = std::fstream(targetPath + "timing_max-planck2.txt", std::fstream::out);
        auto timing_nefertiti2 = std::fstream(targetPath + "timing_nefertiti2.txt", std::fstream::out);
        auto timing_ogre2 = std::fstream(targetPath + "timing_ogre2.txt", std::fstream::out);
        auto timing_spot2 = std::fstream(targetPath + "timing_spot2.txt", std::fstream::out);

        std::map<uint, std::fstream*> indexedTimingFilesMap2{};
        indexedTimingFilesMap2[0] = &timing_armadillo2;
        indexedTimingFilesMap2[1] = &timing_blub2;
        indexedTimingFilesMap2[2] = &timing_bunny2;
        indexedTimingFilesMap2[3] = &timing_max_planck2;
        indexedTimingFilesMap2[4] = &timing_nefertiti2;
        indexedTimingFilesMap2[5] = &timing_ogre2;
        indexedTimingFilesMap2[6] = &timing_spot2;

        constexpr uint nIcoVerts0 = 12;
        constexpr uint nIcoTris0 = 20;
        constexpr uint nIcoEdges0 = 30;
        constexpr double octreeExpansionFactor = 1.1;
        constexpr double sdfGridExpandOffsetFactor = 1.0;
        constexpr double icoStartRadiusFactor = 0.4;

        const uint icoSubdiv = 3;

        const double dt = 0.01;
        const uint NSteps = 200;

        // =============================================
        // Evolver parameters:

        EvolutionParams evolParams;
        evolParams.subdiv = icoSubdiv;
        evolParams.dt = dt;
        evolParams.NSteps = NSteps;
        evolParams.name = "UnnamedEvolution"; // add geom name
        evolParams.saveStates = true;
        evolParams.printStepOutput = true;
        evolParams.writeTimeLog = true;
        evolParams.outputPath = targetPath;
        MeanCurvatureParams mcfParams;
        mcfParams.saveAreaStates = true;
        mcfParams.saveCurvatureStates = true;
        mcfParams.writeMeanAreaLog = true;
        TangentialRedistParams tanRedistParams;
        tanRedistParams.type = 0;
        tanRedistParams.omega_volume = 5.0;
        tanRedistParams.omega_angle = 0.5;

        // 0 - starting surface, 1 - end surface
        //const double stabilityEmphasis = -0.275;
        const double stabilityEmphasis = 0.0;

        uint timingId = 0; // a shitty way to iterate through multiple open timing files
        /**/
        for (const auto& meshFileName : importedFilenames)
        {
            auto geom = objImp.importOBJGeometry(meshFileName);

            // ========== scale factor estim ================

            auto bbox = geom.getBoundingBox();
            auto bboxSize = bbox.getSize();
            const double minBoxSize = std::min<double>({ bboxSize.x, bboxSize.y, bboxSize.z });
            const double limitIcoRadius = icoStartRadiusFactor * minBoxSize;
            const uint expectedVertexCount = (nIcoEdges0 * (pow(4, icoSubdiv) - 1) + 3 * nIcoVerts0) / 3;

            const double maxBoxSize = std::max<double>({ bboxSize.x, bboxSize.y, bboxSize.z });
            const double startingIcoRadius = icoStartRadiusFactor * octreeExpansionFactor * (minBoxSize + 4 * sdfGridExpandOffsetFactor * maxBoxSize);
            std::cout << "----------------------------------------------\n";
            std::cout << "minBoxSize = " << minBoxSize << "\n";
            std::cout << "startingIcoRadius = " << startingIcoRadius << "\n";
            std::cout << "limitIcoRadius = " << limitIcoRadius << "\n";
            std::cout << "----------------------------------------------\n";

            const double startingMeanCoVolArea = 4.0 * M_PI * startingIcoRadius * startingIcoRadius / expectedVertexCount;
            const double expectedMeanCoVolArea = 4.0 * M_PI * limitIcoRadius * limitIcoRadius / expectedVertexCount;
            std::cout << "expectedMeanCoVolArea = " << expectedMeanCoVolArea << "\n";
            std::cout << "startingMeanCoVolArea = " << startingMeanCoVolArea << "\n";
            const double factoredMeanCoVolArea = stabilityEmphasis * expectedMeanCoVolArea + (1.0 - stabilityEmphasis) * startingMeanCoVolArea;
            std::cout << "factoredMeanCoVolArea = " << factoredMeanCoVolArea << "\n";

            const double scaleFactor = pow(dt / factoredMeanCoVolArea, 1.0 / 3.0);
            std::cout << "scaleFactor = " << scaleFactor << "\n";
            std::cout << "----------------------------------------------\n";

            const auto bboxCenter = bbox.getCenter();
            std::cout << "bboxCenter = " << bboxCenter << "\n";
            Matrix4 geomTransformationMatrix(
                scaleFactor, 0.0, 0.0, -bboxCenter.x,
                0.0, scaleFactor, 0.0, -bboxCenter.y,
                0.0, 0.0, scaleFactor, -bboxCenter.z,
                0.0, 0.0, 0.0, 1.0
            );

            geom.applyMatrix(geomTransformationMatrix);

            evolParams.name = geom.name;

            performLagrangianEvolutionTest(geom, baselineOctreeRes, *indexedTimingFilesMap2[timingId],
                vtkExp, evolParams, mcfParams, tanRedistParams);
            timingId++;
        }

    }

    /*
    uint res = 40; // octree resolution

    auto startObjLoad = std::chrono::high_resolution_clock::now();
    // === Timed code ============
    OBJImporter obj = OBJImporter();
    obj.pathPrefix = sourcePath;
    Geometry geom = obj.importOBJGeometry("nefertiti.obj");
    const double scaleF = 0.03 / 4.0;
    std::cout << "scaleFactor = " << scaleF << "\n";
    geom.applyMatrix(Matrix4().setToScale(scaleF, scaleF, scaleF));
    VTKExporter e;
    e.pathPrefix = targetPath;
    e.initExport(&geom, geom.name);
    // === Timed code ============
    auto endObjLoad = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedObj = (endObjLoad - startObjLoad);
    std::cout << "Model loaded after " << elapsedObj.count() << " seconds" << std::endl;


    SDF geom_sdf = SDF(geom, res, targetPath, true);

    std::cout << geom_sdf.getComputationProperties();

    //geom_sdf.exportGrid(&e, "geomSDF");
    // geom_sdf.exportGradientField(&e, "geomSDF_grad");

    // ====== geom Evolution =============================
    EvolutionParams eParams;
    eParams.name = geom.name;
    eParams.dt = 0.01; eParams.NSteps = 200; eParams.subdiv = (uint)3; eParams.elType = ElementType::tri;
    eParams.saveStates = true; eParams.printStepOutput = true; eParams.writeTimeLog = true;
    eParams.outputPath = targetPath;
    MeanCurvatureParams mcfParams;
    mcfParams.saveAreaStates = true;
    mcfParams.saveCurvatureStates = true;
    mcfParams.writeMeanAreaLog = true;
    GradDistanceParams sdfParams;
    sdfParams.targetGeom = &geom; sdfParams.sdfGrid = geom_sdf.grid.get();
    sdfParams.saveDistanceStates = true;
    // sdfParams.saveGradientStates = true;
    //mcfParams.smoothSteps = 10;s
    TangentialRedistParams tRedistParams;
    tRedistParams.type = 0;
    tRedistParams.omega_volume = 0.5;
    tRedistParams.omega_angle = 0.5;
    //tRedistParams.saveTangentialVelocityStates = true;

    Evolver evolver(eParams, mcfParams, sdfParams, &tRedistParams);*/

    return 1;
}


