#ifndef EVOLVER_H_
#define EVOLVER_H_

#include <chrono>
#include "Parameters.h"
#include "../Utils/LinearSolver.h"
#include "../SDF/Grid.h"
#include "../Geometry/Geometry.h"
#include "../GeometryObject/Icosphere.h"
#include "../GeometryObject/CubeSphere.h"
#include "../ExportImport/VTKExporter.h"

// This is a solver for a Surface Evolution equation \partial_t F = v_N + v_T
// where v_N is the normal term and v_T is the tangential redistribution term
class Evolver
{
private:
	// linear systems
	double** SysMatrix = nullptr;
	double* sysRhsX = nullptr;
	double* sysRhsY = nullptr;
	double* sysRhsZ = nullptr;

	LinearSolver* solveX = nullptr;
	LinearSolver* solveY = nullptr;
	LinearSolver* solveZ = nullptr;

	Vector3 center = Vector3(); // center of a test sphere geometry
	
	// mesh data
	std::vector<float> fvAreas = {}; // vertex finite volume areas
	std::vector<float> vCurvatures = {}; // vertex curvature scalars
	std::vector<Vector3> vNormals = {}; // vertex normals
	std::vector<float> vDistances = {}; // vertex distances to target mesh
	std::vector<Vector3> vGradients = {}; // vertex gradients of distance func to target mesh
	std::vector<float> vDotProducts = {}; // dot(-grad(SDF), N) for each vertex

	// ctrl params:
	float rDecay = 1.0f; // radius for the mean curvature flow exponential decay parameter (C2 in eta(SDF))
	float C1 = 1.0f;
	float C2 = rDecay;
	float C = -1.0f; // C < 0 (-grad(SDF)); C > 0 (+grad(SDF))

	bool epsConstant = false; // Laplace-Beltrami func admits a constant value C1;
	bool etaConstant = false; // eta ctrl func (SDF) admits a constant value C;

	// ====== evolution methods ==========
	void init();
	void evolve();

	void clearSystem();
	void initSystem();

	void getInterpolatedSDFValuesforVertex(Vector3* V, float* SDF_V, Vector3* gradSDF_V, std::vector<Vector3>& positionBuffer, std::vector<float>& valueBuffer);
	float tangentialVelocitForVertex(Vector3& V, uint i);

	void saveFVAreaScalars();
	void saveInterpolatedSDFValues();
	void saveInterpolatedDotValues();
	void saveInterpolatedSDFGradients();

	float laplaceBeltramiCtrlFunc(float& SDF_V);
	float etaCtrlFunc(float& SDF_V, Vector3& gradSDF_V, Vector3& nV);

	void getTriangleEvolutionSystem(
		std::vector<std::vector<Vector3>>& fvVerts, std::vector<std::vector<std::vector<uint>>>& adjacentPolys, float& meanArea
	);
	void getQuadEvolutionSystem(
		std::vector<std::vector<Vector3>>& fvVerts, std::vector<std::vector<std::vector<uint>>>& adjacentPolys, float& meanArea
	);

	// postproc
	void updateGeometry(double* Fx, double* Fy, double* Fz);
	void exportGeometry(int step);
	void exportVectorStates(int step);
	void exportTestGeometry(int step, float t); // exports an ico/quad sphere determined by r(t) = sqrt(r0 * r0 - 4 * t) for comparison

	// returns a specialized error compared to a mean-curvature contracting sphere
	// with radius r(t) = sqrt(r0 * r0 - 4 * t)
	float getSphereStepL2Error(float t);
	float getSphereStepError(float t);
public:
	// params:
	uint subdiv = 2; // initial sphere subdivision detail
	uint NSteps = 10; 
	size_t N = 0;

	// output flags:
	bool saveStates = false; 
	bool saveAreaStates = false;  
	bool saveDistanceStates = false; 
	bool saveGradientStates = false; 
	bool saveCurvatureStates = false;

	bool printHappenings = true; // general things happening output
	bool printSolution = false; // whether to print Bi-CGStab solution output for each (xyz) component
	bool printStepOutput = true; // time step output with time measurements etc.

	// specific logs
	bool writeGenericLog = true;
	bool writeErrorLog = false;
	bool writeMeanAreaLog = false;
	bool writeTimeLog = false;

	// whether to compare evolution result with a mean-curvature contracting sphere and return an L2 error:
	bool sphereTest = false;
	float r0 = 1.0f; // test sphere initial radius

	// flag whether to consider a signed distance function for evolution equation (automatically true when performing sphere test)
	bool meanCurvatureFlow = false;

	float dt = 0.01f; // time step
	float tStop = 1.0f; // evolution stopping time
	Grid* sdfGrid = nullptr; // signed distance function grid
	ElementType elType = ElementType::tri;

	Geometry* targetGeom = nullptr;
	// result (iterated):
	Geometry* evolvedSurface = nullptr;

	// L2 error for numerical tests:
	float sphereTestL2Error = 0.0f;

	// id of the sphere test
	int testId = -1;

	std::string geomName = "";
	std::string log_header = "";
	std::string time_log = "";

	Evolver();
	// Evolver(const Evolver& other);
	// ----- test and applied variants of evolver constructor respectively ---------
	// sphere test evolution variant:
	Evolver(EvolutionParams& eParams, SphereTestParams& stParams);
	// applied variant:
	Evolver(EvolutionParams& eParams, MeanCurvatureParams& mcfParams, GradDistanceParams& sdfParams);
	~Evolver();
};

#endif