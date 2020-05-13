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

	double* psi = nullptr;

	LinearSolver* solveX = nullptr;
	LinearSolver* solveY = nullptr;
	LinearSolver* solveZ = nullptr;

	Vector3 center = Vector3(); // center of a test sphere geometry
	
	// mesh data
	std::vector<double> fvAreas = {}; // vertex finite volume areas
	std::vector<double> vCurvatures = {}; // vertex curvature scalars
	std::vector<Vector3> vCurvatureVectors = {}; // vertex curvature vectors
	std::vector<Vector3> vNormals = {}; // vertex normals
	std::vector<double> vDistances = {}; // vertex distances to target mesh
	std::vector<Vector3> vGradients = {}; // vertex gradients of distance func to target mesh
	std::vector<double> vDotProducts = {}; // dot(-grad(SDF), N) for each vertex
	std::vector<Vector3> vNormalVelocityVectors = {}; // normal velocity vectors: v_N = eps(d) * h + eta(d) * N
	std::vector<double> vNormalVelocities = {}; // normal velocity scalars: ||v_N|| = ||eps(d) * h + eta(d) * N||
	std::vector<Vector3> vTangentialVelocities = {}; // tangential velocity data
	//std::vector<Vector3> fvNormals = {}; // finite volume normals

	//Geometry* dummyFVGeom = new Geometry();

	double totalArea;
	double totalCurvatureSqWeightedArea;

	// essential mesh data
	std::vector<std::vector<Vector3>> fvVerts = {};
	std::vector<std::vector<std::vector<uint>>> adjacentPolys = {};

	// log files:
	std::fstream log;
	std::fstream meanAreaLog;
	std::fstream errorLog; // for sphereTest
	std::fstream timingLog;

	std::fstream psi_log; // for redistribution potential
	std::fstream psi_Sys_log; // for lin. system data of the redist. potential
	std::fstream sum_rhs_log; // for verifying that redistribution potential system rhs has zero sum

	// timer
	std::chrono::high_resolution_clock::time_point startGlobalTime;
	std::chrono::duration<double> elapsedTimeTotal;

	// ctrl params:
	double rDecay = 1.0; // radius for the mean curvature flow exponential decay parameter (C2 in eta(SDF))
	double C1 = 1.0;
	double C2 = rDecay;
	double C = -1.0; // multiplies eta(SDF) function
	double D = 0.0; // multiplies sqrt(1 - dot(grad(SDF), N)^2)
	double initSmoothRate = 0.5;
	double smoothDecay = 1.0;

	// tangential redistribution:
	double omega_volume = 100.0;
	double omega_angle = 2.0;
	int redistribution_type = -1;

	bool epsConstant = false; // Laplace-Beltrami func admits a constant value C1;
	bool etaConstant = false; // eta ctrl func (SDF) admits a constant value C;

	bool start_flag = true;
	bool smoothing_flag = false;
	bool interruptEvolution = false;

	int tBegin = 1;
	int degenerateCount = 0;

	// ====== evolution methods ==========
	void init();
	void evolve();

	void clearSystem();
	void initSystem();

	void getInterpolatedSDFValuesforVertex(Vector3* V, double* SDF_V, Vector3* gradSDF_V, std::vector<Vector3>& positionBuffer, std::vector<double>& valueBuffer);
	Vector3 getVolumeTangentialVelocityForVertex(Vector3& V, uint i);
	Vector3 getAngleTangentialVelocityForVertex(Vector3& V, uint i);

	void saveFVAreaScalars();
	void saveInterpolatedSDFValues();
	void saveInterpolatedDotValues();
	void saveInterpolatedSDFGradients();
	void saveMeanCurvatureScalars();
	void saveNormalVelocityScalars();
	void saveMeanCurvatureVectors(int step);
	void saveTangentialVelocityVectors(int step);
	void saveNormalVelocityVectors(int step);
	//void saveFVNormals(int step);
	void saveRedistributionPotential();

	double laplaceBeltramiCtrlFunc(double& SDF_V);
	double laplaceBeltramiSmoothFunc(double t);
	double etaCtrlFunc(double& SDF_V, Vector3& gradSDF_V, Vector3& nV);
	double tangentialRedistDecayFunction(double& SDF_V);
	double tangentialRedistCurvatureFunction(double& H);

	void computeSurfaceNormalsAndCoVolumes();
	void getTriangleEvolutionSystem(double smoothStep, double& meanArea);
	void getQuadEvolutionSystem(double smoothStep, double& meanArea);
	void getTriangleRedistributionSystem();
	void getQuadRedistributionSystem();
	void getCurvaturesAndNormalVelocities();

	// logging
	void openLogs();
	void closeLogs();

	// postproc
	void updateGeometry(double* Fx, double* Fy, double* Fz);
	void exportGeometry(int step);
	void exportResultGeometry(std::string suffix = "_result");
	void exportVectorStates(int step);
	void exportTestGeometry(int step, double t); // exports an ico/quad sphere determined by r(t) = sqrt(r0 * r0 - 4 * t) for comparison

	// returns a specialized error compared to a mean-curvature contracting sphere
	// with radius r(t) = sqrt(r0 * r0 - 4 * t)
	double getSphereStepL2Error(double t);
	double getSphereStepError(double t);

	// params:
	uint subdiv = 2; // initial sphere subdivision detail
	uint NSteps = 10;
	size_t N = 0;

	int smoothSteps = -1; // MCF smooth steps

	// output flags:
	bool saveStates = false;
	bool saveResult = true;

	bool saveAreaStates = false;
	bool saveDistanceStates = false;
	bool saveGradientStates = false;
	bool saveCurvatureStates = false;
	bool saveNormalVelocityStates = false;
	bool saveCurvatureVectors = false;
	bool saveTangentialVelocityStates = false;

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
	double r0 = 1.0; // test sphere initial radius

	// flag whether to consider a signed distance function for evolution equation (automatically true when performing sphere test)
	bool meanCurvatureFlow = false;

	double dt = 0.01; // time step
	double tStop = 1.0; // evolution stopping time
	Grid* sdfGrid = nullptr; // signed distance function grid
	ElementType elType = ElementType::tri;

	Geometry* targetGeom = nullptr;
	// result (iterated):
	Geometry* evolvedSurface = nullptr;
	Box3 bbox;

	// L2 error for numerical tests:
	double sphereTestL2Error = 0.0;

	// id of the sphere test
	int testId = -1;

	std::string geomName = "";
	std::string log_header = "";
	std::string time_log = "";
public:
	Evolver();
	// Evolver(const Evolver& other);
	// ----- test and applied variants of evolver constructor respectively ---------
	// sphere test evolution variant:
	Evolver(EvolutionParams& eParams, SphereTestParams& stParams, TangentialRedistParams* tanParams = nullptr);
	// applied variant:
	Evolver(EvolutionParams& eParams, MeanCurvatureParams& mcfParams, GradDistanceParams& sdfParams, TangentialRedistParams* tanParams = nullptr);
	~Evolver();

	// getters:
	inline double testL2Error() { return sphereTestL2Error; };
	inline uint nSteps() { return NSteps; };
	inline uint nVerts() { return N; };
};

#endif