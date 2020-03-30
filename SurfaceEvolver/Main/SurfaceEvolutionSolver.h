#ifndef SURFACEEVOLUTIONSOLVER_H_
#define SURFACEEVOLUTIONSOLVER_H_

#include <chrono>
#include "../SDF/Grid.h"
#include "../Geometry/Geometry.h"
#include "../GeometryObject/Icosphere.h"
#include "../GeometryObject/CubeSphere.h"
#include "../ExportImport/VTKExporter.h"

enum class ElementType {
	tri = 0,
	quad = 1
};

// This is a solver for a Surface Evolution equation \partial_t F = v_N + v_T
// where v_N is the normal term and v_T is the tangential redistribution term
class SurfaceEvolutionSolver
{
private:
	// linear systems
	double** SysMatrix = nullptr;
	double* sysRhsX = nullptr;
	double* sysRhsY = nullptr;
	double* sysRhsZ = nullptr;

	// finite volume areas per vertex (vertices)
	std::vector<float> fvAreas = {};
	Vector3 center = Vector3(); // center of a test sphere geometry

	void clearSystem();
	void initSystem();

	void getInterpolatedSDFValuesforVertex(Vector3* V, float* SDF_V, Vector3* gradSDF_V, std::vector<Vector3>& positionBuffer, std::vector<float>& valueBuffer);
	// TODO: Add redistribution terms (low prio)
	float tangentialVelocitForVertex(Vector3& V, uint i);
	float laplaceBeltramiCtrlFunc(float& SDF_V);
	float etaCtrlFunc(float& SDF_V, Vector3& gradSDF_V, Vector3& nV);

	void getTriangleEvolutionSystem(
		std::vector<Vector3>& vNormals,	std::vector<std::vector<Vector3>>& fvVerts, std::vector<std::vector<std::vector<uint>>>& adjacentPolys
	);
	void getQuadEvolutionSystem(
		std::vector<Vector3>& vNormals, std::vector<std::vector<Vector3>>& fvVerts, std::vector<std::vector<std::vector<uint>>>& adjacentPolys
	);
	// solver helpers
	void printArray1(std::string name, double* a, int printLim, bool inRow = true);
	void printArray2(std::string name, double** A, int printLim);
	double vectorDot(double* a, double* b);
	double vectorNorm(double* a);

	// solver
	void Bi_CGSTAB_solve(double** A, double* b, double* x, bool print = false);

	// postproc
	void updateGeometry(double* Fx, double* Fy, double* Fz);
	void exportGeometry(int step);
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

	bool saveStates = false;
	bool printHappenings = true; // general things happening output
	bool printSolution = false; // whether to print Bi-CGStab solution output for each (xyz) component
	bool printStepOutput = true; // time step output with time measurements etc.

	// whether to compare evolution result with a mean-curvature contracting sphere and return an L2 error:
	bool sphereTest = false;
	float r0 = 1.0f; // test sphere initial radius

	// flag whether to consider a signed distance function for evolution equation (automatically true when performing sphere test)
	bool meanCurvatureFlow = false;

	float dt = 0.01f; // time step
	float tStop = 1.0f; // evolution stopping time
	Grid* sdfGrid = nullptr; // signed distance function grid
	ElementType type = ElementType::tri;

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

	SurfaceEvolutionSolver();
	// SurfaceEvolutionSolver(const SurfaceEvolutionSolver& other);

	// ----- test and applied variants of evolver constructor respectively ---------
	// sphere test evolution variant:
	SurfaceEvolutionSolver(
		float dt = 0.01f, float tStop = 1.0f, uint subdiv = 2, ElementType type = ElementType::tri, std::string name = "Sphere", int testId = -1,
		bool saveStates = false, bool printHappenings = false, bool printStepOutput = false, bool printSolution = false);

	// applied variant:
	SurfaceEvolutionSolver(
		float dt, int NSteps, uint subdiv, ElementType type = ElementType::tri,
		Geometry* targetGeom = nullptr, Grid* sdfGrid = nullptr, std::string name = "Sphere", float r0 = 1.0f,
		bool saveStates = false, bool printHappenings = false, bool printStepOutput = false, bool printSolution = false);

	~SurfaceEvolutionSolver();

	void init();
	void evolve();
};

#endif