#ifndef SURFACEEVOLUTIONSOLVER_H_
#define SURFACEEVOLUTIONSOLVER_H_

#include <chrono>
#include "../SDF/Grid.h"
#include "../Geometry/Geometry.h"
#include "../GeometryObject/Icosphere.h"
#include "../GeometryObject/CubeSphere.h"

enum class ElementType {
	tri = 0,
	quad = 1
};

// This is a solver for a Surface Evolution equation \partial_t F = v_N + v_T
// where v_N is the normal term and v_T is the tangential redistribution term
class SurfaceEvolutionSolver
{
private:
	double** SysMatrix = nullptr;
	double* sysRhsX = nullptr;
	double* sysRhsY = nullptr;
	double* sysRhsZ = nullptr;

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
public:
	// params:
	uint NSteps = 10; 
	size_t N;

	bool saveStates = false;
	bool printHappenings = true; // general things happening output
	bool printStepOutput = true; // time step output with time measurements etc.

	float dt = 0.01f; // time step
	Grid* sdfGrid = nullptr; // signed distance function grid
	ElementType type = ElementType::quad;

	// result (iterated):
	Geometry* evolvedSurface = nullptr;

	std::string geomName = "";
	std::vector<std::string> time_logs = {};

	SurfaceEvolutionSolver();
	// TODO: Finish copy
	SurfaceEvolutionSolver(const SurfaceEvolutionSolver& other);
	SurfaceEvolutionSolver(Grid* sdfGrid, uint NSteps = 10, float dt = 0.01f, ElementType type = ElementType::quad, std::string name = "Sphere", bool saveStates = false);
	// TODO: Finish delete
	~SurfaceEvolutionSolver();

	void init();
	void evolve();

	// returns a specialized L2error compared to a mean-curvature contracting sphere
	// with radius r(t) = sqrt(1 - 4 * t)
	float getSphereL2Error();
};

#endif