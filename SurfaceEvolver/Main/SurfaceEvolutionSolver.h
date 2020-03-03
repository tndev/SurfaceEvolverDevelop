#ifndef SURFACEEVOLUTIONSOLVER_H_
#define SURFACEEVOLUTIONSOLVER_H_

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
	void getTriangleVertexLaplaceBeltrami(uint i, std::vector<std::vector<uint>>* adjacentPolys);
	void getQuadVertexLaplaceBeltrami(uint i, std::vector<std::vector<uint>>* adjacentPolys);
	// solver helpers
	void printArray1(std::string name, double* a, int printLim, bool inRow = true);
	void printArray2(std::string name, double** A, int printLim);
	double vectorDot(double* a, double* b);
	double vectorNorm(double* a);

	// solver
	void Bi_CGSTAB_solve(double** A, double* b, double* x, bool print = false);
public:
	// params:
	uint NSteps = 10; bool saveStates = false;
	size_t N;

	float dt = 0.01f; // time step
	Grid* sdfGrid = nullptr; // signed distance function grid
	ElementType type = ElementType::quad;

	// result (iterated):
	Geometry* evolvedSurface = nullptr;

	SurfaceEvolutionSolver();
	SurfaceEvolutionSolver(const SurfaceEvolutionSolver& other);
	SurfaceEvolutionSolver(Grid* sdfGrid, uint NSteps = 10, float dt = 0.01f, ElementType type = ElementType::quad, bool saveStates = false);
	~SurfaceEvolutionSolver();

	void init();
	void evolve();
};

#endif