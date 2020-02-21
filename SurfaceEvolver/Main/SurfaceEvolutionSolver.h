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
public:
	// params:
	uint NSteps = 10;
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
};

#endif