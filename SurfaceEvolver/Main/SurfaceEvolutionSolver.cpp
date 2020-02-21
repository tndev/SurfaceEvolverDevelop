#include "SurfaceEvolutionSolver.h"

SurfaceEvolutionSolver::SurfaceEvolutionSolver()
{
}

SurfaceEvolutionSolver::SurfaceEvolutionSolver(const SurfaceEvolutionSolver& other)
{
}

SurfaceEvolutionSolver::SurfaceEvolutionSolver(Grid* sdfGrid, uint NSteps, float dt, ElementType type, bool saveStates)
{
	// ===== Evolution params ==============================
	this->NSteps = NSteps;
	this->dt = dt;
	this->sdfGrid = sdfGrid;
	this->type = type;
	// -----------------------------------------------------

	this->init();


}

SurfaceEvolutionSolver::~SurfaceEvolutionSolver()
{
}

void SurfaceEvolutionSolver::init()
{
	// ===== Preparing a source geometry for immersion =====
	// -----------------------------------------------------
	// radius
	float r = 0.499f * std::min({ sdfGrid->scale.x, sdfGrid->scale.y, sdfGrid->scale.z });
	if (type == ElementType::tri) {
		uint n = 2;
		evolvedSurface = new IcoSphere(n, r);
	}
	else {
		uint n = 3;
		evolvedSurface = new CubeSphere(n, r);
	}
	Vector3 center = sdfGrid->bbox.getCenter();
	Matrix4 M = Matrix4().makeTranslation(center.x, center.y, center.z);
	evolvedSurface->applyMatrix(M);
	// -------------------------------------------------------
}
