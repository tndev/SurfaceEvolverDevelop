#ifndef FASTSWEEP3D_H_
#define FASTSWEEP3D_H_

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include "Grid.h"
#include "../ExportImport/VTKExporter.h"

// This is a solver for Eikonal equation |grad(u(x))| = 1 (i.e. the distance function)
// an Eikonal equation |grad(u(x))| = f(x) with a general rhs might require more than
// 2^n sweeps for dimension n since its characteristics are not straight lines.
class FastSweep3D
{
public:
	Grid* grid = nullptr;

	float f = 1.0f; // rhs
	float h = 1.0f; // dx, dy, dz

	uint Nsweeps = 0;

	FastSweep3D();
	FastSweep3D(const FastSweep3D& other);
	FastSweep3D(Grid* grid, uint Nsweeps, bool saveGridStates = false, bool blur = false);
	~FastSweep3D();

	void sweep(bool saveGridStates);
};

#endif
