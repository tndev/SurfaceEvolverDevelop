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
	std::shared_ptr<Grid> grid = nullptr;

	double f = 1.0; // rhs
	double h = 1.0; // dx, dy, dz

	uint Nsweeps = 0;

	FastSweep3D() = default;
	FastSweep3D(const FastSweep3D& other);
	FastSweep3D(std::shared_ptr<Grid>& grid, uint Nsweeps, bool saveGridStates = false, bool blur = false);
	~FastSweep3D() = default;

	void sweep(bool saveGridStates);
};

#endif
