#ifndef FASTSWEEP3D_H_
#define FASTSWEEP3D_H_

#include <algorithm>
#include <chrono>
#include "Grid.h"

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

	// 8 sweeping directions
	int sweepDir[8][3] = {
		{-1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1},
		{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1}
	};

	FastSweep3D();
	FastSweep3D(Grid* grid, uint Nsweeps);
	~FastSweep3D();

	float EikonalSolveInDim(std::vector<float>& aValues, uint dim);
	void sweep(int dir[]);
};

#endif
