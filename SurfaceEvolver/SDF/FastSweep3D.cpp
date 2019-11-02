#include "FastSweep3D.h"

FastSweep3D::FastSweep3D()
{
}

FastSweep3D::FastSweep3D(Grid* grid)
{
	this->grid = grid;
}

FastSweep3D::~FastSweep3D()
{
}

void FastSweep3D::sweep(int dir[])
{
	uint Nx = this->grid->Nx; uint Ny = this->grid->Ny; uint Nz = this->grid->Nz;
	
	uint iXmin = dir[0] ? 0 : Nx;
	uint iXmax = dir[0] ? Nx : 0;

	uint iYmin = dir[1] ? 0 : Ny;
	uint iYmax = dir[1] ? Ny : 0;

	uint iZmin = dir[2] ? 0 : Nz;
	uint iZmax = dir[2] ? Nz : 0;

	uint gridPos, gridPosXprev, gridPosXnext, gridPosYprev, gridPosYnext, gridPosZprev, gridPosZnext;
	float u, uXprev, uXnext, uYprev, uYnext, uZprev, uZnext;
	float a, b, c;
	std::vector<float> uMins;

	for (uint iz = iZmin; iz < iZmax; iz += dir[2]) {
		for (uint iy = iYmin; iy < iYmax; iy += dir[1]) {
			for (uint ix = iXmin; ix < iXmax; ix += dir[0]) {
				
				gridPos = Nx * Ny * iz + Nx * iy + ix;

				gridPosXprev = Nx * Ny * iz + Nx * iy + std::max((int)ix - 1, 0);
				gridPosXnext = Nx * Ny * iz + Nx * iy + std::min((int)ix + 1, (int)Nx - 1);

				gridPosYprev = Nx * Ny * iz + Nx * std::max((int)iy - 1, 0) + ix;
				gridPosYnext = Nx * Ny * iz + Nx * std::min((int)iy + 1, (int)Ny - 1) + ix;

				gridPosZprev = Nx * Ny * std::max((int)iz - 1, 0) + Nx * iy + ix;
				gridPosZnext = Nx * Ny * std::min((int)iz + 1, (int)Nz - 1) + Nx * iy + ix;

				uXprev = this->grid->field[gridPosXprev];
				uXnext = this->grid->field[gridPosXnext];

				uYprev = this->grid->field[gridPosYprev];
				uYnext = this->grid->field[gridPosYnext];

				uZprev = this->grid->field[gridPosZprev];
				uZnext = this->grid->field[gridPosZnext];

				u = this->grid->field[gridPos];

				a = std::min(uXprev, uXnext);
				b = std::min(uYprev, uYnext);
				c = std::min(uZprev, uZnext);

				uMins = { a, b, c };
				std::sort_heap(uMins.begin(), uMins.end());
			}
		}
	}
}
