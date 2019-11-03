#include "FastSweep3D.h"

FastSweep3D::FastSweep3D()
{
}

FastSweep3D::FastSweep3D(Grid* grid, uint Nsweeps, bool blur)
{
	this->grid = grid;
	this->h = this->grid->scale.x / this->grid->Nx;
	this->Nsweeps = Nsweeps;

	auto startFastSweep = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	std::cout << "Initiating FastSweep3D..." << std::endl;
	for (uint i = 0; i < this->Nsweeps; i++) {
		std::cout << "sweep " << i << " ... " << std::endl;
		this->sweep(this->sweepDir[i]);
	}
	if (blur) this->grid->blur();
	// === Timed code ============
	auto endFastSweep = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedFastSweep = (endFastSweep - startFastSweep);
	std::cout << "FastSweep3D finished after " << elapsedFastSweep.count() << " seconds" << std::endl;
}

FastSweep3D::~FastSweep3D()
{
}

float FastSweep3D::EikonalSolveInDim(std::vector<float>& aValues, uint dim)
{
	if (dim == 1) return aValues[0] + h * f;

	float sumA = 0.0f, sumASq = 0.0f;
	for (uint i = 0; i < dim; i++) {
		sumA += aValues[i];
		sumASq += aValues[i] * aValues[i];
	}

	float a = dim;
	float b = -2.0f * sumA;
	float c = sumASq - h * h * f * f;
	float discriminant = b * b - 4.0f * a * c;
	if (discriminant < 0.0f) {
		return LARGE_VAL;
	}
	else {
		return (-b + sqrt(discriminant)) / (2.0f * a);
	}
}

void FastSweep3D::sweep(int dir[])
{
	uint Nx = this->grid->Nx; uint Ny = this->grid->Ny; uint Nz = this->grid->Nz;
	
	int iXmin = dir[0] > 0 ? 0 : Nx - 1;
	int iXmax = dir[0] > 0 ? Nx : -1;

	int iYmin = dir[1] > 0 ? 0 : Ny - 1;
	int iYmax = dir[1] > 0 ? Ny : -1;

	int iZmin = dir[2] > 0 ? 0 : Nz - 1;
	int iZmax = dir[2] > 0 ? Nz : -1;

	std::cout << "sweep ranges: iX = [" << iXmin << ", " << iXmax << ", step = " << dir[0] << "]" << std::endl;
	std::cout << "sweep ranges: iY = [" << iYmin << ", " << iYmax << ", step = " << dir[1] << "]" << std::endl;
	std::cout << "sweep ranges: iZ = [" << iZmin << ", " << iZmax << ", step = " << dir[2] << "]" << std::endl;

	uint adim, gridPos, gridPosXprev, gridPosXnext, gridPosYprev, gridPosYnext, gridPosZprev, gridPosZnext;
	float u, uNew, uXprev, uXnext, uYprev, uYnext, uZprev, uZnext;
	float uMin, uMinX, uMinY, uMinZ;
	std::vector<float> neighborVals;
	std::vector<float> uMins;

	auto SolveEikonal = [&](uint ix, uint iy, uint iz) -> float {
		adim = 3;

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

		neighborVals.clear();

		neighborVals.push_back(std::min(uXprev, uXnext));
		neighborVals.push_back(std::min(uYprev, uYnext));
		neighborVals.push_back(std::min(uZprev, uZnext));

		uMins.clear();
		for (uint d = 0; d < 3; d++) {
			if (neighborVals[d] < LARGE_VAL && neighborVals[d] < this->grid->field[gridPos]) {
				uMins.push_back(neighborVals[d]);
			}
			else {
				adim--;
			}
		}
		if (adim == 0) return LARGE_VAL;
		std::sort(uMins.begin(), uMins.end());
		float sol;
		for (uint i = 1; i <= adim; i++) {
			sol = this->EikonalSolveInDim(uMins, i);
			if (i == adim || fabs(sol - uMins[i]) < 1e5 * FLT_EPSILON) break;
		}
		return sol;
	};


	for (int iz = iZmin; iz != iZmax; iz += dir[2]) {
		for (int iy = iYmin; iy != iYmax; iy += dir[1]) {
			for (int ix = iXmin; ix != iXmax; ix += dir[0]) {
				gridPos = Nx * Ny * iz + Nx * iy + ix;
				if (!this->grid->frozenCells[gridPos]) {
					u = this->grid->field[gridPos];
					uNew = SolveEikonal(ix, iy, iz);
					if (uNew + 1e5 * FLT_EPSILON < u) {
						this->grid->field[gridPos] = uNew;
						grid->max = ((uNew < LARGE_VAL && uNew > u) ? uNew : grid->max);
					}
				}
			}
		}
	}
}
