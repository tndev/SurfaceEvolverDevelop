#include "FastSweep3D.h"

FastSweep3D::FastSweep3D()
{
}

FastSweep3D::FastSweep3D(Grid* grid)
{
	this->grid = grid;
	this->h = this->grid->scale.x / this->grid->Nx;

	auto startFastSweep = std::chrono::high_resolution_clock::now();
	// === Timed code ============
	std::cout << "Initiating FastSweep3D..." << std::endl;
	for (uint i = 0; i < this->Nsweeps; i++) {
		std::cout << "sweep " << i << " ... " << std::endl;
		this->sweep(i);
	}
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

void FastSweep3D::sweep(uint dirId)
{
	uint Nx = this->grid->Nx; uint Ny = this->grid->Ny; uint Nz = this->grid->Nz;
	
	uint iXmin = this->sweepDir[dirId][0] > 0 ? 0 : Nx;
	uint iXmax = this->sweepDir[dirId][0] > 0 ? Nx : 0;

	uint iYmin = this->sweepDir[dirId][1] > 0 ? 0 : Ny;
	uint iYmax = this->sweepDir[dirId][1] > 0 ? Ny : 0;

	uint iZmin = this->sweepDir[dirId][2] > 0 ? 0 : Nz;
	uint iZmax = this->sweepDir[dirId][2] > 0 ? Nz : 0;

	std::cout << "sweep ranges: iX = [" << iXmin << ", " << iXmax << ", step = " << this->sweepDir[dirId][0] << "]" << std::endl;
	std::cout << "sweep ranges: iY = [" << iYmin << ", " << iYmax << ", step = " << this->sweepDir[dirId][1] << "]" << std::endl;
	std::cout << "sweep ranges: iZ = [" << iZmin << ", " << iZmax << ", step = " << this->sweepDir[dirId][2] << "]" << std::endl;

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

	for (uint iz = iZmin; iz < iZmax; iz += this->sweepDir[dirId][2]) {
		for (uint iy = iYmin; iy < iYmax; iy += this->sweepDir[dirId][1]) {
			for (uint ix = iXmin; ix < iXmax; ix += this->sweepDir[dirId][0]) {
				gridPos = Nx * Ny * iz + Nx * iy + ix;
				u = this->grid->field[gridPos];
				uNew = SolveEikonal(ix, iy, iz);
				if (uNew + 1e5 * FLT_EPSILON < u) {
					this->grid->field[gridPos] = uNew;
					grid->max = (uNew < LARGE_VAL ? uNew : grid->max);
				}
			}
		}
	}
}
