#include "FastSweep3D.h"

FastSweep3D::FastSweep3D()
{
}

FastSweep3D::FastSweep3D(const FastSweep3D& other)
{
	grid = other.grid;
	f = other.f;
	h = other.h;
	Nsweeps = other.Nsweeps;
}

FastSweep3D::FastSweep3D(Grid* grid, uint Nsweeps, bool saveGridStates, bool blur)
{
	this->grid = grid;
	this->h = this->grid->scale.x / this->grid->Nx;
	this->Nsweeps = Nsweeps;

	this->sweep(saveGridStates);
	if (blur) this->grid->blur();
}

FastSweep3D::~FastSweep3D()
{
}

void FastSweep3D::sweep(bool saveGridStates)
{
	const uint Nx = grid->Nx, Ny = grid->Ny, Nz = grid->Nz;
	int s, i, j, k, gridPos;
	float aa[3], tmp, eps = 1e-6;
	float d_curr, d_new, a, b, c, D;

	// sweep directions { start, end, step }
	const int dirX[8][3] = { { 0, Nx - 1, 1 }, { Nx - 1, 0, -1 }, { Nx - 1, 0, -1 }, { Nx - 1, 0, -1 }, { Nx - 1, 0, -1 }, { 0, Nx - 1, 1 }, { 0, Nx - 1, 1 }, { 0, Nx - 1, 1 } };
	const int dirY[8][3] = { { 0, Ny - 1, 1 }, { 0, Ny - 1, 1 }, { Ny - 1, 0, -1 }, { Ny - 1, 0, -1 }, { 0, Ny - 1, 1 }, { 0, Ny - 1, 1 }, { Ny - 1, 0, -1 }, { Ny - 1, 0, -1 } };
	const int dirZ[8][3] = { { 0, Nz - 1, 1 }, { 0, Nz - 1, 1 }, { 0, Nz - 1, 1 }, { Nz - 1, 0, -1 }, { Nz - 1, 0, -1 }, { Nz - 1, 0, -1 }, { Nz - 1, 0, -1 }, { 0, Nz - 1, 1 } };

	if (saveGridStates) { // save grid state
		this->grid->exportToVTI("FSGrid_sweep-0");
	}

	for (s = 0; s < Nsweeps; s++) {
		// std::cout << "sweep " << s << " ... " << std::endl;
		for (i = dirX[s][0]; dirX[s][2] * i <= dirX[s][1]; i += dirX[s][2]) {
			for (j = dirY[s][0]; dirY[s][2] * j <= dirY[s][1]; j += dirY[s][2]) {
				for (k = dirZ[s][0]; dirZ[s][2] * k <= dirZ[s][1]; k += dirZ[s][2]) {
					gridPos = ((i * Ny + j) * Nz + k);
					if (!grid->frozenCells[gridPos]) {

						// === neighboring cells (Upwind Godunov) ===
						// using ternary operator ( ? : ) is faster most of the time
						if (i == 0 || i == (Nx - 1)) {
							if (i == 0) {
								aa[0] = grid->field[gridPos] < grid->field[(((i + 1) * Ny + j) * Nz + k)] ? grid->field[gridPos] : grid->field[(((i + 1) * Ny + j) * Nz + k)];
							}
							if (i == (Nx - 1)) {
								aa[0] = grid->field[(((i - 1) * Ny + j) * Nz + k)] < grid->field[gridPos] ? grid->field[(((i - 1) * Ny + j) * Nz + k)] : grid->field[gridPos];
							}
						}
						else {
							aa[0] = grid->field[(((i - 1) * Ny + j) * Nz + k)] < grid->field[(((i + 1) * Ny + j) * Nz + k)] ? grid->field[(((i - 1) * Ny + j) * Nz + k)] : grid->field[(((i + 1) * Ny + j) * Nz + k)];
						}

						if (j == 0 || j == (Ny - 1)) {
							if (j == 0) {
								aa[1] = grid->field[gridPos] < grid->field[((i * Ny + (j + 1)) * Nz + k)] ? grid->field[gridPos] : grid->field[((i * Ny + (j + 1)) * Nz + k)];
							}
							if (j == (Ny - 1)) {
								aa[1] = grid->field[((i * Ny + (j - 1)) * Nz + k)] < grid->field[gridPos] ? grid->field[((i * Ny + (j - 1)) * Nz + k)] : grid->field[gridPos];
							}
						}
						else {
							aa[1] = grid->field[((i * Ny + (j - 1)) * Nz + k)] < grid->field[((i * Ny + (j + 1)) * Nz + k)] ? grid->field[((i * Ny + (j - 1)) * Nz + k)] : grid->field[((i * Ny + (j + 1)) * Nz + k)];
						}

						if (k == 0 || k == (Nz - 1)) {
							if (k == 0) {
								aa[2] = grid->field[gridPos] < grid->field[((i * Ny + j) * Nz + (k + 1))] ? grid->field[gridPos] : grid->field[((i * Ny + j) * Nz + (k + 1))];
							}
							if (k == (Nz - 1)) {
								aa[2] = grid->field[((i * Ny + j) * Nz + (k - 1))] < grid->field[gridPos] ? grid->field[((i * Ny + j) * Nz + (k - 1))] : grid->field[gridPos];
							}
						}
						else {
							aa[2] = grid->field[((i * Ny + j) * Nz + (k - 1))] < grid->field[((i * Ny + j) * Nz + (k + 1))] ? grid->field[((i * Ny + j) * Nz + (k - 1))] : grid->field[((i * Ny + j) * Nz + (k + 1))];
						}

						// simple bubble sort
						if (aa[0] > aa[1]) { tmp = aa[0]; aa[0] = aa[1]; aa[1] = tmp; }
						if (aa[1] > aa[2]) { tmp = aa[1]; aa[1] = aa[2]; aa[2] = tmp; }
						if (aa[0] > aa[1]) { tmp = aa[0]; aa[0] = aa[1]; aa[1] = tmp; }

						d_curr = aa[0] + h * f;
						if (d_curr <= (aa[1] + eps)) {
							d_new = d_curr;
						}
						else {
							a = 2.0f; b = -2.0f * (aa[0] + aa[1]);
							c = aa[0] * aa[0] + aa[1] * aa[1] - h * h * f * f;
							D = sqrt(b * b - 4.0f * a * c);

							d_curr = ((-b + D) > (-b - D) ? (-b + D) : (-b - D)) / (2.0f * a);

							if (d_curr <= (aa[2] + eps))
								d_new = d_curr;
							else {
								a = 3.0f;
								b = -2.0f * (aa[0] + aa[1] + aa[2]);
								c = aa[0] * aa[0] + aa[1] * aa[1] + aa[2] * aa[2] - h * h * f * f;
								D = sqrt(b * b - 4.0f * a * c);

								d_new = ((-b + D) > (-b - D) ? (-b + D) : (-b - D)) / (2.0f * a);
							}
						}

						grid->field[gridPos] = grid->field[gridPos] < d_new ? grid->field[gridPos] : d_new;
					}
				}
			}
		}

		if (saveGridStates) { // save grid state
			this->grid->exportToVTI("FSGrid_sweep-" + std::to_string(s + 1));
		}
	}
	grid->frozenCells.clear();
}

