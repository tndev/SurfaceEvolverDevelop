#include "Grid.h"

Grid::Grid()
{
}

Grid::Grid(uint Nx, uint Ny, uint Nz, Box3 bbox)
{
	this->dimX = Nx; this->dimY = Ny; this->dimZ = Nz;
	this->bbox = bbox;
	this->scale = this->bbox.getSize();
}

Grid::~Grid()
{
}
