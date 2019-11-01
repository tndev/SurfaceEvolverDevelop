#include "Grid.h"

Grid::Grid()
{
}

Grid::Grid(uint Nx, uint Ny, uint Nz, Box3 bbox)
{
	this->Nx = Nx; this->Ny = Ny; this->Nz = Nz;
	this->bbox = bbox;
	this->scale = this->bbox.getSize();
	this->field = std::vector<float>((size_t)this->Nx * this->Ny * this->Nz); // init field
}

Grid::~Grid()
{
}

void Grid::exportToRawBinary(std::string filename)
{
	size_t fieldSize = (size_t)Nx * Ny * Nz;
	size_t bytes = sizeof(char) * fieldSize;
	char* data = new char[fieldSize];
	for (uint i = 0; i < fieldSize; i++) {
		data[i] = (char)field[i];
	}

	auto raw = std::fstream(filename + ".raw", std::ios::out | std::ios::binary);
	raw.write((char*)&data[0], bytes);
	raw.close();
}
