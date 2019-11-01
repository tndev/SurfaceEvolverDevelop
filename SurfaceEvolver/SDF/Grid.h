#ifndef GRID_H_
#define GRID_H_

#include <vector>
#include <fstream>
#include "../Geometry/Vector3.h"
#include "../Geometry/Box3.h"

#define uint unsigned int

class Grid
{
public:
	std::vector<float> field;
	uint Nx, Ny, Nz; // index dims
	Vector3 scale = Vector3(1.0f, 1.0f, 1.0f);
	Box3 bbox;

	Grid();
	Grid(uint Nx, uint Ny, uint Nz, Box3 bbox);
	~Grid();

	void exportToRawBinary(std::string filename);
};

#endif
