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
	// std::vector<uint> frozenCellIds = {}; // ids of cells that are not to be manipulated with
	uint Nx = 0, Ny = 0, Nz = 0; // index dims
	Vector3 scale = Vector3(1.0f, 1.0f, 1.0f);
	Box3 bbox;

	float min = 0.0f;
	float max = 100.0f;

	Grid();
	Grid(const Grid& other);
	Grid(uint Nx, uint Ny, uint Nz, Box3 bbox);
	~Grid();

	void exportToVTI(std::string filename);

	void clean();
private:
	// fraction of the scale with which the grid should exceed the mesh bbox
	float max_offset_factor = 0.25;
};

#endif
