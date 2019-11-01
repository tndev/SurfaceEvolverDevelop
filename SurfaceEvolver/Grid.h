#ifndef GRID_H_
#define GRID_H_

#include<vector>

#define uint unsigned int

class Grid
{
public:
	std::vector<float>* field;
	uint dimX, dimY, dimZ;
};

#endif
