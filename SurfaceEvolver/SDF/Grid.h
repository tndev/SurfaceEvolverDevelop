#ifndef GRID_H_
#define GRID_H_

#include <vector>
#include <fstream>
#include <algorithm>
#include "../BVH/AABBTree.h"
#include "../Geometry/Vector3.h"
#include "../Geometry/Box3.h"

#define uint unsigned int

#define LARGE_VAL 10000000.0f

class Grid
{
public:
	std::vector<float> field;
	std::vector<bool> frozenCells;
	uint Nx = 0, Ny = 0, Nz = 0; // index dims
	Vector3 scale = Vector3(1.0f, 1.0f, 1.0f);
	Box3 bbox;

	float min = 0.0f;
	float max = 100.0f;

	Grid();
	Grid(const Grid& other);
	Grid(uint Nx, uint Ny, uint Nz, Box3 bbox, bool addOffset = true, float initVal = LARGE_VAL);
	~Grid();

	void exportToVTI(std::string filename);

	void initToVal(float val);
	void blur();
	bool equalInDimTo(Grid& other);
	void add(Grid& other);
	void sub(Grid& other);
	void absField();
	void computeSignField(AABBTree* v_aabb, AABBTree* e_aabb, AABBTree* t_aabb);
	void bruteForceDistanceField(Geometry* geom);
	void aabbDistanceField(AABBTree* aabb);
	void clean();
	void scaleBy(Vector3& s);
	// fraction of the scale with which the grid should exceed the mesh bbox
	float max_offset_factor = 0.25;
	void getSurroundingCells(Vector3& pos,
		uint oldNx, uint oldNy, uint oldNz, std::vector<float>* oldField,
		std::vector<Vector3>* positionBuffer, std::vector<float>* valueBuffer);

	float getL2Norm();
};

Grid subGrids(Grid g0, Grid g1);
Grid absGrid(Grid g);

float trilinearInterpolate(Vector3& P, std::vector<Vector3>& X, std::vector<float>& f);

#endif
