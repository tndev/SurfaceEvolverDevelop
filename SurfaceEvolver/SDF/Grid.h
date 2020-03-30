#ifndef GRID_H_
#define GRID_H_

#include <vector>
#include <stack>
#include <fstream>
#include <algorithm>
#include <stack>
#include <nmmintrin.h>
#include "../BVH/AABBTree.h"
#include "../Geometry/Vector3.h"
#include "../Geometry/Box3.h"

#define uint unsigned int

#define LARGE_VAL 10000000.0f

class Grid
{
public:
	float* field = nullptr;
	bool* frozenCells = nullptr;
	// gradient dims are 2 less than field dims
	// because we're using central differences for derivatives
	float* gradFieldX = nullptr;
	float* gradFieldY = nullptr;
	float* gradFieldZ = nullptr;

	uint gridExtent = 0; uint gradExtent = 0;
	uint Nx, Ny, Nz; // index dims
	Vector3 scale = Vector3(1.0f, 1.0f, 1.0f);
	Box3 bbox;
	Box3 cubeBox;

	float min = 0.0f;
	float max = 100.0f;

	Grid();
	Grid(const Grid& other);
	Grid(uint Nx, uint Ny, uint Nz, Box3 bbox, Box3 cubeBox, float initVal = LARGE_VAL);
	~Grid();
	bool equalInDimTo(Grid& other);

	void exportToVTI(std::string filename);
	void exportGradientToVTK(std::string filename);

	void initToVal(float val);
	void blur();
	

	void add(Grid& other);
	void sub(Grid& other);
	void absField();
	void negate();
	void computeSignField(AABBTree* aabb);
	void computeGradient();
	void expand(float initVal = LARGE_VAL);
	void clip(Box3& targetBox);

	bool hasGradient();

	void bruteForceDistanceField(Geometry* geom);
	void aabbDistanceField(AABBTree* aabb);
	void cleanField();
	void cleanGrad();

	void scaleBy(Vector3& s);

	// fraction of the scale with which the grid should exceed the mesh bbox
	float max_offset_factor = 1.0f;
	void getSurroundingCells(Vector3& pos,
		uint oldNx, uint oldNy, uint oldNz, float* oldField,
		std::vector<Vector3>* positionBuffer, std::vector<float>* valueBuffer);

	float getL2Norm();

	void clearFrozenCells();
	void clearField();
	Vector3 grad(Vector3& p, Vector3& dXYZ, std::vector<Vector3>* positionBuffer, std::vector<float>* valueBuffer);
};

Grid subGrids(Grid g0, Grid g1);
Grid absGrid(Grid g);

float trilinearInterpolate(Vector3& P, std::vector<Vector3>& X, std::vector<float>& f);

#endif
