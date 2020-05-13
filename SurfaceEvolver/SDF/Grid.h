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

#define LARGE_VAL DBL_MAX

class Grid
{
public:
	double* field = nullptr;
	bool* frozenCells = nullptr;
	// gradient dims are 2 less than field dims
	// because we're using central differences for derivatives
	double* gradFieldX = nullptr;
	double* gradFieldY = nullptr;
	double* gradFieldZ = nullptr;

	uint gridExtent = 0; uint gradExtent = 0;
	uint Nx, Ny, Nz; // index dims
	Vector3 scale = Vector3(1.0, 1.0, 1.0);
	Box3 bbox;
	Box3 cubeBox;

	double min = 0.0;
	double max = 100.0;

	Grid();
	Grid(const Grid& other);
	Grid(uint Nx, uint Ny, uint Nz, Box3 bbox, Box3 cubeBox, double initVal = LARGE_VAL);
	~Grid();
	bool equalInDimTo(Grid& other);

	void exportToVTI(std::string filename);
	void exportGradientToVTK(std::string filename);

	void initToVal(double val);
	void blur();
	

	void add(Grid& other);
	void sub(Grid& other);
	void absField();
	void negate();
	void computeSignField(AABBTree* aabb);
	void computeGradient();
	void expand(double initVal = LARGE_VAL);
	void clip(Box3& targetBox);

	bool hasGradient();

	void bruteForceDistanceField(Geometry* geom);
	void aabbDistanceField(AABBTree* aabb);
	void cleanField();
	void cleanGrad();

	void scaleBy(Vector3& s);

	// fraction of the scale with which the grid should exceed the mesh bbox
	double max_offset_factor = 1.0;
	void getSurroundingCells(Vector3& pos,
		uint oldNx, uint oldNy, uint oldNz, double* oldField,
		std::vector<Vector3>* positionBuffer, std::vector<double>* valueBuffer);

	double getL2Norm();

	void clearFrozenCells();
	void clearField();
	Vector3 grad(Vector3& p, Vector3& dXYZ, std::vector<Vector3>* positionBuffer, std::vector<double>* valueBuffer);
};

Grid subGrids(Grid g0, Grid g1);
Grid absGrid(Grid g);

double trilinearInterpolate(Vector3& P, std::vector<Vector3>& X, std::vector<double>& f);

#endif
