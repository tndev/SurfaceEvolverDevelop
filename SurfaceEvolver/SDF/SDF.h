#ifndef SDF_H_
#define SDF_H_

#include "Grid.h"
#include "FastSweep3D.h"
#include "../BVH/Octree.h"
#include "../ExportImport/VTKExporter.h"

enum class SDF_Method {
	fast_sweeping = 0,
	aabb_dist = 1,
	brute_force = 2,
};

class SDF
{
public:
	Geometry* geom;
	AABBTree* tri_aabb;
	AABBTree* edge_aabb;
	AABBTree* vert_aabb;
	Octree* octree;
	Grid* grid;
	FastSweep3D* fastSweep;

	uint resolution;
	SDF_Method method = SDF_Method::fast_sweeping;

	std::string geom_properties = "";
	std::string time_log = "";

	SDF();
	~SDF();
	SDF(const SDF& other);
	SDF(Geometry* geom, uint resolution, bool saveGridStates = false, bool scaleAndInterpolate = false, SDF_Method method = SDF_Method::fast_sweeping);

	void exportGrid(VTKExporter* e, std::string export_name = "");
	std::string getComputationProperties();
public:
	uint resolution_limit = 15;
};

#endif