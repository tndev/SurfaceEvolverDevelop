#ifndef SDF_H_
#define SDF_H_

#include "Grid.h"
#include "FastSweep3D.h"
#include "../BVH/Octree.h"
#include "../ExportImport/VTKExporter.h"

class SDF
{
public:
	Geometry* geom;
	AABBTree* aabb;
	Octree* octree;
	Grid* grid;
	FastSweep3D* fastSweep;

	uint resolution;

	std::string geom_properties = "";
	std::string time_log = "";

	SDF();
	~SDF();
	SDF(const SDF& other);
	SDF(Geometry* geom, uint resolution);

	void exportGrid(VTKExporter* e, std::string export_name = "");
	std::string getComputationProperties();
};

#endif