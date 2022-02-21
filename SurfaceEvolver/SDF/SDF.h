#ifndef SDF_H_
#define SDF_H_

#include "Grid.h"
#include "FastSweep3D.h"
#include "../BVH/Octree.h"
#include "../ExportImport/VTKExporter.h"

struct SDFTimeLog
{
	double AABBTreeTime;
	double OctreeTime;
	double FastSweepingTime;
	double SignFloodFillTime;
	double TotalTime;
	uint GridRes;

	const SDFTimeLog& operator+= (const SDFTimeLog& other)
	{
		AABBTreeTime += other.AABBTreeTime;
		OctreeTime += other.OctreeTime;
		FastSweepingTime += other.FastSweepingTime;
		SignFloodFillTime += other.SignFloodFillTime;
		TotalTime += other.TotalTime;
		return *this;
	}

	const SDFTimeLog& operator/= (const double value)
	{
		AABBTreeTime /= value;
		OctreeTime /= value;
		FastSweepingTime /= value;
		SignFloodFillTime /= value;
		TotalTime /= value;
		return *this;
	}
};

enum class SDF_Method {
	fast_sweeping = 0,
	aabb_dist = 1,
	brute_force = 2,
};

class SDF
{
public:
	std::shared_ptr<Geometry> geom = nullptr;
	std::shared_ptr<AABBTree> tri_aabb = nullptr;
	std::shared_ptr<AABBTree> edge_aabb = nullptr;
	std::shared_ptr<AABBTree> vert_aabb = nullptr;
	std::shared_ptr<Octree> octree = nullptr;
	std::shared_ptr<Grid> grid = nullptr; // contains scalar field (& its gradient)
	std::shared_ptr<FastSweep3D> fastSweep = nullptr;

	uint resolution;
	SDF_Method method = SDF_Method::fast_sweeping;

	std::string geom_properties = "";
	std::string time_log = "";
	std::string last_transform = "";

	std::string pathPrefix;

	SDFTimeLog timeLog;

	SDF();
	SDF(const SDF& other);
	SDF(const Geometry& geom, uint resolution, std::string path, bool computeSign = false, bool computeGradient = true,
		bool saveGridStates = false, bool scaleAndInterpolate = false, SDF_Method method = SDF_Method::fast_sweeping);
	~SDF() = default;

	void exportGrid(VTKExporter* e, std::string export_name = "");
	void exportGradientField(VTKExporter* e, std::string export_name = "");
	std::string getComputationProperties();

	void applyMatrix(Matrix4& m);
public:
	uint resolution_limit = 20;
};

#endif