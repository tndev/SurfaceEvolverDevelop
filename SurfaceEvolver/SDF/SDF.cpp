#include "SDF.h"

SDF::SDF()
{
}

SDF::~SDF()
{
}

SDF::SDF(const SDF& other)
{
	geom = other.geom;
	aabb = other.aabb;
	octree = other.octree;
	grid = other.grid;
	fastSweep = other.fastSweep;

	resolution = other.resolution;
	geom_properties = other.geom_properties;
	time_log = other.time_log;
}

SDF::SDF(Geometry* geom, uint resolution)
{
	this->resolution = resolution;
	this->geom = geom;

	auto startSDF = std::chrono::high_resolution_clock::now();
	// === Timed code ============

	auto startSDF_AABB = std::chrono::high_resolution_clock::now();

	this->aabb = new AABBTree(this->geom);

	auto endSDF_AABB = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedSDF_AABB = (endSDF_AABB - startSDF_AABB);

	auto startSDF_Octree = std::chrono::high_resolution_clock::now();

	this->octree = new Octree(this->aabb, this->aabb->bbox, resolution);

	auto endSDF_Octree = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedSDF_Octree = (endSDF_Octree - startSDF_Octree);

	auto startSDF_FS = std::chrono::high_resolution_clock::now();

	this->grid = new Grid(resolution, resolution, resolution, this->aabb->bbox);
	this->octree->setLeafValueToScalarGrid(this->grid);
	this->fastSweep = new FastSweep3D(this->grid, 8);

	auto endSDF_FS = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedSDF_FS = (endSDF_FS - startSDF_FS);

	// TODO: Sign computation

	/* 
	auto startSDF_Sign = std::chrono::high_resolution_clock::now();
	this->grid->computeSignField(this->aabb);
	auto endSDF_Sign = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedSDF_Sign = (endSDF_Sign - startSDF_Sign);
	*/

	auto endSDF = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedSDF = (endSDF - startSDF);

	this->geom_properties = "=== " + this->geom->name + " === \n" + "verts: " + std::to_string(this->geom->uniqueVertices.size()) +
		", triangles: " + std::to_string(this->aabb->triangles.size()) + ", grid resolution: " + std::to_string(resolution) + "^3 \n";
	this->time_log = 
		"computation times:  AABB: " + std::to_string(elapsedSDF_AABB.count()) +
		" s, Octree: (build: " + std::to_string(elapsedSDF_Octree.count()) + " s, get_leaves: " + std::to_string(this->octree->leaf_retrieve_time) + 
		" s), FastSweep3D: " + std::to_string(elapsedSDF_FS.count()) +
		// ", Sign: " + std::to_string(elapsedSDF_Sign.count()) +
		" s, TOTAL: " + std::to_string(elapsedSDF.count()) + " s" + "\n\n";
}

void SDF::exportGrid(VTKExporter* e, std::string export_name)
{
	if (export_name.empty()) {
		this->grid->exportToVTI("voxFieldSDF" + std::to_string(this->resolution));
	}
	else {
		this->grid->exportToVTI(export_name);
	}		
}

std::string SDF::getComputationProperties()
{
	return this->geom_properties + this->time_log;
}
