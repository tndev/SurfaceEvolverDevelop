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
	tri_aabb = other.tri_aabb;
	vert_aabb = other.vert_aabb;
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

	this->tri_aabb = new AABBTree(this->geom);

	auto endSDF_AABB = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedSDF_AABB = (endSDF_AABB - startSDF_AABB);

	auto startSDF_Octree = std::chrono::high_resolution_clock::now();

	this->octree = new Octree(this->tri_aabb, this->tri_aabb->bbox, resolution);

	auto endSDF_Octree = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedSDF_Octree = (endSDF_Octree - startSDF_Octree);

	auto startSDF_FS = std::chrono::high_resolution_clock::now();

	this->grid = new Grid(resolution, resolution, resolution, this->octree->cubeBox);
	this->octree->setLeafValueToScalarGrid(this->grid);
	this->fastSweep = new FastSweep3D(this->grid, 8);

	auto endSDF_FS = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedSDF_FS = (endSDF_FS - startSDF_FS);

	// TODO: Sign computation
	/*
	This approach requires one to find the closest point on a mesh feature (vertex, edge, triangle).
	After ditching the FastSweep3D more viable methods would work something like this:

	Vector3 v = v_aabb->getClosestVertexToPoint(p);
	Vector3 ep = e_aabb->getClosestEdgePointToPoint(p);
	Vector3 tp = t_aabb->getClosestTrianglePointToPoint(p);

	and then pick the one closest to p.
	This method would be much slower than the original DF computation using just the Triangle AABB, Octree and FastSweep3D.
	Since it would already produce distances, this approach would not need to use an Eikonal solver such as FastSweep.
	*/

	/*
	auto startSDF_VertTree = std::chrono::high_resolution_clock::now();
	this->vert_aabb = new AABBTree(geom, PrimitiveType::vert);
	auto endSDF_VertTree = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedSDF_VertTree = (endSDF_VertTree - startSDF_VertTree);

	auto startSDF_EdgeTree = std::chrono::high_resolution_clock::now();
	this->edge_aabb = new AABBTree(geom, PrimitiveType::edge);
	auto endSDF_EdgeTree = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedSDF_EdgeTree = (endSDF_EdgeTree - startSDF_EdgeTree);

	auto startSDF_Sign = std::chrono::high_resolution_clock::now();
	this->grid->computeSignField(this->vert_aabb, this->edge_aabb, this->tri_aabb);
	auto endSDF_Sign = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedSDF_Sign = (endSDF_Sign - startSDF_Sign);
	*/

	auto endSDF = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsedSDF = (endSDF - startSDF);

	this->geom_properties = "=== " + this->geom->name + " === \n" + "verts: " + std::to_string(this->geom->uniqueVertices.size()) +
		", triangles: " + std::to_string(this->tri_aabb->primitives.size()) + ", grid resolution: " + std::to_string(resolution) + "^3 \n";
	this->time_log = 
		"computation times:  AABBTree: " + std::to_string(elapsedSDF_AABB.count()) +
		" s, Octree: (build: " + std::to_string(elapsedSDF_Octree.count()) + " s, get_leaves: " + std::to_string(this->octree->leaf_retrieve_time) + 
		" s), FastSweep3D: " + std::to_string(elapsedSDF_FS.count()) + " s, \n" +
		//" vertex AABBTree: " + std::to_string(elapsedSDF_VertTree.count()) + " s, edge AABBTree: " + std::to_string(elapsedSDF_EdgeTree.count()) + " s, grid sign computation: " + std::to_string(elapsedSDF_Sign.count()) + "s\n" +
		"====> TOTAL: " + std::to_string(elapsedSDF.count()) + " s" + "\n\n";
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
