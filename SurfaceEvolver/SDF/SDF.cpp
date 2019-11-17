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

SDF::SDF(Geometry* geom, uint resolution, bool scaleAndInterpolate, SDF_Method method)
{
	this->resolution = resolution;
	this->geom = geom;
	this->method = method;

	if (method == SDF_Method::fast_sweeping) {
		auto startSDF = std::chrono::high_resolution_clock::now();
		// === Timed code ============

		auto startSDF_AABB = std::chrono::high_resolution_clock::now();

		this->tri_aabb = new AABBTree(this->geom);

		auto endSDF_AABB = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedSDF_AABB = (endSDF_AABB - startSDF_AABB);

		auto startSDF_Octree = std::chrono::high_resolution_clock::now();

		if (scaleAndInterpolate && resolution > this->resolution_limit) {
			this->resolution = this->resolution_limit;
			this->octree = new Octree(this->tri_aabb, this->tri_aabb->bbox, this->resolution_limit);
		}
		else {
			this->octree = new Octree(this->tri_aabb, this->tri_aabb->bbox, resolution);
		}		

		auto endSDF_Octree = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedSDF_Octree = (endSDF_Octree - startSDF_Octree);

		std::chrono::duration<float> elapsedGridScale;
		std::chrono::duration<float> elapsedSDF_FS;

		if (scaleAndInterpolate && resolution > this->resolution_limit) {
			auto startSDF_FS = std::chrono::high_resolution_clock::now();

			this->grid = new Grid(this->resolution_limit, this->resolution_limit, this->resolution_limit, this->octree->cubeBox);
			this->octree->setLeafValueToScalarGrid(this->grid);
			this->fastSweep = new FastSweep3D(this->grid, 8);

			auto endSDF_FS = std::chrono::high_resolution_clock::now();
			elapsedSDF_FS = (endSDF_FS - startSDF_FS);

			float origRes = std::floor((1.0f + 2.0f * this->grid->max_offset_factor) * resolution);

			Vector3 scaleFactor = Vector3(
				origRes / std::floor((1.0f + 2.0f * this->grid->max_offset_factor) * this->resolution_limit),
				origRes / std::floor((1.0f + 2.0f * this->grid->max_offset_factor) * this->resolution_limit),
				origRes / std::floor((1.0f + 2.0f * this->grid->max_offset_factor) * this->resolution_limit)
			);

			auto startGridScale = std::chrono::high_resolution_clock::now();
			this->grid->scaleBy(scaleFactor);
			auto endGridScale = std::chrono::high_resolution_clock::now();
			elapsedGridScale = (endGridScale - startGridScale);
		}
		else {
			auto startSDF_FS = std::chrono::high_resolution_clock::now();

			this->grid = new Grid(resolution, resolution, resolution, this->octree->cubeBox);
			this->octree->setLeafValueToScalarGrid(this->grid);
			this->fastSweep = new FastSweep3D(this->grid, 8);

			auto endSDF_FS = std::chrono::high_resolution_clock::now();
			elapsedSDF_FS = (endSDF_FS - startSDF_FS);
		}

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
			", triangles: " + std::to_string(this->tri_aabb->primitives.size()) + ", octree resolution: " + std::to_string(resolution) + "^3, grid resolution: " + std::to_string(this->grid->Nx) + "^3 " +
			((scaleAndInterpolate && resolution > this->resolution_limit) ? ", rescaled from: " + std::to_string(this->resolution) + "^2\n" : "\n");
		this->time_log =
			"computation times:  AABBTree: " + std::to_string(elapsedSDF_AABB.count()) +
			" s, Octree: (build: " + std::to_string(elapsedSDF_Octree.count()) + " s, get_leaves: " + std::to_string(this->octree->leaf_retrieve_time) +
			" s), FastSweep3D: " + std::to_string(elapsedSDF_FS.count()) + " s, \n" +
			((scaleAndInterpolate && resolution > this->resolution_limit) ? "Grid scale: " + std::to_string(elapsedGridScale.count()) + " s \n": "") +
			//" vertex AABBTree: " + std::to_string(elapsedSDF_VertTree.count()) + " s, edge AABBTree: " + std::to_string(elapsedSDF_EdgeTree.count()) + " s, grid sign computation: " + std::to_string(elapsedSDF_Sign.count()) + "s\n" +
			"====> TOTAL: " + std::to_string(elapsedSDF.count()) + " s" + "\n\n";
	}
	else if (method == SDF_Method::aabb_dist) {

		auto startSDF = std::chrono::high_resolution_clock::now();
		// === Timed code ============

		auto startSDF_AABB = std::chrono::high_resolution_clock::now();

		this->tri_aabb = new AABBTree(this->geom);

		auto endSDF_AABB = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedSDF_AABB = (endSDF_AABB - startSDF_AABB);

		Box3 bbox = geom->getBoundingBox();
		Vector3 size = bbox.getSize();
		float maxDim = std::max({ size.x, size.y, size.z });
		Box3 cubeBox = Box3(bbox.min, bbox.min + Vector3(maxDim, maxDim, maxDim));

		this->grid = new Grid(resolution, resolution, resolution, cubeBox);
		auto startSDF_Lookup = std::chrono::high_resolution_clock::now();
		this->grid->aabbDistanceField(this->tri_aabb);
		auto endSDF_Lookup = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedSDF_Lookup = (endSDF_Lookup - startSDF_Lookup);

		auto endSDF = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedSDF = (endSDF - startSDF);

		this->geom_properties = "=== " + this->geom->name + " === \n" + "verts: " + std::to_string(this->geom->uniqueVertices.size()) +
			", triangles: " + std::to_string(this->tri_aabb->primitives.size()) + ", grid resolution: " + std::to_string(resolution) + "^3 \n";
		this->time_log =
			"computation times:  AABBTree build: " + std::to_string(elapsedSDF_AABB.count()) + " s, DF \w AABBTree query: " + std::to_string(elapsedSDF_Lookup.count()) + " s\n" +
			"====> TOTAL: " + std::to_string(elapsedSDF.count()) + " s" + "\n\n";
	}
	else if (method == SDF_Method::brute_force) {
		auto startSDF = std::chrono::high_resolution_clock::now();
		// === Timed code ============

		Box3 bbox = geom->getBoundingBox();
		Vector3 size = bbox.getSize();
		float maxDim = std::max({ size.x, size.y, size.z });
		Box3 cubeBox = Box3(bbox.min, bbox.min + Vector3(maxDim, maxDim, maxDim));

		this->grid = new Grid(resolution, resolution, resolution, cubeBox);
		this->grid->bruteForceDistanceField(geom);

		auto endSDF = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedSDF = (endSDF - startSDF);

		this->geom_properties = "=== " + this->geom->name + " === \n" + "verts: " + std::to_string(this->geom->uniqueVertices.size()) +
			", triangles: " + std::to_string(this->geom->vertexIndices.size() / 3) + ", grid resolution: " + std::to_string(resolution) + "^3 \n";
		this->time_log =
			"BRUTE FORCE ::: computation times: ====> TOTAL: " + std::to_string(elapsedSDF.count()) + " s" + "\n\n";
	}
}

void SDF::exportGrid(VTKExporter* e, std::string export_name)
{
	if (export_name.empty()) {
		std::string method_name;
		if (this->method == SDF_Method::fast_sweeping) {
			method_name = "_FSM_";
		}
		else if (this->method == SDF_Method::aabb_dist) {
			method_name = "_AABB_";
		}
		else if (this->method == SDF_Method::brute_force) {
			method_name = "_BRUTE_F_";
		}
		this->grid->exportToVTI("voxFieldSDF" + std::to_string(this->resolution) + method_name + std::to_string(this->resolution));
	}
	else {
		this->grid->exportToVTI(export_name);
	}		
}

std::string SDF::getComputationProperties()
{
	return this->geom_properties + this->time_log;
}
