#include "SDF.h"

SDF::SDF()
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

SDF::SDF(Geometry* geom, uint resolution, bool computeSign, bool computeGradient, 
	bool saveGridStates, bool scaleAndInterpolate, SDF_Method method)
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

		std::chrono::duration<float> elapsedSDF_Sign;

		if (scaleAndInterpolate && resolution > this->resolution_limit) {
			auto startSDF_FS = std::chrono::high_resolution_clock::now();

			this->grid = new Grid(this->resolution_limit, this->resolution_limit, this->resolution_limit, this->octree->bbox, this->octree->cubeBox);
			this->octree->setLeafValueToScalarGrid(this->grid);
			this->grid->expand();
			this->fastSweep = new FastSweep3D(this->grid, 8, saveGridStates);

			auto endSDF_FS = std::chrono::high_resolution_clock::now();
			elapsedSDF_FS = (endSDF_FS - startSDF_FS);

			if (computeSign) {
				auto startSDF_Sign = std::chrono::high_resolution_clock::now();
				this->grid->computeSignField(this->tri_aabb);
				auto endSDF_Sign = std::chrono::high_resolution_clock::now();
				elapsedSDF_Sign = (endSDF_Sign - startSDF_Sign);
			}

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
			
			VTKExporter exporter = VTKExporter();

			this->grid = new Grid(resolution, resolution, resolution, this->octree->bbox, this->octree->cubeBox);
			this->octree->setLeafValueToScalarGrid(this->grid);
			
			this->grid->expand();
			this->fastSweep = new FastSweep3D(this->grid, 8, saveGridStates);

			auto endSDF_FS = std::chrono::high_resolution_clock::now();
			elapsedSDF_FS = (endSDF_FS - startSDF_FS);

			if (computeSign) {
				auto startSDF_Sign = std::chrono::high_resolution_clock::now();
				this->grid->computeSignField(this->tri_aabb);
				auto endSDF_Sign = std::chrono::high_resolution_clock::now();
				elapsedSDF_Sign = (endSDF_Sign - startSDF_Sign);
			}
		}

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
			(computeSign ? ("grid sign computation: " + std::to_string(elapsedSDF_Sign.count()) + "s\n") : "") +
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

		this->grid = new Grid(resolution, resolution, resolution, bbox, cubeBox);
		auto startSDF_Lookup = std::chrono::high_resolution_clock::now();
		this->grid->expand();
		this->grid->aabbDistanceField(this->tri_aabb);
		auto endSDF_Lookup = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedSDF_Lookup = (endSDF_Lookup - startSDF_Lookup);

		auto endSDF = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedSDF = (endSDF - startSDF);

		this->geom_properties = "=== " + this->geom->name + " === \n" + "verts: " + std::to_string(this->geom->uniqueVertices.size()) +
			", triangles: " + std::to_string(this->tri_aabb->primitives.size()) + ", grid resolution: " + std::to_string(this->grid->Nx) + "^3 \n";
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

		this->grid = new Grid(resolution, resolution, resolution, bbox, cubeBox);
		this->grid->expand();
		this->grid->bruteForceDistanceField(geom);

		auto endSDF = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedSDF = (endSDF - startSDF);

		this->geom_properties = "=== " + this->geom->name + " === \n" + "verts: " + std::to_string(this->geom->uniqueVertices.size()) +
			", triangles: " + std::to_string(this->geom->vertexIndices.size() / 3) + ", grid resolution: " + std::to_string(this->grid->Nx) + "^3 \n";
		this->time_log =
			"BRUTE FORCE ::: computation times: ====> TOTAL: " + std::to_string(elapsedSDF.count()) + " s" + "\n\n";
	}

	if (computeGradient) {
		this->grid->computeGradient();
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

void SDF::exportGradientField(VTKExporter* e, std::string export_name)
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
		this->grid->exportGradientToVTK("voxGradSDF" + std::to_string(this->resolution) + method_name + std::to_string(this->resolution));
	}
	else {
		this->grid->exportGradientToVTK(export_name);
	}
}

std::string SDF::getComputationProperties()
{
	return this->geom_properties + this->time_log;
}

void SDF::applyMatrix(Matrix4& m)
{
	if (!this->tri_aabb && !this->octree) {
		std::cout << "SDF not initiated" << std::endl;
		return;
	}

	auto startSDFTransform = std::chrono::high_resolution_clock::now();
	this->geom->applyMatrix(m);
	Matrix4 mInverse = Matrix4().getInverse(m);
	this->tri_aabb->applyMatrix(m);
	auto endSDFTransform = std::chrono::high_resolution_clock::now();

	auto startSDFRecompute = std::chrono::high_resolution_clock::now();

	delete this->octree;
	this->octree = new Octree(this->tri_aabb, this->tri_aabb->bbox, resolution);

	delete this->grid;
	this->grid = new Grid(resolution, resolution, resolution, this->octree->bbox, this->octree->cubeBox);

	this->octree->setLeafValueToScalarGrid(this->grid);
	
	delete this->fastSweep;
	this->fastSweep = new FastSweep3D(this->grid, 8);

	auto endSDFRecompute = std::chrono::high_resolution_clock::now();

	std::chrono::duration<float> elapsedSDF_Transform = (endSDFTransform - startSDFTransform);
	std::chrono::duration<float> elapsedSDF_Recompute = (endSDFRecompute - startSDFRecompute);

	this->last_transform = "SDF Transform: " + std::to_string(elapsedSDF_Transform.count()) +
		" s, SDF Recompute: " + std::to_string(elapsedSDF_Recompute.count()) +
		" s, TOTAL: " + std::to_string(elapsedSDF_Recompute.count() + elapsedSDF_Transform.count()) + " s\n";
}
