#include "Grid.h"

Grid::Grid()
{
}

Grid::Grid(const Grid& other)
{
	this->field = std::vector<float>(other.field);
	this->frozenCells = std::vector<bool>(other.frozenCells);
	this->Nx = other.Nx; this->Ny = other.Ny; this->Nz = other.Nz;
	this->scale = other.scale;
	this->bbox = other.bbox;

	this->min = other.min;
	this->max = other.max;
}

Grid::Grid(uint Nx, uint Ny, uint Nz, Box3 bbox, float initVal)
{
	// re-scale to fit the rest of the field
	this->bbox = bbox;
	Vector3 origScale = bbox.getSize();
	float offset = max_offset_factor * std::fmaxf(origScale.x, std::fmaxf(origScale.y, origScale.z));
	this->bbox.expandByOffset(offset);
	Vector3 newScale = this->bbox.getSize();
	this->Nx = (uint)std::floor(newScale.x / origScale.x * Nx); 
	this->Ny = (uint)std::floor(newScale.y / origScale.y * Ny);
	this->Nz = (uint)std::floor(newScale.z / origScale.z * Nz);
	this->scale = newScale;
	this->field = std::vector<float>((size_t)this->Nx * this->Ny * this->Nz, initVal); // init field
	this->frozenCells = std::vector<bool>((size_t)this->Nx * this->Ny * this->Nz, false); // unfreeze all
}

Grid::~Grid()
{
}

void Grid::exportToVTI(std::string filename)
{
	std::fstream vti(filename + ".vti", std::fstream::out);

	Vector3 o = bbox.min; // origin
	uint nx = Nx - 1;
	uint ny = Ny - 1;
	uint nz = Nz - 1;
	float dx = scale.x / nx;
	float dy = scale.y / ny;
	float dz = scale.z / nz;

	vti << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
	vti << "	<ImageData WholeExtent=\"0 " << nx << " 0 " << ny << " 0 " << nz << "\" Origin=\"" << o.x << " " << o.y << " " << o.z << "\" Spacing=\"" << dx << " " << dy << " " << dz << "\">" << std::endl;
	vti << "		<Piece Extent=\"0 " << nx << " 0 " << ny << " 0 " << nz << "\">" << std::endl;
	vti << "			<PointData Scalars=\"Scalars_\">" << std::endl;
	vti << "				<DataArray type=\"Float32\" Name=\"Scalars_\" format=\"ascii\" RangeMin=\"" << min << "\" RangeMax=\"" << max << "\">" << std::endl;

	for (uint i = 0; i < field.size(); i ++) {
		vti << this->field[i] << std::endl;
	}

	vti << "				</DataArray>" << std::endl;
	vti << "			</PointData>" << std::endl;
	vti << "		<CellData>" << std::endl;
	vti << "		</CellData>" << std::endl;
	vti << "	</Piece>" << std::endl;
	vti << "	</ImageData>" << std::endl;
	vti << "</VTKFile>";

	vti.close();
}

void Grid::initToVal(float val)
{
	this->field = std::vector<float>((size_t)Nx * Ny * Nz, val);
}

void Grid::blur()
{
	auto applyGrid3x3Kernel = [&](std::vector<float>* bfield, uint ix, uint iy, uint iz) {
		float kernelSum = 0.0f;
		uint gridPos = 0;
		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				for (int k = -1; k <= 1; k++) {
					gridPos =
						Nx * Ny * std::min(std::max((int)iz + i, 0), (int)Nz - 1) +
						Nx * std::min(std::max((int)iy + j, 0), (int)Ny - 1) +
						std::min(std::max((int)ix + k, 0), (int)Nx - 1);
					kernelSum += bfield->at(gridPos);
				}
			}
		}
		gridPos = Nx * Ny * iz + Nx * iy + ix;
		float result = kernelSum / 27.0f;
		field[gridPos] = result;
	};

	std::vector<float> bfield = std::vector<float>(field);
	for (uint iz = 1; iz < Nz - 1; iz++) {
		for (uint iy = 1; iy < Ny - 1; iy++) {
			for (uint ix = 1; ix < Nx - 1; ix++) applyGrid3x3Kernel(&bfield, ix, iy, iz);
		}
	}
}

bool Grid::equalInDimTo(Grid& other)
{
	return (
		(this->Nx == other.Nx && this->Ny == other.Ny && this->Nz == other.Nz) // &&
		// this->bbox.equals(other.bbox) // this seems like an unnecessary criterion
	);
}

void Grid::add(Grid& other)
{
	if (this->equalInDimTo(other)) {
		for (uint i = 0; i < this->field.size(); i++) {
			this->field[i] += other.field[i];
		}
	}
	else {
		std::cout << "invalid grid dimensions" << std::endl;
	}
}

void Grid::sub(Grid& other)
{
	if (this->equalInDimTo(other)) {
		for (uint i = 0; i < this->field.size(); i++) {
			this->field[i] -= other.field[i];
		}
	}
	else {
		std::cout << "invalid grid dimensions" << std::endl;
	}
}

void Grid::absField()
{
	for (uint i = 0; i < this->field.size(); i++) {
		this->field[i] = fabs(this->field[i]);
	}
}

void Grid::computeSignField(AABBTree* v_aabb, AABBTree* e_aabb, AABBTree* t_aabb)
{
	Vector3 p = Vector3();
	int sign = 1; uint gridPos;
	int vId, eId, tId;

	Vector3 o = bbox.min; // origin
	Vector3 pn; float rv, re, rt;
	uint nx = Nx - 1;
	uint ny = Ny - 1;
	uint nz = Nz - 1;
	float dx = scale.x / nx;
	float dy = scale.y / ny;
	float dz = scale.z / nz;

	std::vector<Vector3> v_awpn = v_aabb->geom->getAngleWeightedVertexPseudoNormals();
	std::vector<Vector3> e_awpn = e_aabb->geom->getAngleWeightedEdgePseudoNormals();
	std::vector<Vector3> t_n = t_aabb->geom->getTriangleNormals();

	for (uint iz = 0; iz < Nz; iz++) {
		for (uint iy = 0; iy < Ny; iy++) {
			for (uint ix = 0; ix < Nx; ix++) {
				p.set(
					o.x + ix * dx,
					o.y + iy * dy,
					o.z + iz * dz
				);
				vId = v_aabb->getClosestPrimitiveId(p);
				if (vId < 0) {
					continue;
				}

				rv = (p - v_aabb->geom->uniqueVertices[vId]).lengthSq();

				eId = e_aabb->getClosestPrimitiveId(p);
				if (eId < 0) {
					continue;
				}

				re = getDistanceToAnEdgeSq(&e_aabb->primitives[eId].vertices, p);
				
				if (rv < re) {
					sign = (dot((p - v_aabb->geom->uniqueVertices[vId]), v_awpn[vId]) < 0.0f ? -1 : 1);
					gridPos = Nx * Ny * iz + Nx * iy + ix;
					this->field[gridPos] *= sign;
					continue;
				}
				else {
					tId = t_aabb->getClosestPrimitiveId(p);
					if (tId < 0) {
						continue;
					}

					rt = getDistanceToATriangleSq(&t_aabb->primitives[tId].vertices, p);

					if (re < rt) {
						sign = (dot((p - getClosestPtOnAnEdge(&e_aabb->primitives[eId].vertices, p)), e_awpn[eId]) < 0.0f ? -1 : 1);
						gridPos = Nx * Ny * iz + Nx * iy + ix;
						this->field[gridPos] *= sign;
						continue;
					}
					else {
						sign = (dot((p - getClosestPtOnATriangle(&t_aabb->primitives[tId].vertices, p)), t_n[tId]) < 0.0f ? -1 : 1);
						gridPos = Nx * Ny * iz + Nx * iy + ix;
						this->field[gridPos] *= sign;
					}
				}
			}
		}
	}
}

void Grid::bruteForceDistanceField(Geometry* geom)
{
	const uint Nx = this->Nx, Ny = this->Ny, Nz = this->Nz;
	int ix, iy, iz, i, gridPos;
	Vector3 p = Vector3(); Tri T;

	Vector3 o = bbox.min; // origin
	uint nx = Nx - 1;
	uint ny = Ny - 1;
	uint nz = Nz - 1;
	float dx = scale.x / nx;
	float dy = scale.y / ny;
	float dz = scale.z / nz;
	float result_distSq, distSq;

	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {

				result_distSq = FLT_MAX;
				p.set(
					o.x + ix * dx,
					o.y + iy * dy,
					o.z + iz * dz
				);

				for (i = 0; i < geom->vertexIndices.size(); i += 3) {
					T = { 
						&geom->uniqueVertices[geom->vertexIndices[i]],
						&geom->uniqueVertices[geom->vertexIndices[i + 1]],
						&geom->uniqueVertices[geom->vertexIndices[i + 2]] 
					};
					distSq = getDistanceToATriangleSq(&T, p);

					result_distSq = distSq < result_distSq ? distSq : result_distSq;
				}

				gridPos = Nx * Ny * iz + Nx * iy + ix;
				this->field[gridPos] = sqrt(result_distSq);
			}
		}
	}
}

void Grid::aabbDistanceField(AABBTree* aabb)
{
	const uint Nx = this->Nx, Ny = this->Ny, Nz = this->Nz;
	int ix, iy, iz, i, gridPos;
	Vector3 p = Vector3(); Tri T;

	Vector3 o = bbox.min; // origin
	uint nx = Nx - 1;
	uint ny = Ny - 1;
	uint nz = Nz - 1;
	float dx = scale.x / nx;
	float dy = scale.y / ny;
	float dz = scale.z / nz;
	float distSq;

	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {

				p.set(
					o.x + ix * dx,
					o.y + iy * dy,
					o.z + iz * dz
				);

				i = aabb->getClosestPrimitiveId(p);
				if (i < 0) continue;

				distSq = getDistanceToAPrimitiveSq(aabb->primitives[i], p);

				gridPos = Nx * Ny * iz + Nx * iy + ix;
				this->field[gridPos] = sqrt(distSq);
			}
		}
	}
}

void Grid::clean()
{
	Nx = 0, Ny = 0, Nz = 0;
	scale = Vector3(1.0f, 1.0f, 1.0f);
	field.clear();
}

Grid subGrids(Grid g0, Grid g1)
{
	Grid result = g0;
	result.sub(g1);
	return result;
}

Grid absGrid(Grid g)
{
	g.absField();
	return g;
}
