#include "Grid.h"


Grid::Grid(const Grid& other)
{	
	this->Nx = other.Nx; this->Ny = other.Ny; this->Nz = other.Nz;

	this->gridExtent = other.gridExtent;
	this->gradExtent = other.gradExtent;

	this->field = other.field;
	this->frozenCells = other.frozenCells;

	this->gradFieldX = other.gradFieldX;
	this->gradFieldY = other.gradFieldY;
	this->gradFieldZ = other.gradFieldZ;

	this->scale = other.scale;
	this->bbox = other.bbox;
	this->cubeBox = other.cubeBox;

	this->min = other.min;
	this->max = other.max;
}

Grid::Grid(uint Nx, uint Ny, uint Nz, const Box3& bbox, const Box3& cubeBox, std::string geomName, double initVal)
{
	this->geomName = geomName;
	this->cubeBox = cubeBox;
	this->bbox = bbox;
	this->scale = this->cubeBox.getSize();
	this->Nx = Nx; this->Ny = Ny; this->Nz = Nz;

	this->gridExtent = this->Nx * this->Ny * this->Nz;
	this->field = std::vector<double>(gridExtent, initVal);
	this->frozenCells = std::vector<bool>(gridExtent, false);
}

void Grid::exportToVTI(std::string filename)
{
	std::fstream vti(geomName + filename + ".vti", std::fstream::out);

	Vector3 o = bbox.min; // origin
	uint nx = Nx - 1;
	uint ny = Ny - 1;
	uint nz = Nz - 1;
	double dx = scale.x / Nx;
	double dy = scale.y / Ny;
	double dz = scale.z / Nz;

	vti << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
	vti << "	<ImageData WholeExtent=\"0 " << nx << " 0 " << ny << " 0 " << nz << "\" Origin=\"" << o.x + 0.5 * dx << " " << o.y + 0.5 * dy << " " << o.z + 0.5 * dz << "\" Spacing=\"" << dx << " " << dy << " " << dz << "\">" << std::endl;
	vti << "		<Piece Extent=\"0 " << nx << " 0 " << ny << " 0 " << nz << "\">" << std::endl;
	vti << "			<PointData Scalars=\"Scalars_\">" << std::endl;
	vti << "				<DataArray type=\"Float32\" Name=\"Scalars_\" format=\"ascii\" RangeMin=\"" << min << "\" RangeMax=\"" << max << "\">" << std::endl;

	for (const auto& val : field) {
		vti << static_cast<float>(val) << std::endl;
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

void Grid::exportGradientToVTK(std::string filename)
{
	std::fstream vtk(filename + ".vtk", std::fstream::out);

	uint nx = Nx - 2;
	uint ny = Ny - 2;
	uint nz = Nz - 2;
	double dx = scale.x / Nx;
	double dy = scale.y / Ny;
	double dz = scale.z / Nz;
	Vector3 o = bbox.min; // origin
	o.x += dx;	o.y += dy; o.z += dz;

	vtk << "# vtk DataFile Version 3.0" << std::endl;
	vtk << "vtk output" << std::endl;
	vtk << "ASCII" << std::endl;
	vtk << "DATASET STRUCTURED_GRID" << std::endl;
	vtk << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;

	vtk << "POINTS " << gradExtent << " double" << std::endl;

	for (uint iz = 0; iz < nz; iz++) {
		for (uint iy = 0; iy < ny; iy++) {
			for (uint ix = 0; ix < nx; ix++) {
				vtk << o.x + (ix + 0.5) * dx << " " << o.y + (iy + 0.5) * dy << " " << o.z + (iz + 0.5) * dz << std::endl;
			}
		}
	}

	vtk << "POINT_DATA " << gradExtent << std::endl;
	vtk << "VECTORS grad double" << std::endl;

	for (uint i = 0; i < gradExtent; i++) {
		vtk << gradFieldX[i] << " " << gradFieldY[i] << " " << gradFieldZ[i] << std::endl;
	}

	vtk.close();
}

void Grid::initToVal(double val)
{
	field.clear();
	frozenCells.clear();
	this->field = std::vector<double>(gridExtent, val);
	this->frozenCells = std::vector<bool>(gridExtent, false);
}

void Grid::blur()
{
	auto applyGrid3x3Kernel = [&](std::vector<double>& bfield, uint ix, uint iy, uint iz) {
		double kernelSum = 0.0;
		uint gridPos = 0;
		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				for (int k = -1; k <= 1; k++) {
					gridPos =
						Nx * Ny * std::min(std::max((int)iz + i, 0), (int)Nz - 1) +
						Nx * std::min(std::max((int)iy + j, 0), (int)Ny - 1) +
						std::min(std::max((int)ix + k, 0), (int)Nx - 1);
					kernelSum += bfield[gridPos];
				}
			}
		}
		gridPos = Nx * Ny * iz + Nx * iy + ix;
		double result = kernelSum / 27.0;
		field[gridPos] = result;
	};

	auto bfield = field;
	for (uint i = 0; i < gridExtent; i++) bfield[i] = this->field[i];
	for (uint iz = 1; iz < Nz - 1; iz++) {
		for (uint iy = 1; iy < Ny - 1; iy++) {
			for (uint ix = 1; ix < Nx - 1; ix++) applyGrid3x3Kernel(bfield, ix, iy, iz);
		}
	}
	field = bfield;
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
		for (uint i = 0; i < gridExtent; i++) {
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
		for (uint i = 0; i < gridExtent; i++) {
			this->field[i] -= other.field[i];
		}
	}
	else {
		std::cout << "invalid grid dimensions" << std::endl;
	}
}

void Grid::absField()
{
	for (auto& elem : field) {
		elem = std::fabs(elem);
	}
}

void Grid::negate()
{
	double val;
	uint gridPos;

	for (gridPos = 0; gridPos < this->gridExtent; gridPos++) {
		val = -1.0 * this->field[gridPos];
		this->field[gridPos] = val;
	}
}

void Grid::computeSignField()
{
	/**/
	this->negate();

	double val; uint gridPos;
	uint nx = Nx - 1;
	uint ny = Ny - 1;
	uint nz = Nz - 1;
	
	uint iz = 0, iy = 0, ix = 0;

	union {	__m128i idsTriple; uint ids[3];	};
	union { __m128i idsMask; uint imask[3]; };
	std::stack<__m128i> stack = {};

	idsTriple = _mm_setr_epi32(ix, iy, iz, INT_MAX);

	// find the first unfrozen cell

	gridPos = 0;
	while (this->frozenCells[gridPos]) {
		idsMask = _mm_cmplt_epi32(idsTriple, _mm_setr_epi32(nx, ny, nz, INT_MAX));
		idsTriple = _mm_add_epi32(
			idsTriple, _mm_setr_epi32(
			(imask[0] > 0) * 1, (imask[1] > 0) * 1, (imask[2] > 0) * 1, 0)
		);
		ix = ids[0]; iy = ids[1]; iz = ids[2];
		gridPos = Nx * Ny * iz + Nx * iy + ix;
	}

	ids[0] = ix; ids[1] = iy; ids[2] = iz;
	stack.push(idsTriple);

	// a simple voxel flood
	while (!stack.empty()) {
		idsTriple = stack.top();
		stack.pop();

		ix = ids[0]; iy = ids[1]; iz = ids[2];
		gridPos = Nx * Ny * iz + Nx * iy + ix;

		if (!this->frozenCells[gridPos]) {
			val = -1.0 * this->field[gridPos];
			this->field[gridPos] = val;
			this->frozenCells[gridPos] = true; // freeze cell when done

			idsMask = _mm_cmpgt_epi32(idsTriple, _mm_set1_epi32(0)); // lower bounds

			if (imask[0] > 0) {
				stack.push(_mm_setr_epi32(ix - 1, iy, iz, INT_MAX));
			}
			if (imask[1] > 0) {
				stack.push(_mm_setr_epi32(ix, iy - 1, iz, INT_MAX));
			}
			if (imask[2] > 0) {
				stack.push(_mm_setr_epi32(ix, iy, iz - 1, INT_MAX));
			}

			idsMask = _mm_cmplt_epi32(idsTriple, _mm_setr_epi32(nx, ny, nz, INT_MAX)); // upper bounds

			if (imask[0] > 0) {
				stack.push(_mm_setr_epi32(ix + 1, iy, iz, INT_MAX));
			}
			if (imask[1] > 0) {
				stack.push(_mm_setr_epi32(ix, iy + 1, iz, INT_MAX));
			}
			if (imask[2] > 0) {
				stack.push(_mm_setr_epi32(ix, iy, iz + 1, INT_MAX));
			}
		}
	}

	/*
	this->negate();

	double val; uint gridPos;
	uint nx = Nx - 1;
	uint ny = Ny - 1;
	uint nz = Nz - 1;

	uint iz = 0, iy = 0, ix = 0;

	std::stack<std::tuple<uint, uint, uint>> stack = {};
	std::tuple<uint, uint, uint> idsTriple;

	// find the first unfrozen cell
	gridPos = 0;
	while (this->frozenCells[gridPos]) {
		ix += (ix < nx ? 1 : 0);
		iy += (iy < ny ? 1 : 0);
		iz += (iz < nz ? 1 : 0);
		gridPos = Nx * Ny * iz + Nx * iy + ix;
	}

	stack.push({ ix, iy, iz });

	// a simple voxel flood
	while (stack.size()) {
		idsTriple = stack.top();
		stack.pop();

		ix = std::get<0>(idsTriple);
		iy = std::get<1>(idsTriple);
		iz = std::get<2>(idsTriple);

		gridPos = Nx * Ny * iz + Nx * iy + ix;

		if (!this->frozenCells[gridPos]) {
			val = -1.0 * this->field[gridPos];
			this->field[gridPos] = val;
			this->frozenCells[gridPos] = true; // freeze cell when done

			if (ix > 0) {
				stack.push({ ix - 1, iy, iz });
			}
			if (ix < nx) {
				stack.push({ ix + 1, iy, iz });
			}
			if (iy > 0) {
				stack.push({ ix, iy - 1, iz });
			}
			if (iy < ny) {
				stack.push({ ix, iy + 1, iz });
			}
			if (iz > 0) {
				stack.push({ ix, iy, iz - 1 });
			}
			if (iz < nz) {
				stack.push({ ix, iy, iz + 1 });
			}
		}
	}*/
}

void Grid::computeGradient()
{
	int ix, iy, iz, i, gradPos;
	
	double grad_f_x, grad_f_y, grad_f_z, norm;
	uint ix0, iy0, iz0, ix1, iy1, iz1;
	uint gridPosPrevX, gridPosNextX;
	uint gridPosPrevY, gridPosNextY;
	uint gridPosPrevZ, gridPosNextZ;
	uint nx = Nx - 1;
	uint ny = Ny - 1;
	uint nz = Nz - 1;
	double dx = scale.x / nx;
	double dy = scale.y / ny;
	double dz = scale.z / nz;
	this->gradExtent = (nx - 1) * (ny - 1) * (nz - 1);

	this->gradFieldX = std::vector<double>(gradExtent); // init vect field
	this->gradFieldY = std::vector<double>(gradExtent);
	this->gradFieldZ = std::vector<double>(gradExtent);

	for (iz = 1; iz < nz; iz++) {
		for (iy = 1; iy < ny; iy++) {
			for (ix = 1; ix < nx; ix++) {
				gradPos = (nx - 1) * (ny - 1) * (iz - 1) + (nx - 1) * (iy - 1) + ix - 1;

				ix0 = ix - 1; iy0 = iy - 1;	iz0 = iz - 1;
				ix1 = ix + 1; iy1 = iy + 1;	iz1 = iz + 1;

				gridPosPrevX = Nx * Ny * iz + Nx * iy + ix0;
				gridPosNextX = Nx * Ny * iz + Nx * iy + ix1;

				gridPosPrevY = Nx * Ny * iz + Nx * iy0 + ix;
				gridPosNextY = Nx * Ny * iz + Nx * iy1 + ix;

				gridPosPrevZ = Nx * Ny * iz0 + Nx * iy + ix;
				gridPosNextZ = Nx * Ny * iz1 + Nx * iy + ix;

				// central difference
				grad_f_x = (field[gridPosNextX] - field[gridPosPrevX]) / (2.0 * dx);
				grad_f_y = (field[gridPosNextY] - field[gridPosPrevY]) / (2.0 * dy);
				grad_f_z = (field[gridPosNextZ] - field[gridPosPrevZ]) / (2.0 * dz);

				// Eikonal gradient |grad(d)|=1 has to be normalized
				norm = sqrt(grad_f_x * grad_f_x + grad_f_y * grad_f_y + grad_f_z * grad_f_z);

				if (norm > 10 * FLT_EPSILON) {
					gradFieldX[gradPos] = grad_f_x / norm;
					gradFieldY[gradPos] = grad_f_y / norm;
					gradFieldZ[gradPos] = grad_f_z / norm;
				}
				else {
					gradFieldX[gradPos] = grad_f_x;
					gradFieldY[gradPos] = grad_f_y;
					gradFieldZ[gradPos] = grad_f_z;
				}
			}
		}
	}
}

void Grid::expand(double initVal)
{
	// re-scale to fit the rest of the field
	Vector3 origScale = bbox.getSize();
	double offset = max_offset_factor * std::fmaxf(origScale.x, std::fmaxf(origScale.y, origScale.z));
	Box3 newBBox = Box3(this->bbox);
	newBBox.expandByOffset(offset);
	Vector3 newScale = newBBox.getSize();

	double cellSize = origScale.x / Nx;
	uint Nx_new = (uint)std::floor(newScale.x / cellSize);
	uint Ny_new = (uint)std::floor(newScale.y / cellSize);
	uint Nz_new = (uint)std::floor(newScale.z / cellSize);

	uint gridExtent_new = Nx_new * Ny_new * Nz_new;

	int ix, iy, iz, i, j, k, gridPosOld, gridPosNew;

	auto newField = std::vector<double>(gridExtent_new);
	auto newFrozenCells = std::vector<bool>(gridExtent_new);

	// initialize larger grid
	for (i = 0; i < gridExtent_new; i++) {
		newField[i] = initVal;
		newFrozenCells[i] = false;
	}

	// find bounds of the old grid:
	uint ix_min = (uint)std::floor(offset * Nx_new / newScale.x);
	uint iy_min = (uint)std::floor(offset * Ny_new / newScale.y);
	uint iz_min = (uint)std::floor(offset * Nz_new / newScale.z);

	uint ix_Max = ix_min + Nx;
	uint iy_Max = iy_min + Ny;
	uint iz_Max = iz_min + Nz;

	i = 0;
	for (iz = iz_min; iz < iz_Max; iz++) {
		j = 0;
		for (iy = iy_min; iy < iy_Max; iy++) {
			k = 0;
			for (ix = ix_min; ix < ix_Max; ix++) {

				gridPosOld = Nx * Ny * i + Nx * j + k;
				gridPosNew = Nx_new * Ny_new * iz + Nx_new * iy + ix;

				newField[gridPosNew] = this->field[gridPosOld];
				newFrozenCells[gridPosNew] = this->frozenCells[gridPosOld];
				k++;
			}
			j++;
		}
		i++;
	}

	this->scale = newScale;

	clearField();
	clearFrozenCells();

	this->Nx = Nx_new; this->Ny = Ny_new; this->Nz = Nz_new;
	this->gridExtent = Nx * Ny * Nz;
	this->bbox = newBBox;

	this->field = newField;
	this->frozenCells = newFrozenCells;
}

void Grid::clip(const Box3& targetBox)
{
	Vector3 size = targetBox.getSize();
	// resize within bounds
	size.x = std::fminf(size.x, scale.x);
	size.y = std::fminf(size.y, scale.y);
	size.z = std::fminf(size.z, scale.z);
	auto newTargetBox = targetBox;
	newTargetBox.setToSize(&size);

	double cellSize = scale.x / Nx;

	uint Nx_new = std::floor(size.x / cellSize);
	uint Ny_new = std::floor(size.y / cellSize);
	uint Nz_new = std::floor(size.z / cellSize);
	uint gridExtent_new = Nx_new * Ny_new * Nz_new;

	int ix, iy, iz, gridPosOld, gridPosNew;

	auto newField = std::vector<double>(gridExtent_new);
	auto newFrozenCells = std::vector<bool>(gridExtent_new);

	for (iz = 0; iz < Nz_new; iz++) {
		for (iy = 0; iy < Ny_new; iy++) {
			for (ix = 0; ix < Nx_new; ix++) {

				gridPosOld = Nx * Ny * iz + Nx * iy + ix;
				gridPosNew = Nx_new * Ny_new * iz + Nx_new * iy + ix;

				newField[gridPosNew] = this->field[gridPosOld];
				newFrozenCells[gridPosNew] = this->frozenCells[gridPosOld];
			}
		}
	}

	clearField();
	clearFrozenCells();

	this->Nx = Nx_new; this->Ny = Ny_new; this->Nz = Nz_new;
	this->gridExtent = Nx * Ny * Nz;
	this->bbox = Box3(newTargetBox);

	this->field = newField;
	this->frozenCells = newFrozenCells;

	this->scale = size;
}

bool Grid::hasGradient()
{
	return (!gradFieldX.empty() && !gradFieldY.empty() && !gradFieldZ.empty() && gradExtent > 0);
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
	double dx = scale.x / nx;
	double dy = scale.y / ny;
	double dz = scale.z / nz;
	double result_distSq, distSq;

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
					Vector3** t = new Vector3 * [3];
					t[0] = &geom->uniqueVertices[geom->vertexIndices[i]];
					t[1] = &geom->uniqueVertices[geom->vertexIndices[i + 1]];
					t[2] = &geom->uniqueVertices[geom->vertexIndices[i + 2]];
					distSq = getDistanceToATriangleSq(t, &p);
					delete[] t;

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
	double dx = scale.x / nx;
	double dy = scale.y / ny;
	double dz = scale.z / nz;
	double distSq;

	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {

				p.set(
					o.x + ix * dx,
					o.y + iy * dy,
					o.z + iz * dz
				);

				i = aabb->getClosestPrimitiveIdAndDist(p, &distSq);
				if (i < 0) continue;

				gridPos = Nx * Ny * iz + Nx * iy + ix;
				this->field[gridPos] = sqrt(distSq);
			}
		}
	}
}

void Grid::cleanField()
{
	field.clear();
	frozenCells.clear();
	Nx = 0, Ny = 0, Nz = 0;
	scale = Vector3(1.0, 1.0, 1.0);
}

void Grid::cleanGrad()
{
	gradFieldX.clear();
	gradFieldY.clear();
	gradFieldZ.clear();
}

void Grid::scaleBy(Vector3& s)
{
	uint oldNx = Nx, oldNy = Ny, oldNz = Nz;

	this->Nx = (uint)std::floor(s.x * Nx);
	this->Ny = (uint)std::floor(s.y * Ny);
	this->Nz = (uint)std::floor(s.z * Nz);
	auto oldField = field;

	this->gridExtent = this->Nx * this->Ny * this->Nz;
	this->field = std::vector<double>(this->gridExtent);
	Vector3 p = Vector3(), o = bbox.min;
	uint nx = Nx - 1;
	uint ny = Ny - 1;
	uint nz = Nz - 1;
	uint gridPos;
	double dx = this->scale.x / nx;
	double dy = this->scale.y / ny;
	double dz = this->scale.z / nz;
	std::vector<Vector3> positionBuffer = {};
	std::vector<double> valueBuffer = {};

	uint ix, iy, iz;
	
	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {

				p.set(
					o.x + ix * dx,
					o.y + iy * dy,
					o.z + iz * dz
				);				

				positionBuffer.clear(); // for old min and max positions
				valueBuffer.clear(); // for old cell vertex values
				
				this->getSurroundingCells(p, oldNx, oldNy, oldNz, oldField, &positionBuffer, &valueBuffer);

				gridPos = Nx * Ny * iz + Nx * iy + ix;
				field[gridPos] = trilinearInterpolate(p, positionBuffer, valueBuffer);
			}
		}
	}

}

void Grid::getSurroundingCells(Vector3& pos,
	uint oldNx, uint oldNy, uint oldNz, std::vector<double>& oldField,
	std::vector<Vector3>* positionBuffer, std::vector<double>* valueBuffer)
{
	uint ix = std::min((uint)std::floor((pos.x - bbox.min.x) * oldNx / this->scale.x), oldNx - 1);
	uint iy = std::min((uint)std::floor((pos.y - bbox.min.y) * oldNy / this->scale.y), oldNy - 1);
	uint iz = std::min((uint)std::floor((pos.z - bbox.min.z) * oldNz / this->scale.z), oldNz - 1);

	uint ix1 = std::min(ix + 1, oldNx - 1);
	uint iy1 = std::min(iy + 1, oldNy - 1);
	uint iz1 = std::min(iz + 1, oldNz - 1);

	uint i000 = oldNx * oldNy * iz + oldNx * iy + ix;
	uint i100 = oldNx * oldNy * iz + oldNx * iy + ix1;
	uint i010 = oldNx * oldNy * iz + oldNx * iy1 + ix;
	uint i110 = oldNx * oldNy * iz + oldNx * iy1 + ix1;
	uint i001 = oldNx * oldNy * iz1 + oldNx * iy + ix;
	uint i101 = oldNx * oldNy * iz1 + oldNx * iy + ix1;
	uint i011 = oldNx * oldNy * iz1 + oldNx * iy1 + ix;
	uint i111 = oldNx * oldNy * iz1 + oldNx * iy1 + ix1;

	valueBuffer->reserve(8);
	valueBuffer->push_back(oldField[i000]);
	valueBuffer->push_back(oldField[i100]);
	valueBuffer->push_back(oldField[i010]);
	valueBuffer->push_back(oldField[i110]);
	valueBuffer->push_back(oldField[i001]);
	valueBuffer->push_back(oldField[i101]);
	valueBuffer->push_back(oldField[i011]);
	valueBuffer->push_back(oldField[i111]);

	positionBuffer->reserve(2);
	positionBuffer->push_back(Vector3(
		bbox.min.x + (double)ix / (double)oldNx * this->scale.x,
		bbox.min.y + (double)iy / (double)oldNy * this->scale.y,
		bbox.min.z + (double)iz / (double)oldNz * this->scale.z
	));

	positionBuffer->push_back(Vector3(
		bbox.min.x + (double)(ix + 1) / (double)oldNx * this->scale.x,
		bbox.min.y + (double)(iy + 1) / (double)oldNy * this->scale.y,
		bbox.min.z + (double)(iz + 1) / (double)oldNz * this->scale.z
	));
}

double Grid::getL2Norm()
{
	double result = 0.0;

	uint gridPos, ix, iy, iz;

	for (iz = 0; iz < Nz; iz++) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {

				gridPos = Nx * Ny * iz + Nx * iy + ix;
				result += this->field[gridPos] * this->field[gridPos];
			}
		}
	}

	return sqrt(result / (Nx * Ny * Nz));
}

void Grid::clearFrozenCells()
{
	frozenCells.clear();
}

void Grid::clearField()
{
	field.clear();
}

Vector3 Grid::grad(Vector3& p, Vector3& dXYZ, std::vector<Vector3>* positionBuffer, std::vector<double>* valueBuffer)
{
	Vector3 gradf_p = Vector3();
	
	uint ix = std::min((uint)std::floor((p.x - bbox.min.x) * Nx / this->scale.x), Nx - 1);
	uint iy = std::min((uint)std::floor((p.y - bbox.min.y) * Ny / this->scale.y), Ny - 1);
	uint iz = std::min((uint)std::floor((p.z - bbox.min.z) * Nz / this->scale.z), Nz - 1);

	uint ix1 = std::min(ix + 1, Nx - 1);
	uint iy1 = std::min(iy + 1, Ny - 1);
	uint iz1 = std::min(iz + 1, Nz - 1);

	uint i000 = Nx * Ny * iz + Nx * iy + ix;
	uint i100 = Nx * Ny * iz + Nx * iy + ix1;
	uint i010 = Nx * Ny * iz + Nx * iy1 + ix;
	uint i110 = Nx * Ny * iz + Nx * iy1 + ix1;
	uint i001 = Nx * Ny * iz1 + Nx * iy + ix;
	uint i101 = Nx * Ny * iz1 + Nx * iy + ix1;
	uint i011 = Nx * Ny * iz1 + Nx * iy1 + ix;
	uint i111 = Nx * Ny * iz1 + Nx * iy1 + ix1;

	valueBuffer->reserve(8);
	valueBuffer->push_back(field[i000]);
	valueBuffer->push_back(field[i100]);
	valueBuffer->push_back(field[i010]);
	valueBuffer->push_back(field[i110]);
	valueBuffer->push_back(field[i001]);
	valueBuffer->push_back(field[i101]);
	valueBuffer->push_back(field[i011]);
	valueBuffer->push_back(field[i111]);

	return gradf_p;
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

double trilinearInterpolate(Vector3& P, std::vector<Vector3>& X, std::vector<double>& f)
{
	double x = P.x, y = P.y, z = P.z;
	double x0 = X[0].x, y0 = X[0].y, z0 = X[0].z; // cell min
	double x1 = X[1].x, y1 = X[1].y, z1 = X[1].z; // cell max

	// cell values
	double c000 = f[0], c100 = f[1], c010 = f[2], c110 = f[3];
	double c001 = f[4], c101 = f[5], c011 = f[6], c111 = f[7];

	double det = (x0 - x1) * (y0 - y1) * (z0 - z1);

	double a0 =
		(c111 * x0 * y0 * z0 - c011 * x1 * y0 * z0 - c101 * x0 * y1 * z0 + c001 * x1 * y1 * z0 -
		 c110 * x0 * y0 * z1 + c010 * x1 * y0 * z1 + c100 * x0 * y1 * z1 - c000 * x1 * y1 * z1) / det;

	double a1 =
		(c011 * y0 * z0 - c111 * y0 * z0 - c001 * y1 * z0 + c101 * y1 * z0 - c010 * y0 * z1 +
		 c110 * y0 * z1 + c000 * y1 * z1 - c100 * y1 * z1) / det;

	double a2 =
		(c101 * x0 * z0 - c111 * x0 * z0 - c001 * x1 * z0 + c011 * x1 * z0 - c100 * x0 * z1 +
		 c110 * x0 * z1 + c000 * x1 * z1 - c010 * x1 * z1) / det;

	double a3 =
		(c110 * x0 * y0 - c111 * x0 * y0 - c010 * x1 * y0 + c011 * x1 * y0 - c100 * x0 * y1 +
		 c101 * x0 * y1 + c000 * x1 * y1 - c001 * x1 * y1) / det;

	double a4 =
		(c001 * z0 - c011 * z0 - c101 * z0 + c111 * z0 - c000 * z1 + c010 * z1 + c100 * z1 -
		 c110 * z1) / det;

	double a5 =
		(c010 * y0 - c011 * y0 - c110 * y0 + c111 * y0 - c000 * y1 + c001 * y1 + c100 * y1 -
		 c101 * y1) / det;

	double a6 =
		(c100 * x0 - c101 * x0 - c110 * x0 + c111 * x0 - c000 * x1 + c001 * x1 + c010 * x1 -
		 c011 * x1) / det;

	double a7 =
		(c000 - c001 - c010 + c011 - c100 + c101 + c110 - c111) / det;


	return a0 + a1 * x + a2 * y + a3 * z + a4 * x * y + a5 * x * z + a6 * y * z + a7 * x * y * z;
}
