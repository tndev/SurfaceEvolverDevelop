#include "Grid.h"

Grid::Grid()
{
}

Grid::Grid(const Grid& other)
{
	this->field = std::vector<float>(other.field);
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

	for (uint i = 0; i < field.size(); i += 9) {
		vti << "					";
		for (uint j = 0; j < 9; j++) {
			vti << this->field[(size_t)i + j] << " ";
		}
		vti << std::endl;
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

void Grid::clean()
{
	Nx = 0, Ny = 0, Nz = 0;
	scale = Vector3(1.0f, 1.0f, 1.0f);
	field.clear();
}
