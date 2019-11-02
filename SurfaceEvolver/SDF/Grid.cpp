#include "Grid.h"

Grid::Grid()
{
}

Grid::Grid(uint Nx, uint Ny, uint Nz, Box3 bbox)
{
	this->Nx = Nx; this->Ny = Ny; this->Nz = Nz;
	this->bbox = bbox;
	this->scale = this->bbox.getSize();
	this->field = std::vector<float>((size_t)this->Nx * this->Ny * this->Nz); // init field
}

Grid::~Grid()
{
}

void Grid::exportToRawBinary(std::string filename)
{
	size_t fieldSize = (size_t)Nx * Ny * Nz;
	size_t bytes = sizeof(char) * fieldSize;
	char* data = new char[fieldSize];
	for (uint i = 0; i < fieldSize; i++) {
		data[i] = (char)field[i];
	}

	auto raw = std::fstream(filename + ".raw", std::ios::out | std::ios::binary);
	raw.write((char*)&data[0], bytes);
	raw.close();
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
