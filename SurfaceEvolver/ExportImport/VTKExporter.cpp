#include "VTKExporter.h"

VTKExporter::VTKExporter(std::string outputType)
{
	this->outputType = outputType;
}

VTKExporter::~VTKExporter()
{
}

void VTKExporter::initExport(Geometry object, std::string filename)
{
	if (object.triangulations.size() > 0) {
		// analyze input geometry
		std::pair<std::vector<BufferGeom::Triangulation>, std::vector<size_t>> triangulationsAndSizes = object.getSortedPolygonTriangulationsAndSizes();
		unsigned int NPolyTypes = triangulationsAndSizes.second.size();
		this->outputType = NPolyTypes == 1 ? "POLYDATA" : "UNSTRUCTURED_GRID";
	}
	else {
		std::cout << filename << ".vtk: attempting to write non-triangulated geometry" << std::endl;
	}

	// std::string suffix = this->outputType == "POLYDATA" ? ".vtk" : ".vtk";
	std::string suffix = ".vtk";

	std::fstream vtk(pathPrefix + filename + suffix, std::fstream::out);

	std::vector<Vector3> uniqueVertices = object.uniqueVertices;
	size_t pointCount = uniqueVertices.size();

	vtk << "# vtk DataFile Version 4.2" << std::endl;
	vtk << "vtk output" << std::endl;
	vtk << "ASCII" << std::endl;
	vtk << "DATASET " << this->outputType << std::endl;
	vtk << "POINTS " << pointCount << " float" << std::endl;

	// std::string newline = this->outputType == "POLYDATA" ? "\n" : "\n";
	std::string newline = "\n";

	if (pointCount > 0) {
		for (int i = 0; i < pointCount; i++) {
			vtk << uniqueVertices[i].x << " " << uniqueVertices[i].y << " " << uniqueVertices[i].z << newline;
		}
	}

	vtk << std::endl << std::endl;

	auto writePolygonIndices = [&](std::string cellType) {
		size_t NTriIds = this->countTriangulationIndices(object.triangulations);
		size_t NTriangulations = object.triangulations.size();

		std::vector<unsigned int> polySizes = {};

		for (unsigned int i = 0; i < NTriangulations; i++) {
			unsigned int tSize = object.triangulations[i].size();
			unsigned int vtkRowLength = 3 + tSize;

			if (i == 0 && this->outputType == "POLYDATA") {
				vtk << cellType << " " << NTriangulations << " " << vtkRowLength * NTriangulations << " " << std::endl;
			}
			else if (i == 0 && this->outputType == "UNSTRUCTURED_GRID") {
				vtk << cellType << " " << NTriangulations << " " << NTriIds + NTriangulations << " " << std::endl;
			}

			std::vector<unsigned int> polygonIds = object.getPolygonIndicesFromTriangulation(object.triangulations[i]);
			polySizes.push_back(vtkRowLength - 1);
			vtk << vtkRowLength - 1 << " ";
			for (unsigned int k = 0; k < polygonIds.size(); k++) {
				vtk << polygonIds[k] << (k < polygonIds.size() - 1 ? " " : "\n");
			}
		}

		if (this->outputType == "UNSTRUCTURED_GRID") {
			vtk << std::endl;
			vtk << "CELL_TYPES " << NTriangulations << std::endl;

			for (unsigned int i = 0; i < NTriangulations; i++) {
				unsigned int polyType;
				unsigned int ps = polySizes[i];
				if (ps == 3) {
					polyType = 5; // VTK_TRIANGLE
				}
				else if (ps == 4) {
					polyType = 9; // VTK_QUAD
				}
				else {
					polyType = 7; // VTK_POLYGON
				}
				vtk << polyType << std::endl;
			}
		}
	};

	std::string cellType = (this->outputType == "POLYDATA" ? "POLYGONS" : (this->outputType == "UNSTRUCTURED_GRID" ? "CELLS" : ""));

	if (object.hasTriangulations()) {
		writePolygonIndices(cellType);
	}

	vtk.close();
}

void VTKExporter::exportPointData(std::vector<Vector3> points, std::string filename)
{
	std::string suffix = ".vtk";

	std::fstream vtk(pathPrefix + filename + suffix, std::fstream::out);

	vtk << "# vtk DataFile Version 4.2" << std::endl;
	vtk << "vtk output" << std::endl;
	vtk << "ASCII" << std::endl;
	vtk << "DATASET " << this->outputType << std::endl;
	vtk << "POINTS " << points.size() << " float" << std::endl << std::endl;

	for (auto&& p : points) {
		vtk << p.x << " " << p.y << " " << p.z << std::endl;
	}

	vtk.close();
}

size_t VTKExporter::countTriangulationIndices(std::vector<BufferGeom::Triangulation>& triangulations)
{
	std::vector<BufferGeom::Triangulation> T = triangulations;
	size_t nCellIds = 0;

	for (std::vector<BufferGeom::Triangulation>::iterator it = T.begin(); it < T.end(); ++it) {
		size_t currentSize = it->size();
		nCellIds += currentSize;
	}

	return nCellIds;
}
