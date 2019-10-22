#include "VTKExporter.h"

VTKExporter::VTKExporter()
{
}

VTKExporter::~VTKExporter()
{
}

void VTKExporter::initExport(Geometry object, std::string filename)
{
	size_t pointCount = object.uniqueVertices.size() / 3;
	std::fstream vtk(pathPrefix + filename + ".vtk", std::fstream::out);

	vtk << "# vtk DataFile Version 4.2" << std::endl;
	vtk << "vtk output" << std::endl;
	vtk << "ASCII" << std::endl;
	vtk << "DATASET POLYDATA" << std::endl;
	vtk << "POINTS " << pointCount << " float" << std::endl;

	if (pointCount > 0) {
		for (int i = 0; i < pointCount; i++) {
			unsigned int i0 = 3 * i, i1 = 3 * i + 1, i2 = 3 * i + 2;
			vtk << object.uniqueVertices[i0] << " " << object.uniqueVertices[i1] << " " << object.uniqueVertices[i2] << std::endl;
		}
	}

	// TODO: divide all polygons into groups according to their number of sides and write those as separate index groups in the file
	if (object.hasTriangulations()) {
		size_t polyCount = object.triangulations.size();
		unsigned int vtkRowLength = 3 + object.triangulations[0].size(); // assuming all faces have the same number of sides

		vtk << "POLYGONS" << " " << polyCount << " " << vtkRowLength * polyCount << " " << std::endl;

		for (unsigned int i = 0; i < polyCount; i++) {
			std::vector<unsigned int> t = object.triangulations[i];
			std::vector<std::vector<unsigned int>> triangles = std::vector<std::vector<unsigned int>>();
			for (unsigned int j = 0; j < t.size(); j++) {
				triangles.push_back({object.vertexIndices[3 * t[j]], object.vertexIndices[3 * t[j] + 1], object.vertexIndices[3 * t[j] + 2]});
			}
			std::vector<unsigned int> polygonIds = object.getPolygonVerticesFromTriangulation(triangles);

			vtk << vtkRowLength - 1 << " ";
			for (unsigned int j = 0; j < polygonIds.size(); j++) {
				vtk << polygonIds[j] << (j < polygonIds.size() - 1 ? " " : "\n");
			}
		}
	}

	vtk.close();
}
