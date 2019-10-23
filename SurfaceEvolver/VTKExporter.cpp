#include "VTKExporter.h"

VTKExporter::VTKExporter()
{
}

VTKExporter::~VTKExporter()
{
}

void VTKExporter::initExport(Geometry object, std::string filename)
{
	std::fstream vtk(pathPrefix + filename + ".vtk", std::fstream::out);

	std::vector<Vector3> uniqueVertices = object.uniqueVertices;
	size_t pointCount = uniqueVertices.size();

	vtk << "# vtk DataFile Version 4.2" << std::endl;
	vtk << "vtk output" << std::endl;
	vtk << "ASCII" << std::endl;
	vtk << "DATASET POLYDATA" << std::endl;
	vtk << "POINTS " << pointCount << " float" << std::endl;

	if (pointCount > 0) {
		for (int i = 0; i < pointCount; i++) {
			vtk << uniqueVertices[i].x << " " << uniqueVertices[i].y << " " << uniqueVertices[i].z << std::endl;
		}
	}

	// divide all polygons into groups according to their number of sides and write those as separate index groups in the file
	if (object.hasTriangulations()) {
		auto triAndSizes = getSortedPolygonTriangulationsAndSizes(object.triangulations);
		std::vector<Triangulation> polygonTriangulations = triAndSizes.first;
		std::vector<size_t> sizes = triAndSizes.second;
		size_t NPolyTypes = sizes.size();
		unsigned int pos = 0;

		vtk << std::endl;

		for (unsigned int i = 0; i < NPolyTypes; i++) {
			size_t polyCount = sizes[i];
			unsigned int tSize = polygonTriangulations[pos].size();
			unsigned int vtkRowLength = 3 + tSize;

			vtk << "POLYGONS" << " " << polyCount << " " << vtkRowLength * polyCount << " " << std::endl;
			
			for (unsigned int j = 0; j < polyCount; j++) {
				unsigned int ptj = pos + j;
				std::vector<unsigned int> polygonIds = object.getPolygonIndicesFromTriangulation(polygonTriangulations[ptj]);
				vtk << vtkRowLength - 1 << " ";
				for (unsigned int k = 0; k < polygonIds.size(); k++) {
					vtk << polygonIds[k] << (k < polygonIds.size() - 1 ? " " : "\n");
				}
			}

			pos += polyCount;
		}
	}

	vtk.close();
}

std::pair<std::vector<Triangulation>, std::vector<size_t>> VTKExporter::getSortedPolygonTriangulationsAndSizes(std::vector<Triangulation>& triangulations)
{
	std::vector<Triangulation> T = triangulations;
	std::vector<size_t> sizes = std::vector<size_t>();
	std::sort(T.begin(), T.end(), [](const Triangulation& a, const Triangulation& b) { return a.size() < b.size(); });
	unsigned int tCount = 0; // triangulation count
	size_t currentSize = T.begin()->size();

	for (std::vector<Triangulation>::iterator it = T.begin(); it < T.end(); ++it) {
		if (it->size() > currentSize) {
			currentSize = it->size();
			sizes.push_back(tCount);
			tCount = 0;
			continue;
		}
		tCount++;
	}
	sizes.push_back(tCount);

	std::pair<std::vector<Triangulation>, std::vector<size_t>> result = { T, sizes };
	return result;
}
