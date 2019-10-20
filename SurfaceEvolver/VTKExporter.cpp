#include "VTKExporter.h"

VTKExporter::VTKExporter()
{
}

VTKExporter::~VTKExporter()
{
}

void VTKExporter::initExport(Geometry object, std::string filename)
{
	unsigned int pointCount = object.vertices.size() / 3;
	std::fstream vtk(pathPrefix + filename, std::fstream::out);

	vtk << "# vtk DataFile Version 4.1" << std::endl;
	vtk << "vtk output" << std::endl;
	vtk << "ASCII" << std::endl;
	vtk << "DATASET POLYDATA" << std::endl;
	vtk << "POINTS " << pointCount << " float" << std::endl;

	if (pointCount > 0) {
		for (int i = 0; i < pointCount; i++) {
			unsigned int i0 = 3 * i, i1 = 3 * i + 1, i2 = 3 * i + 2;
			vtk << object.vertices[i0] << " " << object.vertices[i1] << " " << object.vertices[i2] << std::endl;
		}
	}

	// Temporary solution: hasTriangulations corresponds to quads;
	bool hasTriangulations = object.hasTriangulations();
	unsigned int step = 3 + (hasTriangulations ? 3 : 0);
	unsigned int polyCount = object.vertexIndices.size() / step;
	unsigned int vtkRowLength = 4 + (hasTriangulations ? 1 : 0);

	vtk << "POLYGONS" << " " << polyCount << " " << vtkRowLength * polyCount << " " << std::endl;

	if (polyCount > 0 && hasTriangulations) {
		unsigned int NPoly = object.triangulations.size();
		for (unsigned int i = 0; i < polyCount; i++)	{
			unsigned int i0 = 3 * i, i1 = 3 * i + 1, i2 = 3 * i + 2, i3 = 3 * (i + 1) + 2;
			vtk << 4 << " " << object.vertexIndices[i0] << " " << object.vertexIndices[i1] << " " << object.vertexIndices[i2] << " " << object.vertexIndices[i3] << std::endl;
		}
	} else if (polyCount > 0 && !hasTriangulations) {
		for (int i = 0; i < polyCount; i++) {
			unsigned int i0 = 3 * i, i1 = 3 * i + 1, i2 = 3 * i + 2;
			vtk << 4 << " " << object.vertexIndices[i0] << " " << object.vertexIndices[i1] << " " << object.vertexIndices[i2] << std::endl;
		}
	}

	vtk.close();
}
