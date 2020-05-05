#ifndef VTKEXPORTER_H_
#define VTKEXPORTER_H_

#include <iostream>
#include <fstream>
#include <string>
#include "../Geometry/Geometry.h"


class VTKExporter
{
public:
	std::string outputType = "POLYDATA";
	std::string pathPrefix = ".\\";
	VTKExporter(std::string outputType = "POLYDATA");
	~VTKExporter();

	void initExport(Geometry* object, std::string filename);
	void exportPointData(std::vector<Vector3> points, std::string filename);
	void exportGeometryVertexNormals(Geometry* object, std::string filename);
	// vertId = -1 exports all FV geometries, triId = -1 all triangles per fv
	void exportGeometryFiniteVolumeGrid(
		Geometry* object, std::vector<std::vector<Vector3>>& fvVerts, std::vector<std::vector<std::vector<uint>>>& adjacentPolys, 
		std::string filename, int vertId = -1, int triId = -1, bool fromStartToVertId = true);
	void exportVectorDataOnGeometry(Geometry* object, std::vector<Vector3>* data, std::string filename);
private:
	size_t countTriangulationIndices(std::vector<BufferGeom::Triangulation>& triangulations);
};

#endif