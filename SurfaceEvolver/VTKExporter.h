#ifndef VTKEXPORTER_H_
#define VTKEXPORTER_H_

#include <iostream>
#include <fstream>
#include <string>
#include "Geometry.h"


class VTKExporter
{
public:
	std::string outputType = "POLYDATA";
	std::string pathPrefix = ".\\";
	VTKExporter(std::string outputType = "POLYDATA");
	~VTKExporter();

	void initExport(Geometry object, std::string filename);
private:
	size_t countTriangulationIndices(std::vector<BufferGeom::Triangulation>& triangulations);
};

#endif