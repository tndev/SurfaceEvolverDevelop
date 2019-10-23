#ifndef VTKEXPORTER_H_
#define VTKEXPORTER_H_

#include <iostream>
#include <fstream>
#include <string>
#include "Geometry.h"


class VTKExporter
{
public:
	std::string pathPrefix = ".\\";
	VTKExporter();
	~VTKExporter();

	void initExport(Geometry object, std::string filename);
private:
	std::pair<std::vector<Triangulation>, std::vector<size_t>> getSortedPolygonTriangulationsAndSizes(std::vector<Triangulation>& triangulations);
};

#endif