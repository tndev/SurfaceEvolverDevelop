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
};

#endif