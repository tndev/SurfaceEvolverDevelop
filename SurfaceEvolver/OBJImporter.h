#ifndef OBJIMPORTER_H_
#define OBJIMPORTER_H_

#include <iostream>
#include <string>
#include <fstream>
#include "Geometry.h"
#include "Vector3.h"

class OBJImporter
{
public:
	std::string pathPrefix = ".\\";
	OBJImporter();
	~OBJImporter();

	Geometry importOBJGeometry(std::string filename);
};

#endif
