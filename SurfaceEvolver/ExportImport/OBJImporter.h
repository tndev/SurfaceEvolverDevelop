#ifndef OBJIMPORTER_H_
#define OBJIMPORTER_H_

#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include "../Geometry/Geometry.h"
#include "../Geometry/Vector3.h"

class OBJImporter
{
public:
	std::string pathPrefix = ".\\";
	OBJImporter();
	~OBJImporter();

	Geometry importOBJGeometry(std::string filename);
private:
	void setGeometry(
		Geometry& geom, std::vector<Vector3>& vertices, std::vector<Vector3>& normals,
		std::vector<unsigned int>& vertexIndices, std::vector<unsigned int>& normalIndices,
		unsigned int nTriangles, std::vector<unsigned int>& faceVertexCounts);
};

#endif
