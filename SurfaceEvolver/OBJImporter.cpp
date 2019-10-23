#include "OBJImporter.h"

OBJImporter::OBJImporter()
{
}

OBJImporter::~OBJImporter()
{
}

Geometry OBJImporter::importOBJGeometry(std::string filename)
{
	Geometry result = Geometry();

	if (filename.compare(".obj") != 0) {
		std::cout << "File " << filename << " has invalid suffix!" << std::endl;
		return result;
	}

	std::fstream obj(pathPrefix + filename, std::fstream::in);

	if (!obj.is_open()) {
		std::cout << "File " << filename << " was not found!" << std::endl;
		return result;
	}

	std::cout << filename << " opened successfully" << std::endl;

	std::string line;

	auto getTokens = [&](std::string str, std::vector<std::string>& tokens, std::string delimiter, size_t buffSize) {
		size_t pos = 0;
		while (pos <= buffSize) {
			pos = str.find(delimiter);
			tokens.push_back(str.substr(0, str.find(delimiter)));
			str = str.erase(0, pos + delimiter.length());
		}
	};

	while (std::getline(obj, line)) {
		if (line.compare(0, 2, "v ") == 0) {
			std::vector<std::string> tokens = std::vector<std::string>();
			getTokens(line, tokens, " ", line.size());

			float vx = std::stof(tokens[1]);
			float vy = std::stof(tokens[2]);
			float vz = std::stof(tokens[3]);

			Vector3 vertex = Vector3(vx, vy, vz);
			result.uniqueVertices.push_back(vertex);
		}
		else if (line.compare(0, 3, "vt ") == 0) {
			// TODO : tangents
		}
		else if (line.compare(0, 3, "vn ") == 0) {
			std::vector<std::string> tokens = std::vector<std::string>();
			getTokens(line, tokens, " ", line.size());

			float nx = std::stof(tokens[1]);
			float ny = std::stof(tokens[2]);
			float nz = std::stof(tokens[3]);

			result.normals.push_back(nx);
			result.normals.push_back(ny);
			result.normals.push_back(nz);
		}
		else if (line.compare(0, 2, "f ") == 0) {
			std::vector<std::string> tokens = std::vector<std::string>();
			getTokens(line, tokens, " ", line.size());

			std::vector<unsigned int> polygonIndices = std::vector<unsigned int>();

			for (unsigned int i = 1; i < tokens.size(); i++) {
				std::vector<std::string> subTokens = std::vector<std::string>();
				getTokens(tokens[i], subTokens, "/", tokens[i].size());

				unsigned int vi = std::stoi(subTokens[0]);
				polygonIndices.push_back(vi);

				// TODO: Triangulate vertex indices
				// TODO: what to do with vertex normal indices and tangent indices?
			}
		}
	}

	std::string line;
	int i = 0;

	return result;
}
