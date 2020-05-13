#include "OBJImporter.h"

bool isNaN(unsigned int x) {
	return x != x;
}

OBJImporter::OBJImporter()
{
}

OBJImporter::~OBJImporter()
{
}

Geometry OBJImporter::importOBJGeometry(std::string filename)
{
	if (filename.compare(".obj") == 0) {
		std::cout << "File " << filename << " has invalid suffix!" << std::endl;
		return Geometry();
	}

	std::fstream obj(pathPrefix + filename, std::fstream::in);

	if (!obj.is_open()) {
		std::cout << "File " << filename << " was not found!" << std::endl;
		return Geometry();
	}

	std::cout << filename << " opened successfully" << std::endl;

	int suffix_pos = filename.find(".obj");
	filename = filename.erase(suffix_pos, filename.size() - 1);
	Geometry result = Geometry(filename);

	std::string line;

	auto getTokens = [&](std::string str, std::vector<std::string>& tokens, std::string delimiter, size_t buffSize) {
		size_t pos = 0;
		std::string query = "";
		while (pos <= buffSize) {
			pos = str.find(delimiter);
			query = str.substr(0, str.find(delimiter));
			if (query != "") {
				tokens.push_back(query);
				str = str.erase(0, pos + delimiter.length());
			}
		}
	};

	std::vector<Vector3> vertices = std::vector<Vector3>();
	std::vector<Vector3> uvVertices = std::vector<Vector3>();
	std::vector<Vector3> normals = std::vector<Vector3>();

	std::vector<unsigned int> _vertexIndices = std::vector<unsigned int>();
	std::vector<unsigned int> _uvIndices = std::vector<unsigned int>();
	std::vector<unsigned int> _normalIndices = std::vector<unsigned int>();
	std::vector<unsigned int> _faceVertexCounts = std::vector<unsigned int>();

	unsigned int nVertices = 0;
	unsigned int nTriangles = 0;
	std::vector<unsigned int> noVertsOfFace = std::vector<unsigned int>();

	// TODO: some exporters use 'o' (blender), some 'g' (3ds max), some both (modo, rhino), and some don't even bother (meshlab)!
	while (std::getline(obj, line)) {
		if (line.compare(0, 2, "v ") == 0) {
			std::vector<std::string> tokens = std::vector<std::string>();
			getTokens(line, tokens, " ", line.size());

			double vx = std::stof(tokens[1]);
			double vy = std::stof(tokens[2]);
			double vz = std::stof(tokens[3]);

			vertices.push_back(Vector3(vx, vy, vz));
		}
		else if (line.compare(0, 3, "vt ") == 0) {
			std::vector<std::string> tokens = std::vector<std::string>();
			getTokens(line, tokens, " ", line.size());

			double vtx = std::stof(tokens[1]);
			double vty = std::stof(tokens[2]);

			uvVertices.push_back(Vector3(vtx, vty, 0));
		}
		else if (line.compare(0, 3, "vn ") == 0) {
			std::vector<std::string> tokens = std::vector<std::string>();
			getTokens(line, tokens, " ", line.size());

			double nx = std::stof(tokens[1]);
			double ny = std::stof(tokens[2]);
			double nz = std::stof(tokens[3]);

			normals.push_back(Vector3(nx, ny, nz));
		}
		else if (line.compare(0, 2, "f ") == 0) {
			std::vector<std::string> tokens = std::vector<std::string>();
			getTokens(line, tokens, " ", line.size());

			nVertices = 0;

			for (unsigned int i = 1; i < tokens.size(); i++) {
				std::vector<std::string> subTokens = std::vector<std::string>();
				getTokens(tokens[i], subTokens, "/", tokens[i].size());

				unsigned int vi, ui, ni;

				if (subTokens.size() > 0) {
					vi = std::stoi(subTokens[0]);
					vi = vi > 0 ? vi - 1 : (unsigned int)vertices.size() + vi;
					_vertexIndices.push_back(vi);
					if (vi < vertices.size()) {
						nVertices++;
					}
					else {
						std::cout << "Error: OBJImporter: Index out of range!" << std::endl;
					}
				}

				if (subTokens.size() > 1) {
					ui = std::stoi(subTokens[1]);
					ui = ui > 0 ? ui - 1 : (unsigned int)uvVertices.size() + ui;
					_uvIndices.push_back(ui);
				}

				if (subTokens.size() > 2) {
					ni = std::stoi(subTokens[2]);
					ni = ni > 0 ? ni - 1 : (unsigned int)normals.size() + ni;
					_normalIndices.push_back(ni);
				}
			}
			_faceVertexCounts.push_back(nVertices);

			if (nVertices >= 3) {
				nTriangles += ((nVertices - 2) * 3);
			}
			else {
				std::cout << "Error: OBJImporter: Invalid number of vertices after 'f'. Must be at least 3!" << std::endl;
			}
		}
	}

	obj.close();

	setGeometry(result, vertices, normals, _vertexIndices, _normalIndices, nTriangles, _faceVertexCounts);

	return result;
}

void OBJImporter::setGeometry(Geometry& geom,
	std::vector<Vector3>& vertices, std::vector<Vector3>& normals,
	std::vector<unsigned int>& vertexIndices, std::vector<unsigned int>& normalIndices,
	unsigned int nTriangles, std::vector<unsigned int>& faceVertexCounts)
{
	unsigned int idxV = 0;
	unsigned int idxN = 0;
	unsigned int idxI = 0;
	unsigned int totalVerts = 0;
	unsigned int failedTriangulationCount = 0;

	geom.uniqueVertices = std::vector<Vector3>(vertices.size());
	geom.vertices = std::vector<double>(((size_t)3 * nTriangles));
	geom.vertexIndices = std::vector<unsigned int>(nTriangles);
	geom.normals = std::vector<double>(((size_t)3 * nTriangles));
	geom.triangulations = std::vector<BufferGeom::Triangulation>();

	std::map<Vector3, unsigned int> vertexToIdx = std::map<Vector3, unsigned int>();

	BufferGeom::Face helperFace = BufferGeom::Face();

	for (auto&& vertCount:faceVertexCounts) {
		BufferGeom::Face faceVerts = BufferGeom::Face();
		BufferGeom::Face faceNorms = BufferGeom::Face();

		// Some .obj files specify uvs or normals only for some of the vertices
		bool faceHasNormals = normalIndices.size() > 0;

		for (unsigned int k = totalVerts; k < vertCount + totalVerts; k++) {
			faceVerts.push_back(&vertices[vertexIndices[k]]);

			if (faceHasNormals) {
				faceNorms.push_back(&normals[normalIndices[k]]);
			}
		}

		if (!faceHasNormals) {
			helperFace = faceVerts;
			Vector3 defaultNormal = geom.getNormal(helperFace);
			for (unsigned int k = 0; k < vertCount; k++) {
				faceNorms.push_back(&defaultNormal);
			}
		}

		// Two or more elements in one face can be the same, but I'm not dealing with this case for now

		std::vector<std::vector<unsigned int>> faces = geom.getTriangulatedIndices(faceVerts);

		if (faces.size() > 0) {
			for (unsigned int k = 0; k < vertCount; k++) {
				auto found = vertexToIdx.find(*faceVerts[k]);
				if (found == vertexToIdx.end()) {
					vertexToIdx.insert({ *faceVerts[k], vertexToIdx.size() });
				}
			}

			BufferGeom::Triangulation triang = {};
			for (auto&& f:faces) {
				for (auto&& i:f) {
					Vector3 v = *faceVerts[i];
					geom.vertices[idxV++] = v.x;
					geom.vertices[idxV++] = v.y;
					geom.vertices[idxV++] = v.z;

					unsigned int vId = vertexToIdx.at(v);
					geom.vertexIndices[idxI++] = vId;
					geom.uniqueVertices[vId] = v;

					Vector3 norm = *faceNorms[i];
					geom.normals[idxN++] = norm.x;
					geom.normals[idxN++] = norm.y;
					geom.normals[idxN++] = norm.z;

				}
				triang.push_back(idxV / 9 - 1);
			}
			geom.triangulations.push_back(triang);

		}
		else {
			// triangulation failed. The face is most probably messed up (huge, non-planar, non-convex, overlapping polygon)
			failedTriangulationCount += (faceVerts.size() - 2) * 3;
		}

		totalVerts += vertCount;
	}

	if (failedTriangulationCount > 0) {
		// triangulation failed and we skipped some triangles. We need to splice the unused part of the buffers..
		geom.vertices = std::vector<double>(geom.vertices.begin(), geom.vertices.end() - ((size_t)failedTriangulationCount * 3));
		geom.vertexIndices = std::vector<unsigned int>(geom.vertexIndices.begin(), geom.vertexIndices.end() - ((size_t)failedTriangulationCount));
		geom.normals = std::vector<double>(geom.normals.begin(), geom.normals.end() - ((size_t)failedTriangulationCount * 3));
	}
}
