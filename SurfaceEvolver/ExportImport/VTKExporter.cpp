#include "VTKExporter.h"

VTKExporter::VTKExporter(std::string outputType)
{
	this->outputType = outputType;
}

VTKExporter::~VTKExporter()
{
}

void VTKExporter::initExport(Geometry* object, std::string filename)
{
	if (object->triangulations.size() > 0) {
		// analyze input geometry
		std::pair<std::vector<BufferGeom::Triangulation>, std::vector<size_t>> triangulationsAndSizes = object->getSortedPolygonTriangulationsAndSizes();
		unsigned int NPolyTypes = triangulationsAndSizes.second.size();
		this->outputType = NPolyTypes == 1 ? "POLYDATA" : "UNSTRUCTURED_GRID";
	}
	else {
		std::cout << filename << ".vtk: attempting to write non-triangulated geometry" << std::endl;
	}

	// std::string suffix = this->outputType == "POLYDATA" ? ".vtk" : ".vtk";
	std::string suffix = ".vtk";

	std::fstream vtk(pathPrefix + filename + suffix, std::fstream::out);

	std::vector<Vector3>* uniqueVertices = &object->uniqueVertices;
	size_t pointCount = uniqueVertices->size();

	vtk << "# vtk DataFile Version 4.2" << std::endl;
	vtk << "vtk output" << std::endl;
	vtk << "ASCII" << std::endl;
	vtk << "DATASET " << this->outputType << std::endl;
	vtk << "POINTS " << pointCount << " float" << std::endl;

	// std::string newline = this->outputType == "POLYDATA" ? "\n" : "\n";
	std::string newline = "\n";

	if (pointCount > 0) {
		for (int i = 0; i < pointCount; i++) {
			vtk << uniqueVertices->at(i).x << " " << uniqueVertices->at(i).y << " " << uniqueVertices->at(i).z << newline;
		}
	}

	vtk << std::endl << std::endl;

	auto writePolygonIndices = [&](std::string cellType) {
		size_t NTriIds = this->countTriangulationIndices(object->triangulations);
		size_t NTriangulations = object->triangulations.size();

		std::vector<unsigned int> polySizes = {};

		for (unsigned int i = 0; i < NTriangulations; i++) {
			unsigned int tSize = object->triangulations[i].size();
			unsigned int vtkRowLength = 3 + tSize;

			if (i == 0 && this->outputType == "POLYDATA") {
				vtk << cellType << " " << NTriangulations << " " << vtkRowLength * NTriangulations << " " << std::endl;
			}
			else if (i == 0 && this->outputType == "UNSTRUCTURED_GRID") {
				vtk << cellType << " " << NTriangulations << " " << NTriIds + NTriangulations << " " << std::endl;
			}

			std::vector<unsigned int> polygonIds = object->getPolygonIndicesFromTriangulation(object->triangulations[i]);
			polySizes.push_back(vtkRowLength - 1);
			vtk << vtkRowLength - 1 << " ";
			for (unsigned int k = 0; k < polygonIds.size(); k++) {
				vtk << polygonIds[k] << (k < polygonIds.size() - 1 ? " " : "\n");
			}
		}

		if (this->outputType == "UNSTRUCTURED_GRID") {
			vtk << std::endl;
			vtk << "CELL_TYPES " << NTriangulations << std::endl;

			for (unsigned int i = 0; i < NTriangulations; i++) {
				unsigned int polyType;
				unsigned int ps = polySizes[i];
				if (ps == 3) {
					polyType = 5; // VTK_TRIANGLE
				}
				else if (ps == 4) {
					polyType = 9; // VTK_QUAD
				}
				else {
					polyType = 7; // VTK_POLYGON
				}
				vtk << polyType << std::endl;
			}
		}
	};

	auto writeScalarData = [&]() {
		size_t NTables = object->scalarTables.size();
		vtk << std::endl;
		vtk << "POINT_DATA " << pointCount << std::endl;

		for (unsigned int t = 0; t < NTables; t++) {
			vtk << "SCALARS " << object->scalarTables[t].name << " float 1" << std::endl;
			vtk << "LOOKUP_TABLE default" << std::endl;

			for (unsigned int i = 0; i < pointCount; i++) {
				vtk << object->scalarTables[t][i] << std::endl;
			}
		}
	};

	std::string cellType = (this->outputType == "POLYDATA" ? "POLYGONS" : (this->outputType == "UNSTRUCTURED_GRID" ? "CELLS" : ""));

	if (object->hasTriangulations()) {
		writePolygonIndices(cellType);
	}

	if (object->hasScalarData()) {
		writeScalarData();
	}

	vtk.close();
}

void VTKExporter::exportPointData(std::vector<Vector3> points, std::string filename)
{
	std::string suffix = ".vtk";

	std::fstream vtk(pathPrefix + filename + suffix, std::fstream::out);

	vtk << "# vtk DataFile Version 4.2" << std::endl;
	vtk << "vtk output" << std::endl;
	vtk << "ASCII" << std::endl;
	vtk << "DATASET " << this->outputType << std::endl;
	vtk << "POINTS " << points.size() << " float" << std::endl << std::endl;

	for (auto&& p : points) {
		vtk << p.x << " " << p.y << " " << p.z << std::endl;
	}

	vtk.close();
}

void VTKExporter::exportGeometryVertexNormals(Geometry* object, std::string filename)
{
	std::string suffix = ".vtk";

	std::fstream vtk(pathPrefix + filename + suffix, std::fstream::out);

	std::vector<Vector3> uniqueVertices = object->uniqueVertices;
	size_t pointCount = uniqueVertices.size();

	vtk << "# vtk DataFile Version 4.2" << std::endl;
	vtk << "vtk output" << std::endl;
	vtk << "ASCII" << std::endl;
	vtk << "DATASET UNSTRUCTURED_GRID" << std::endl;
	vtk << "POINTS " << pointCount << " float" << std::endl;

	if (pointCount > 0) {
		for (int i = 0; i < pointCount; i++) {
			vtk << uniqueVertices[i].x << " " << uniqueVertices[i].y << " " << uniqueVertices[i].z << std::endl;
		}
	}

	vtk << "POINT_DATA " << pointCount << std::endl;
	vtk << "VECTORS normals float" << std::endl;
	if (pointCount > 0) {
		std::vector<Vector3> pseudoNormals = object->getAngleWeightedVertexPseudoNormals();

		for (auto&& n : pseudoNormals) {
			vtk << n.x << " " << n.y << " " << n.z << std::endl;
		}
	}

	vtk.close();
}

void VTKExporter::exportGeometryFiniteVolumeGrid(
	Geometry* object, 
	std::vector<std::vector<Vector3>>& fvVerts, std::vector<std::vector<std::vector<uint>>>& adjacentPolys, 
	std::string filename, int vertId, int triId, bool fromStartToVertId)
{
	std::vector<Vector3>* uniqueVertices = &object->uniqueVertices;
	size_t pointCount = uniqueVertices->size();
	std::vector<Geometry> fvGeometries = {};

	uint i_min = ((vertId == -1 || fromStartToVertId) ? 0 : std::min(vertId, (int)fvVerts.size()));
	uint i_max = (vertId == -1 ? fvVerts.size() : std::min(vertId + 1, (int)fvVerts.size()));
	vertId = (fromStartToVertId && i_max == fvVerts.size() ? i_max : vertId);

	for (uint i = i_min; i < i_max; i++) {
		Geometry fvGeom = Geometry();

		fvGeom.uniqueVertices = { object->uniqueVertices[i] };
		int oddRingCount = fvVerts[i].size() % 2;
		if (triId == -1 || triId == fvVerts[i].size() - 1) {
			for (uint j = 0; j < fvVerts[i].size(); j += 2) {
				fvGeom.uniqueVertices.push_back(fvVerts[i][j]);
				if (oddRingCount && j == fvVerts[i].size() - 1) {
					continue;
				}
				fvGeom.uniqueVertices.push_back(fvVerts[i][j + 1]);
				
				fvGeom.vertexIndices.push_back(0);
				fvGeom.vertexIndices.push_back(j + 1);
				fvGeom.vertexIndices.push_back(j + 2);
				
				fvGeom.vertexIndices.push_back(0);
				fvGeom.vertexIndices.push_back(j + 2);
				uint modLastId = (oddRingCount && j == fvVerts[i].size() - 2 ? (j + 3) % fvVerts[i].size() : j + 3);
				fvGeom.vertexIndices.push_back(modLastId);

				if (triId == -1) {
					fvGeom.triangulations.push_back({ j, j + 1 });
				}
				else {
					fvGeom.triangulations.push_back({ j });
					fvGeom.triangulations.push_back({ j + 1 });
				}
			}
		}
		else {
			fvGeom.uniqueVertices.push_back(fvVerts[i][0]);
			for (uint j = 0; j <= triId; j++) {
				fvGeom.uniqueVertices.push_back(fvVerts[i][j + 1]);

				fvGeom.vertexIndices.push_back(0);
				fvGeom.vertexIndices.push_back(j + 1);
				fvGeom.vertexIndices.push_back(j + 2);

				fvGeom.triangulations.push_back({ j });
			}
		}

		fvGeom.vertices = std::vector<float>(3 * fvGeom.vertexIndices.size());
		for (uint j = 0; j < fvGeom.vertexIndices.size(); j++) {
			fvGeom.vertices[3 * j] = fvGeom.uniqueVertices[fvGeom.vertexIndices[j]].x;
			fvGeom.vertices[3 * j + 1] = fvGeom.uniqueVertices[fvGeom.vertexIndices[j]].y;
			fvGeom.vertices[3 * j + 2] = fvGeom.uniqueVertices[fvGeom.vertexIndices[j]].z;
		}
		fvGeometries.push_back(fvGeom);
	}

	Geometry result = mergeGeometries(fvGeometries);
	this->initExport(&result, filename);
}

void VTKExporter::exportVectorDataOnGeometry(Geometry* object, std::vector<Vector3>* data, std::string filename)
{
	std::string suffix = ".vtk";

	std::fstream vtk(pathPrefix + filename + suffix, std::fstream::out);

	std::vector<Vector3> uniqueVertices = object->uniqueVertices;
	size_t pointCount = uniqueVertices.size();

	vtk << "# vtk DataFile Version 4.2" << std::endl;
	vtk << "vtk output" << std::endl;
	vtk << "ASCII" << std::endl;
	vtk << "DATASET UNSTRUCTURED_GRID" << std::endl;
	vtk << "POINTS " << pointCount << " float" << std::endl;

	if (pointCount > 0) {
		for (int i = 0; i < pointCount; i++) {
			vtk << uniqueVertices[i].x << " " << uniqueVertices[i].y << " " << uniqueVertices[i].z << std::endl;
		}
	}

	vtk << "POINT_DATA " << pointCount << std::endl;
	vtk << "VECTORS normals float" << std::endl;
	if (pointCount > 0) {
		for (auto&& v : *data) {
			vtk << v.x << " " << v.y << " " << v.z << std::endl;
		}
	}

	vtk.close();
}

size_t VTKExporter::countTriangulationIndices(std::vector<BufferGeom::Triangulation>& triangulations)
{
	std::vector<BufferGeom::Triangulation> T = triangulations;
	size_t nCellIds = 0;

	for (std::vector<BufferGeom::Triangulation>::iterator it = T.begin(); it < T.end(); ++it) {
		size_t currentSize = it->size();
		nCellIds += currentSize;
	}

	return nCellIds;
}
