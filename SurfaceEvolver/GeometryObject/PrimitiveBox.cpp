#include "PrimitiveBox.h"

PrimitiveBox::PrimitiveBox()
{
}

PrimitiveBox::PrimitiveBox(const PrimitiveBox& other)
{
	Geometry::Geometry(other);
	this->dimensions[0] = other.dimensions[0];
	this->dimensions[2] = other.dimensions[1];
	this->dimensions[1] = other.dimensions[2];

	this->segments[0] = other.segments[0];
	this->segments[2] = other.segments[1];
	this->segments[1] = other.segments[2];
}

PrimitiveBox::PrimitiveBox(double x, double y, double z, unsigned int sx, unsigned int sy, unsigned int sz, bool quad, std::string name, bool lastWall)
{
	this->dimensions[0] = x;
	this->dimensions[2] = y;
	this->dimensions[1] = z;

	this->segments[0] = sx;
	this->segments[2] = sy;
	this->segments[1] = sz;

	this->quad = quad;
	this->lastWall = lastWall;

	if (name.empty()) {
		this->name = "PrimitiveBox, dimensions: (" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) +
			"), segments: (" + std::to_string(sx) + ", " + std::to_string(sy) + ", " + std::to_string(sz) + ")";
	}
	build();
}

PrimitiveBox::~PrimitiveBox()
{
}


void PrimitiveBox::build()
{
	// building the wall vertices according to the following disjoint 2D covering scheme:
	//           ___>
	//          | 1 |
	//       ...|___|...
	//      | 2 | 3 | 4 |
	//      |___|___|___|
	//          : 5 :
	//          :___:
	//          : 6 :
	//          :...:
	//
	// where 1 is the ground (xy) floor and wall 6 can be omitted

	const double width = dimensions[0];
	const double depth = dimensions[2];
	const double height = dimensions[1];

	const unsigned int sw = segments[0];
	const unsigned int sd = segments[2];
	const unsigned int sh = segments[1];

	const unsigned int nw = sw + 1;
	const unsigned int nd = sd + 1;

	const unsigned int vertexCount =
		nw * nd + // wall 1
		(2 * sd + nw) * sh + // walls 2, 3, and 4
		(sw - 1) * sd + // wall 5
		(sw - 1) * (sh - 1) * lastWall; // wall 6

	const unsigned int faceCount = 2 * (sw * sd + sh * sw + sd * sh);

	const unsigned int vertexBufferCount = 3 * vertexCount;
	const unsigned int normalBufferCount = 3 * 3 * 2 * faceCount;
	const unsigned int vertexIndexBufferCount = 3 * 2 * faceCount;

	vertices = std::vector<double>(vertexBufferCount);
	normals = std::vector<double>(normalBufferCount);
	vertexIndices = std::vector<unsigned int>(vertexIndexBufferCount);

	triangulations = std::vector<std::vector<unsigned int>>();
	unsigned int triId = 0;

	auto bufferQuadFromIds = [&](
		unsigned int i0, unsigned int i1, unsigned int i2, unsigned int i3,
		Vector3 normal
	) {
			vertexIndices[3 * triId] = i0;
			vertexIndices[3 * triId + 1] = i2;
			vertexIndices[3 * triId + 2] = i3;

			normals[9 * triId] = normal.x;
			normals[9 * triId + 1] = normal.y;
			normals[9 * triId + 2] = normal.z;

			normals[9 * triId + 3] = normal.x;
			normals[9 * triId + 4] = normal.y;
			normals[9 * triId + 5] = normal.z;

			normals[9 * triId + 6] = normal.x;
			normals[9 * triId + 7] = normal.y;
			normals[9 * triId + 8] = normal.z;

			vertexIndices[3 * triId + 3] = i0;
			vertexIndices[3 * triId + 4] = i1;
			vertexIndices[3 * triId + 5] = i2;

			normals[9 * (triId + 1)] = normal.x;
			normals[9 * (triId + 1) + 1] = normal.y;
			normals[9 * (triId + 1) + 2] = normal.z;

			normals[9 * (triId + 1) + 3] = normal.x;
			normals[9 * (triId + 1) + 4] = normal.y;
			normals[9 * (triId + 1) + 5] = normal.z;

			normals[9 * (triId + 1) + 6] = normal.x;
			normals[9 * (triId + 1) + 7] = normal.y;
			normals[9 * (triId + 1) + 8] = normal.z;

			if (quad) triangulations.push_back({ triId, triId + 1 });
			else {
				triangulations.push_back({ triId });
				triangulations.push_back({ triId + 1 });
			}
	};

	Vector3 currentNormal = Vector3(0, 0, -1);

	// wall 1:
	for (unsigned int i = 0; i <= sd; i++) {
		for (unsigned int j = 0; j <= sw; j++) {
			vertices[3 * (i * nw + j)] = (1 - (double)j / (double)sw) * width;
			vertices[3 * (i * nw + j) + 1] = (1 - (double)i / (double)sd) * depth;
			vertices[3 * (i * nw + j) + 2] = 0;

			// floor (wall 1) faces:
			if (i > 0 && j > 0) {
				bufferQuadFromIds (
					i * nw + (j - 1),
					i * nw + j,
					(i - 1) * nw + j,
					(i - 1) * nw + (j - 1),
					currentNormal
				);
				triId += 2;
			}
		}
	}

	unsigned int startingId = nw * nd;
	unsigned int rowLength = 2 * sd + nw;

	// walls 2, 3, and 4:
	for (unsigned int i = 0; i < sh; i++) {
		// wall 2:
		currentNormal.set(1, 0, 0);
		for (unsigned int j = 0; j <= sd; j++) {
			vertices[3 * (startingId + i * rowLength + j)] = width;
			vertices[3 * (startingId + i * rowLength + j) + 1] = (1 - (double)j / (double)sd) * depth;
			vertices[3 * (startingId + i * rowLength + j) + 2] = (double)(i + 1) / (double)sh * height;

			// wall 2 faces: first row:
			if (i == 0 && j > 0) {
				bufferQuadFromIds (
					startingId + j - 1,
					startingId + j,
					j * nw,
					(j - 1) * nw,
					currentNormal
				);
				triId += 2;
			}
			// wall 2 faces: next rows:
			else if (i > 0 && j > 0) {
				bufferQuadFromIds(
					startingId + i * rowLength + j - 1,
					startingId + i * rowLength + j,
					startingId + (i - 1) * rowLength + j,
					startingId + (i - 1) * rowLength + j - 1,
					currentNormal
				);
				triId += 2;
			}
		}
		// wall 3:
		currentNormal.set(0, -1, 0);
		for (unsigned int j = 1; j <= sw; j++) {
			vertices[3 * (startingId + i * rowLength + sd + j)] = (1 - (double)j / (double)sw) * width;
			vertices[3 * (startingId + i * rowLength + sd + j) + 1] = 0;
			vertices[3 * (startingId + i * rowLength + sd + j) + 2] = (double)(i + 1) / (double)sh * height;

			// wall 3 faces: first row:
			if (i == 0) {
				bufferQuadFromIds(
					startingId + sd + j - 1,
					startingId + sd + j,
					startingId - nw + j,
					startingId - nw + j - 1,
					currentNormal
				);
				triId += 2;
			}
			// wall 3 faces: next rows:
			else if (i > 0) {
				bufferQuadFromIds(
					startingId + sd + i * rowLength + j - 1,
					startingId + sd + i * rowLength + j,
					startingId + sd + (i - 1) * rowLength + j,
					startingId + sd + (i - 1) * rowLength + j - 1,
					currentNormal
				);
				triId += 2;
			}
		}
		// wall 4:
		currentNormal.set(-1, 0, 0);
		for (unsigned int j = 1; j <= sd; j++) {
			vertices[3 * (startingId + i * rowLength + sd + sw + j)] = 0;
			vertices[3 * (startingId + i * rowLength + sd + sw + j) + 1] = (double)j / (double)sd * depth;
			vertices[3 * (startingId + i * rowLength + sd + sw + j) + 2] = (double)(i + 1) / (double)sh * height;

			// wall 4 faces: first row:
			if (i == 0) {
				bufferQuadFromIds(
					startingId + sd + sw + j - 1,
					startingId + sd + sw + j,
					startingId - j * nw - 1,
					startingId - (j - 1) * nw - 1,
					currentNormal
				);
				triId += 2;
			}
			// wall 4 faces: next rows:
			else if (i > 0) {
				bufferQuadFromIds(
					startingId + sd + sw + i * rowLength + j - 1,
					startingId + sd + sw + i * rowLength + j,
					startingId + sd + sw + (i - 1) * rowLength + j,
					startingId + sd + sw + (i - 1) * rowLength + j - 1,
					currentNormal
				);
				triId += 2;
			}
		}
	}

	startingId += rowLength * sh;

	// wall 5:
	currentNormal.set(0, 0, 1);
	for (unsigned int i = 1; i <= sd; i++) {

		if (sw > 1) {
			for (unsigned int j = 1; j < sw; j++) {
				vertices[3 * (startingId + (i - 1) * (sw - 1) + j - 1)] = (1 - (double)j / (double)sw) * width;
				vertices[3 * (startingId + (i - 1) * (sw - 1) + j - 1) + 1] = (double)i / (double)sd * depth;
				vertices[3 * (startingId + (i - 1) * (sw - 1) + j - 1) + 2] = height;

				// ceiling wall (5) face: upper left corner:
				if (i == 1 && j == 1) {
					bufferQuadFromIds(
						startingId - rowLength + sd - 1,
						startingId,
						startingId - rowLength + sd + 1,
						startingId - rowLength + sd,
						currentNormal
					);
					triId += 2;
				}
				// ceiling wall (5) faces: left edge:
				else if (i > 1 && j == 1) {
					bufferQuadFromIds(
						startingId - rowLength + sd - i,
						startingId + (i - 1) * (sw - 1),
						startingId + (i - 2) * (sw - 1),
						startingId - rowLength + sd - i + 1,
						currentNormal
					);
					triId += 2;
				}
				// ceiling wall (5) faces: first row:
				else if (i == 1) {
					bufferQuadFromIds(
						startingId + j - 2,
						startingId + j - 1,
						startingId - sw - sd + j - 1,
						startingId - sw - sd + (j - 2),
						currentNormal
					);
					triId += 2;
				}
				// ceiling wall (5) interior faces:
				else if (i > 1 && j > 1) {
					bufferQuadFromIds(
						startingId + (i - 1) * (sw - 1) + j - 2,
						startingId + (i - 1) * (sw - 1) + j - 1,
						startingId + (i - 2) * (sw - 1) + j - 1,
						startingId + (i - 2) * (sw - 1) + j - 2,
						currentNormal
					);
					triId += 2;
				}
			}
			// ceiling wall (5) face: upper right corner:
			if (i == 1) {
				bufferQuadFromIds(
					startingId + sw - 2,
					startingId - nd + 1,
					startingId - nd,
					startingId - nd - 1,
					currentNormal
				);
				triId += 2;
			}
			// ceiling wall (5) faces: right edge:
			else {
				bufferQuadFromIds(
					startingId + i * (sw - 1) - 1,
					startingId - nd + i,
					startingId - nd + i - 1,
					startingId + (i - 1) * (sw - 1) - 1,
					currentNormal
				);
				triId += 2;
			}
			// case for sw = 1:
		}
		else {
			if (i == 1) {
				bufferQuadFromIds(
					startingId - rowLength + sd - 1,
					startingId - sd,
					startingId - nd,
					startingId - nd - 1,
					currentNormal
				);
				triId += 2;
			}
			else {
				bufferQuadFromIds(
					startingId - rowLength + sd - i,
					startingId - nd + i,
					startingId - nd + i - 1,
					startingId - rowLength + sd - i + 1,
					currentNormal
				);
				triId += 2;
			}
		}
	}

	startingId += (sw - 1) * sd;

	if (lastWall) {
	// wall 6:
		currentNormal.set(0, 1, 0);
		for (unsigned int i = 1; i < sh; i++) {
			if (sw > 1) {
				for (unsigned int j = 1; j < sw; j++) {
					vertices[3 * (startingId + (i - 1) * (sw - 1) + j - 1)] = (1 - (double)j / (double)sw) * width;
					vertices[3 * (startingId + (i - 1) * (sw - 1) + j - 1) + 1] = depth;
					vertices[3 * (startingId + (i - 1) * (sw - 1) + j - 1) + 2] = (1 - (double)i / (double)sh) * height;

					// wall 6 face: upper left corner:
					if (i == 1 && j == 1) {
						bufferQuadFromIds(
							startingId - (sw - 1) * sd - 2 * rowLength,
							startingId,
							startingId - (sw - 1),
							startingId - (sw - 1) * sd - rowLength,
							currentNormal
						);
						triId += 2;
					}
					// wall 6 faces: first row interior:
					else if (i == 1 && j > 1) {
						bufferQuadFromIds(
							startingId + j - 2,
							startingId + j - 1,
							startingId - (sw - 1) + j - 1,
							startingId - (sw - 1) + j - 2,
							currentNormal
						);
						triId += 2;
					}
					// wall 6 faces: left edge:
					else if (i > 1 && j == 1) {
						bufferQuadFromIds(
							startingId - (sw - 1) * sd - 2 * rowLength - (i - 1) * rowLength,
							startingId + (i - 1) * (sw - 1),
							startingId + (i - 2) * (sw - 1),
							startingId - (sw - 1) * sd - i * rowLength,
							currentNormal
						);
						triId += 2;
					}
					// wall 6 faces: interior:
					else if (i > 1 && j > 1) {
						bufferQuadFromIds(
							startingId + (i - 1) * (sw - 1) + j - 2,
							startingId + (i - 1) * (sw - 1) + j - 1,
							startingId + (i - 2) * (sw - 1) + j - 1,
							startingId + (i - 2) * (sw - 1) + j - 2,
							currentNormal
						);
						triId += 2;
					}
				}
				// wall 6 face: upper right corner:
				if (i == 1) {
					bufferQuadFromIds(
						startingId + (sw - 2),
						startingId - (sw - 1) * sd - rowLength - 1,
						startingId - (sw - 1) * sd - 1,
						startingId - 1,
						currentNormal
					);
					triId += 2;
				}
				// wall 6 faces: right edge:
				else {
					bufferQuadFromIds(
						startingId + i * (sw - 1) - 1,
						startingId - (sw - 1) * sd - i * rowLength - 1,
						startingId - (sw - 1) * sd - (i - 1) * rowLength - 1,
						startingId + (i - 1) * (sw - 1) - 1,
						currentNormal
					);
					triId += 2;
				}
				// case for sw = 1: wall 6 faces:
			}
			else {
				bufferQuadFromIds(
					startingId - (sw - 1) * sd - 2 * rowLength - (i - 1) * rowLength,
					startingId - (sw - 1) * sd - i * rowLength - 1,
					startingId - (sw - 1) * sd - (i - 1) * rowLength - 1,
					startingId - (sw - 1) * sd - i * rowLength,
					currentNormal
				);
				triId += 2;
			}
		}

		if (sw > 1) {
			// wall 6 face: lower left corner:
			bufferQuadFromIds(
				0,
				1,
				startingId + (sw - 1) * (sh - 2),
				nw * nd,
				currentNormal
			);
			triId += 2;
			// wall 6 faces: last row interior:
			for (unsigned int j = 1; j < sw - 1; j++) {
				bufferQuadFromIds(
					j,
					j + 1,
					startingId + (sw - 1) * (sh - 2) + j,
					startingId + (sw - 1) * (sh - 2) + j - 1,
					currentNormal
				);
				triId += 2;
			}
			// wall 6 face: lower right corner:
			bufferQuadFromIds(
				sw - 1,
				sw,
				nw * nd + rowLength - 1,
				vertexCount - 1,
				currentNormal
			);
			// case for sw = 1: wall 6 last face:
		}
		else {
			bufferQuadFromIds(
				0,
				sw,
				nw * nd + rowLength - 1,
				nw * nd,
				currentNormal
			);
		}
	}	

	// duplicate vertices into geometryVertices
	std::vector<double> geometryVertices = std::vector<double>(3 * this->vertexIndices.size());

	for (unsigned int i = 0; i < this->vertexIndices.size(); i++) {
		geometryVertices[i * 3] = this->vertices[this->vertexIndices[i] * 3];
		geometryVertices[i * 3 + 1] = this->vertices[this->vertexIndices[i] * 3 + 1];
		geometryVertices[i * 3 + 2] = this->vertices[this->vertexIndices[i] * 3 + 2];
	}

	// copy unique vertex coords into uniqueVertices
	this->uniqueVertices = std::vector<Vector3>(this->getVertices());
	this->vertices = std::vector<double>(geometryVertices);
}
