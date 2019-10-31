#include "Deform.h"

Deform::Deform()
{
}

Deform::Deform(Geometry* geom)
{
	this->geom = geom;
}

Deform::~Deform()
{
}

void Deform::spherify(float param)
{
	result = Geometry(*geom);
	Box3 bbox = geom->getBoundingBox();
	Vector3 center = bbox.getCenter();
	float radius = 0.5 * bbox.getSize().length();

	Vector3 vertex = Vector3();
	Vector3 targetVertex = Vector3();
	Vector3 normal = Vector3();
	Vector3 targetNormal = Vector3();

	Vector3 helperVector = Vector3();

	for (unsigned int i = 0; i < geom->vertices.size(); i += 3) {
		vertex.set(geom->vertices[i], geom->vertices[i + 1], geom->vertices[i + 2]);
		normal.set(geom->normals[i], geom->normals[i + 1], geom->normals[i + 2]);

		helperVector = vertex;

		targetNormal = normalize(vertex - center);
		if (!targetNormal.lengthSq()) {
			targetNormal = normal;
		}

		targetVertex = radius * targetNormal + center;
		vertex.lerp(targetVertex, param);
		normal.lerp(targetNormal, param);
		if (!normal.lengthSq()) {
			normal.set(geom->normals[i], geom->normals[i + 1], geom->normals[i + 2]);
		}

		result.vertices[i] = vertex.x;
		result.vertices[i + 1] = vertex.y;
		result.vertices[i + 2] = vertex.z;

		result.normals[i] = normal.x;
		result.normals[i + 1] = normal.y;
		result.normals[i + 2] = normal.z;

		result.uniqueVertices[geom->vertexIndices[i / 3]] = vertex;
	}
}


