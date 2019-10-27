#ifndef AABBNODE_H_
#define AABBNODE_H_

#include <vector>
#include "Geometry.h"
#include "Box3.h"

template<class T>
class AABBNode
{
public:
	AABBNode<T>* left;
	AABBNode<T>* right;
	Box3 bbox;
	std::vector<T> content; // could contain vertices, edges, or polygons

	AABBNode();
	~AABBNode();
	AABBNode(std::vector<T>& objects, unsigned int depthLeft);
};

#endif
