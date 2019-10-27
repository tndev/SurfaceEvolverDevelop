#ifndef AABBTREE_H
#define AABBTREE_H

#include "Geometry.h"
#include "Vector3.h"
#include "AABBNode.h"

template<class T>
class AABBTree
{
public:
	AABBNode<T> root;
	unsigned int depth = 0;
	unsigned int nLeaves = 0;

	AABBTree();
	AABBTree(std::vector<T>& objects, unsigned int depthLeft);
};


#endif
