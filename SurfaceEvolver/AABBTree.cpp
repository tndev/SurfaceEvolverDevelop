#include "AABBTree.h"

AABBTree::AABBTree()
{
}

AABBTree::AABBTree(std::vector<StructGeom::Triangle>& triangles, unsigned int depthLeft)
{
	this->buildNode(triangles, depthLeft);
}

AABBTree::~AABBTree()
{
}

void AABBTree::buildNode(std::vector<StructGeom::Triangle>& triangles, unsigned int depthLeft)
{
}
