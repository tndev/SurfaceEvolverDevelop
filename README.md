# SurfaceEvolver

The goal of this app is to re-mesh an arbitrary input triangle mesh via a [method of evolving surfaces with tangential redistribution](http://www.math.sk/mikula/mrss_SISC.pdf) (area based, angle based...).

Done: generate a distance field field around a mesh using the [Fast Sweeping Method](https://graphics.stanford.edu/courses/cs468-03-fall/Papers/zhao_fastsweep1.pdf). The initial condition for the FSM is set by constructing an AABB tree for the mesh triangles (or any mesh features in general), using it for faster lookup during Octree box cell intersection queries. The Octree leaves are then set to the distance to the closest triangle to the cell centroid.

Current WIP: optimize the construction of AABBTree and Octree. Try flat node array + dynamic update under small transformations of the mesh

### 1. Take the input geometry & construct an AABB Tree
![AABBFull](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/BunnyAABBNodes.jpg)

### 2. Using the AABB Tree, construct an Octree with homogeneous disjoint leaf cells which intersect the mesh

![OctreeFull](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/BunnyOctreeFull.jpg)

### 3. Construct a grid and set the value of all cells which intersect the mesh to exact distance from cell centroid to the closest triangle and `LARGE_VAL` elsewhere:
![OctreeLeafCells](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/BunnyOctreeLeafCells.jpg)

### 4. Using the grid from step (3) as an initial condition for the Eikonal equation `|grad(d(x))| = 1` find its solution (distance function) using the Fast Sweeping Method:
![FastSweep](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/FS_resized.gif)

## Result:
![DF](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/bunnySDFContour1.jpg)
### (mean) L2 Error:
![error](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/BunnySDF_FS_Error.jpg)

