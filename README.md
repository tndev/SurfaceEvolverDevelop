# SurfaceEvolver

The goal of this app is to re-mesh an arbitrary input triangle mesh via a [method of evolving surfaces with tangential redistribution](http://www.math.sk/mikula/mrss_SISC.pdf) (area based, angle based...).

Done: generate a distance field field around a mesh using the [Fast Sweeping Method](https://graphics.stanford.edu/courses/cs468-03-fall/Papers/zhao_fastsweep1.pdf). The initial condition for the FSM is set by constructing an AABB tree for the mesh triangles (or any mesh features in general), using it for faster lookup during Octree box cell intersection queries. The Octree leaves are then set to the distance to the closest triangle to the cell centroid.

Current WIP: optimize the construction of AABBTree and Octree. Try flat node array + dynamic update under small transformations of the mesh

![Voxelization](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/BunnyAABBToOctree.jpg){:height="700px" width="400px"}

![FastSweep](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/FastSweepResized.gif){:height="700px" width="400px"}

## Result:

![DF](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/bunnySDFContour1.jpg){:height="700px" width="400px"}


