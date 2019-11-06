# SurfaceEvolver

The goal of this app is to re-mesh an arbitrary input triangle mesh via a [method of evolving surfaces with tangential redistribution](http://www.math.sk/mikula/mrss_SISC.pdf) (area based, angle based...).

Current WIP: generate an SDF field around a mesh using the [Fast Sweeping Method](https://graphics.stanford.edu/courses/cs468-03-fall/Papers/zhao_fastsweep1.pdf). I already have a prototype for an AABBTree structure which speeds up triangle lookup for an Octree, so one can set the proper initial condition near the mesh:

![Voxelization](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/BunnyAABBToOctree.jpg?raw=true)

## Result:

![SDF](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/BunnySDF_2.jpg?raw=true)
