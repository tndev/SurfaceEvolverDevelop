# SurfaceEvolver

The goal of this app is to re-mesh an arbitrary input triangle mesh via a [method of evolving surfaces with tangential redistribution](http://www.math.sk/mikula/mrss_SISC.pdf) (area based, angle based...).

Done: generate a distance field field around a mesh using the [Fast Sweeping Method](https://graphics.stanford.edu/courses/cs468-03-fall/Papers/zhao_fastsweep1.pdf).

Current WIP: optimize the construction of AABBTree and Octree. Try flat node array + dynamic update under small transformations of the mesh

![Voxelization](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/BunnyAABBToOctree.jpg)

## Result:

![DF](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/CorrectBunnyDist.jpg)

### Cube:

![cubeDF](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/CorrectCubeDist.jpg)
