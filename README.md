# SurfaceEvolver

The goal of this app is to re-mesh an arbitrary input triangle mesh via a [method of evolving surfaces with tangential redistribution](http://www.math.sk/mikula/mrss_SISC.pdf) (area based, angle based...).

Done: generate a distance field field around a mesh using the [Fast Sweeping Method](https://graphics.stanford.edu/courses/cs468-03-fall/Papers/zhao_fastsweep1.pdf). The initial condition for the FSM is set by constructing an AABB tree for the mesh triangles (or any mesh features in general), using it for faster lookup during Octree box cell intersection queries. The Octree leaves are then set to the distance to the closest triangle to the cell centroid.

Current WIP: 

- optimize the construction of AABBTree and Octree. Try flat node array
- evolve a simple primitive onto an object contained within (e.g.: `CubeSphere` -> `Cube` )

### ============ Progress So far ===============:

# 1.1 Evolution In the Normal Direction - Distance Function

#### 1. Take the input geometry & construct an AABB Tree
![AABBFull](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/BunnyAABBNodes.jpg)

#### 2. Using the AABB Tree, construct an Octree with homogeneous disjoint leaf cells which intersect the mesh

![OctreeFull](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/BunnyOctreeFull.jpg)

#### 3. Construct a grid and set the value of all cells which intersect the mesh to exact distance from cell centroid to the closest triangle and `LARGE_VAL` elsewhere:
![OctreeLeafCells](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/BunnyOctreeLeafCells.jpg)

#### 4. Using the grid from step (3) as an initial condition for the Eikonal equation `|grad(d(x))| = 1` find its solution (distance function) using the Fast Sweeping Method:
![FastSweep](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/FS_resized.gif)

### Result:
![DF](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/BunnySDF_FS.jpg)
### Error :
|MyMethodGrid(x,y,z) - BruteForceGrid(x,y,z)|
![error](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/BunnySDF_FS_Error.jpg)

#### Sign is computed by negating the grid field d(x,y,z) and flood filling "unfrozen" voxels to set external sign > 0:
![SignComp](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/BunnySDF_Sign.jpg)

- clearly, if the mesh has holes, the flood fill will set all non-boundary voxels to external

#### Evolution will be carried out in the direction of `-grad(SDF(x,y,z))`:
![bunnyDirection](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/EvolutionInBunyDirection3D.jpg)

# 1.2 Evolution in the Normal Direction - Finite Volume Scheme for Laplace-Beltrami Operator and Mean Curvature Flow
### Finite Volumes (cotan scheme):
![icoFVs](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/IcoSphereFVBuilding.gif)
#### Sphere Test:
![sphereTest](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/ShrinkingSphere.gif)
