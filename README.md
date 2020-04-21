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
#### Model:
![modelEqns](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/EvolutionModel.jpg)
#### Default real parameter values:
(see `./SurfaceEvolver/EvolutionCore/Parameters.h`)
- `C1 = 1.0, C2 = rDecay` (exp. decay radius of the initial sphere)
- `C = 1.0, D = 0.0`
- `S0 = 0.05 (initSmoothRate), lambda = 0.1 (smoothDecay)`

Smoothing model: Follows the above evolution equation with parameters:
![smoothEqns](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/SmoothingModel.jpg)

### Finite Volumes (cotan scheme)
![icoFVs](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/IcoSphereFVBuilding.gif)

#### Sphere Test:
```
================================
>>> Evolution error log ........
dt = 0.01, Nsteps = 6, Nverts = 42
L2Error = 0.00306691
dt = 0.0025, Nsteps = 24, Nverts = 162
L2Error = 0.000867556, EOC = 1.82176
dt = 0.000625, Nsteps = 96, Nverts = 642
L2Error = 0.000210643, EOC = 2.04216
dt = 0.00015625, Nsteps = 384, Nverts = 2562
L2Error = 5.1363e-05, EOC = 2.036
```

##### 
![sphereTest](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/ShrinkingSphere.gif)

# 1.3 Evolution in the Normal Direction - Using -grad(d(x)) to Control Evolution:
### (original model, remeshed without tangential redistribution,  after smoothing)
## Scaled Icosahedron
- **octreeResolution** = `40^3`, **SDF_gridResolution** = `120 x 106 x 106`
- **NTimeSteps** = `150`, **NSmoothingSteps** = `10`, **dt** = `0.03`
- **NVerts** = `162`
![icoEllipsoid](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/IcoSphereBasicRemesh.jpg)

## Cube
- **octreeResolution** = `40^3`, **SDF_gridResolution** = `120^3`
- **NTimeSteps** = `100`, **NSmoothingSteps** = `10`, **dt** = `0.03`
- **NVerts** = `162`
![cube](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/CubeBasicRemesh.jpg)

## How Does Evolution Behave For Non-Convex Target Meshes?

## Cube With Holes
- **octreeResolution** = `80^3`, **SDF_gridResolution** = `240 x 239 x 240`
- **NTimeSteps** = `130`, **dt** = `0.018`
- **NVerts** = `2562`

**Graphics annotation:**
- **vectors:** vertex normals (black), grad(distance) (red)
- **scalars:** dot(grad(distance), vNormal)
![cubeWHoles1](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/ToCubeWithHolesWDotProducts.gif)

## Cube With Holes (Larger Time Step)
- **octreeResolution** = `40^3`, **SDF_gridResolution** = `120 x 119 x 120`
- **NTimeSteps** = `100`, **dt** = `0.058`
- **NVerts** = `2562`

**Graphics annotation:**
- **vectors:** vertex normals (black), grad(distance) (red)
- **scalars:** dot(grad(distance), vNormal)
![cubeWHoles](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/CubeWithHolesLargeStep.gif)

## Arc
- **octreeResolution** = `40^3`, **SDF_gridResolution** = `120 x 115 x 96`
- **NTimeSteps** = `100`, **dt** = `0.013`
- **NVerts** = `2562`
### Left View:
![evolveArc1](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/evolUpdateSend1.gif)
### Right View (strange indentation):
![evolveArc2](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/evolUpdateSend2.gif)

## Stanford Bunny
- **octreeResolution** = `40^3`, **SDF_gridResolution** = `119 x 119 x 111`
- **NTimeSteps** = `100` (terminates after `59` steps on a degenerate triangle), **dt** = `0.03`
- **NVerts** = `2562`
![bunnyEvolGif](https://github.com/MCInversion/SurfaceEvolverDevelop/blob/master/SurfaceEvolver/Images/BunnyEvolution.gif)
