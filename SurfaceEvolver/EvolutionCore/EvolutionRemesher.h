#ifndef EVOLUTIONREMESHER_H_
#define EVOLUTIONREMESHER_H_

#include "../Geometry/Geometry.h"
#include "../SDF/SDF.h"
#include "Parameters.h"
#include "Evolver.h"

class EvolutionRemesher
{
public:
	EvolutionRemesher();
	EvolutionRemesher(Geometry* geom);

private:
	Geometry* geom;
	SDF sdf; // signed distance function

	ElementType eType;
	EvolutionParams eParams;
	SphereTestParams stParams;
	MeanCurvatureParams mcfParams;
	GradDistanceParams sdfParams;

};

#endif

