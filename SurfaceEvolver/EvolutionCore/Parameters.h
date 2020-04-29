#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <string>
#include "../Geometry/Geometry.h"
#include "../SDF/Grid.h"

enum class ElementType {
	tri = 0,
	quad = 1
};

struct EvolutionParams {
	float dt = 0.01f;
	float tStop = 1.0f;
	int NSteps = 10;
	unsigned int subdiv = 2;
	ElementType elType = ElementType::tri;

	std::string name = "Sphere";

	Geometry* sourceGeometry = nullptr;

	// ===== output flags ==========
	bool saveStates = false;
	bool saveResult = true;
	bool printHappenings = false;
	bool printStepOutput = false;
	bool printSolution = false;

	bool writeGenericLog = true;
	bool writeTimeLog = false;
	// =============================
};

struct SphereTestParams {
	float r0 = 1.0f;
	int testId = -1;

	// ===== output flags ==========
	bool writeErrorLog = false;
	// =============================
};

struct MeanCurvatureParams {
	// Laplace-Beltrami ctrl constants:
	float rDecay = 1.0f;
	float C1 = 1.0f;
	float C2 = rDecay;

	bool constant = false;

	// ===== output flags ==========
	bool saveAreaStates = false;
	bool saveCurvatureStates = false;
	bool saveNormalVelocityStates = false;
	bool saveCurvatureVectors = false;

	bool writeMeanAreaLog = false;
	// =============================

	int smoothSteps = -1; // -1 - no smooth steps
	float initSmoothRate = 0.05f;
	float smoothDecay = 0.1f;
};

struct GradDistanceParams {
	Geometry* targetGeom = nullptr;
	Grid* sdfGrid = nullptr;

	// grad SDF ctrl constants:
	float C = 1.0f;
	float D = 0.0f;

	bool constant = false;

	// ===== output flags ==========
	bool saveDistanceStates = false;
	bool saveGradientStates = false;
	// =============================
};

struct TangentialRedistParams {
	// ctrl params
	float omega_volume = 100.0f;
	float omega_angle = 2.0f;

	// ===== Redistribution type =======
	int type = 0; // -1 - none, 0 - angle-based, 1 - volume-based, 2 - length-based

	// ===== output flags ======
	bool saveTangentialVelocityStates = false;
	// ================
};

#endif
