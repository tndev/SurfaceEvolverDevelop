#include "Evolver.h"


void Evolver::init()
{
	// ===== Preparing a source geometry for immersion =====
	// -----------------------------------------------------
	// radius
	float r = (meanCurvatureFlow ? r0 : 0.4f * std::min({ sdfGrid->scale.x, sdfGrid->scale.y, sdfGrid->scale.z }));
	rDecay = r;

	if (elType == ElementType::tri) {
		uint n = subdiv;
		evolvedSurface = new IcoSphere(n, r);
	}
	else {
		uint n = std::ceil(1.5 * subdiv);
		evolvedSurface = new CubeSphere(n, r);
	}

	Box3 bbox = (meanCurvatureFlow || targetGeom != nullptr ? targetGeom->getBoundingBox() : sdfGrid->bbox);
	this->center = bbox.getCenter();
	Matrix4 M = Matrix4().makeTranslation(center.x, center.y, center.z);
	evolvedSurface->applyMatrix(M);
	this->N = evolvedSurface->uniqueVertices.size();

	log_header.clear();
	log_header += "=============================================================\n";
	log_header += "----------------------- Initializing ------------------------\n";
	log_header += "------------ S U R F A C E    E V O L U T I O N -------------\n";
	log_header += ">> NSteps = " + std::to_string(NSteps) + ", dt = " + std::to_string(dt) + ", tStop = " + std::to_string(tStop) + "\n";
	log_header += ">> Model: " + geomName + ", NVerts = " + std::to_string(N) + "\n";
	if (!meanCurvatureFlow) {
		log_header += ">> SDF grid bounds: min = (" + std::to_string(sdfGrid->bbox.min.x) + ", " + std::to_string(sdfGrid->bbox.min.y) + ", " + std::to_string(sdfGrid->bbox.min.z)
			+ "), max = (" + std::to_string(sdfGrid->bbox.max.x) + ", " + std::to_string(sdfGrid->bbox.max.y) + ", " + std::to_string(sdfGrid->bbox.max.z) + ")\n";
		log_header += ">> SDF resolution: " + std::to_string(sdfGrid->Nx) + " x " + std::to_string(sdfGrid->Ny) + " x " + std::to_string(sdfGrid->Nz) + "\n";
	}
	else {
		log_header += ">> Performing mean curvature flow test on sphere with r0 = " + std::to_string(r0) + "\n";
	}

	log_header += "-------------------------------------------------------------\n";
	log_header += "starting evolution ...\n";


	std::cout << log_header;

	if (saveStates) {
		if (!sphereTest) {
			evolvedSurface->clearScalarData();
			// save first distance state:
			if (saveDistanceStates) saveInterpolatedSDFValues();
			// save first fvState:
			if (saveAreaStates) saveFVAreaScalars();
			// save first dot(-grad(SDF), N)  state
			if (saveGradientStates) saveInterpolatedDotValues();
			// save first curvature state:
			if (saveCurvatureStates || tangential_redist) {
				computeSurfaceNormalsAndCoVolumes();
				getCurvaturesAndNormalVelocities();
				if (saveCurvatureStates) saveMeanCurvatureScalars();
			}
		}
		exportGeometry(0);
	}

	// init solvers
	solveX = new LinearSolver(N);
	solveY = new LinearSolver(N);
	solveZ = new LinearSolver(N);

	// -------------------------------------------------------
}

void Evolver::evolve()
{
	if (!meanCurvatureFlow) {
		if (!sdfGrid->hasGradient()) {
			sdfGrid->computeGradient();
		}
	}

	sphereTestL2Error = 0.0f;

	if (start_flag) {
		startGlobalTime = std::chrono::high_resolution_clock::now();
		start_flag = false;
	}

	uint TSteps = (smoothing_flag ? smoothSteps : NSteps);

	for (int ti = tBegin; ti < tBegin + TSteps; ti++) {

		float t = (float)ti / (float)NSteps * tStop;
		std::string timeStepLine = "-------- time step ti = " + std::to_string(ti) + ", t = " + std::to_string(t) + " -----------\n";
		if (printStepOutput) std::cout << timeStepLine;
		if (writeGenericLog) log << timeStepLine;

		// clear scalar data from previous step
		if (saveAreaStates || saveDistanceStates || saveCurvatureStates) evolvedSurface->clearScalarData();


		// <===== P S E U D O N O R M A L S   &   C O - V O L U M E S ========
		auto startNandCoVols = std::chrono::high_resolution_clock::now();
		if (printHappenings) std::cout << "computing surface pseudonormals and finite volumes..." << std::endl;

		this->computeSurfaceNormalsAndCoVolumes();

		if (saveCurvatureStates || tangential_redist) this->getCurvaturesAndNormalVelocities();
		if (saveCurvatureStates) this->saveMeanCurvatureScalars();

		if (printHappenings) std::cout << "... computed" << std::endl;
		auto endNandCoVols = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedNandCoVols = (endNandCoVols - startNandCoVols);
		// ===================================================================>



		this->initSystem(); // prep linear system of dim NVerts * NVerts



		// <===== M A T R I X    F I L L ==========
		auto startFillMatrix = std::chrono::high_resolution_clock::now();

		if (printHappenings) std::cout << "filling sys matrix & rhs ..." << std::endl;
		float meanArea;
		if (this->elType == ElementType::tri) {
			this->getTriangleEvolutionSystem(ti - tBegin, meanArea);
		}
		else {
			this->getQuadEvolutionSystem(ti - tBegin, meanArea);
		}

		if (degenerateCount > 0) {
			degenerateCount /= 3;
			std::cout << "WARNING: DETECTED " << degenerateCount << " DEGENERATE ELEMENTS!\n";
		}
		if (std::isnan(meanArea)) {
			std::cout << "meanArea = NaN!\n";
			break;
		}
		if (printHappenings) {
			std::cout << "... matrix & rhs filled" << std::endl;
			std::cout << "finite volume mean area: " << meanArea << std::endl;
		}
		if (writeGenericLog) log << "finite volume mean area: " << meanArea << std::endl;
		if (writeMeanAreaLog) meanAreaLog << meanArea << std::endl;

		auto endFillMatrix = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedFillMatrix = (endFillMatrix - startFillMatrix);
		// ========================================>



		// <===== S A V E   S D F  S T A T E ==============
		double saveVectTime = 0.0;
		if (saveGradientStates) {
			auto startSaveVectors = std::chrono::high_resolution_clock::now();

			exportVectorStates(ti - 1);
			saveInterpolatedDotValues();
			auto endSaveVectors = std::chrono::high_resolution_clock::now();
			std::chrono::duration<float> elapsedSaveVectors = (endSaveVectors - startSaveVectors);
			saveVectTime = elapsedSaveVectors.count();
		}
		// ================================================>




		// ------------ sys alloc ------------
		double* FnextX = new double[N];
		double* FnextY = new double[N];
		double* FnextZ = new double[N];
		// -----------------------------------



		// <===== L I N   S Y S   S O L V E ========
		auto startLinSolve = std::chrono::high_resolution_clock::now();

		if (printHappenings) std::cout << "3 x " << N << " x " << N << "-sys Bi CGSTAB solve ..." << std::endl;
		solveX->Bi_CGSTAB_Solve(SysMatrix, sysRhsX, FnextX, printSolution);
		solveY->Bi_CGSTAB_Solve(SysMatrix, sysRhsY, FnextY, printSolution);
		solveZ->Bi_CGSTAB_Solve(SysMatrix, sysRhsZ, FnextZ, printSolution);
		if (printSolution) {
			solveX->printArray1("Fx", FnextX, 6);
			solveY->printArray1("Fy", FnextY, 6);
			solveZ->printArray1("Fz", FnextZ, 6);
		}

		if (printHappenings) std::cout << "... done" << std::endl;

		auto endLinSolve = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedLinSolve = (endLinSolve - startLinSolve);
		// ========================================>



		// <====== G E O M   U P D A T E ==========
		auto startUpdate = std::chrono::high_resolution_clock::now();

		if (printHappenings) std::cout << "updating surface geometry ..." << std::endl;
		updateGeometry(FnextX, FnextY, FnextZ);
		if (printHappenings) std::cout << "... done" << std::endl;

		if (interruptEvolution) {
			std::cout << "ERROR! Solution explosion outside of SDF bounds! - INTERRUPTING EVOLUTION\n";
			break;
		}

		if (saveStates) {
			std::cout << "exporting step to VTK..." << std::endl;
			exportGeometry(ti);
			if (sphereTest) exportTestGeometry(ti, t);
			std::cout << "... done" << std::endl;
		}

		if (sphereTest) {
			float stepErrorSq = this->getSphereStepL2Error(t);
			float stepError = sqrt(stepErrorSq);

			std::string errorLine = "step " + std::to_string(ti) + " error: " + std::to_string(stepError) + "\n";
			if (printStepOutput) std::cout << errorLine;
			if (writeGenericLog) log << errorLine;
			if (writeErrorLog) errorLog << errorLine;

			sphereTestL2Error += stepErrorSq;
		}

		auto endUpdate = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedUpdate = (endUpdate - startUpdate);
		// ========================================>

		std::string time_log_line =
			"step " + std::to_string(ti) + " time log:\n" +
			"surface pseudonormals & finite volumes: " + std::to_string(elapsedNandCoVols.count()) + " s\n" +
			(saveGradientStates ? "export normals and gradients: " + std::to_string(saveVectTime) + "s\n" : "") +
			"filling sys matrix: " + std::to_string(elapsedFillMatrix.count()) + " s\n" +
			"(vector) linear sys solve: " + std::to_string(elapsedLinSolve.count()) + " s\n" +
			"geom update: " + std::to_string(elapsedUpdate.count()) + " s\n" +
			"evolution step total: " + std::to_string(elapsedNandCoVols.count() + elapsedFillMatrix.count() + elapsedLinSolve.count() + elapsedUpdate.count()) + " s\n";
		if (printStepOutput) std::cout << time_log_line << std::endl;
		if (writeTimeLog) timingLog << time_log_line;

		// ========= E N D =========================
	}

	if (saveResult) {
		exportResultGeometry((smoothing_flag ? "_smoothed" : "_result"));
	}

	if (sphereTest) {
		sphereTestL2Error *= dt;
		sphereTestL2Error = sqrt(sphereTestL2Error);

		std::string total_error_line = "total L2 Error: " + std::to_string(sphereTestL2Error) + "\n";

		if (printStepOutput) std::cout << total_error_line;
		if (writeGenericLog) log << total_error_line;
		if (writeErrorLog) errorLog << total_error_line;
	}

	// save last distance gradient state:
	if (saveGradientStates) {
		saveInterpolatedSDFGradients();
		exportVectorStates(NSteps);
	}

	if (smoothSteps < 0 || smoothing_flag) {
		auto endGlobalTime = std::chrono::high_resolution_clock::now();
		elapsedTimeTotal = (endGlobalTime - startGlobalTime);
		std::string total_time_message = "\n Total time for " + std::to_string(NSteps + std::max(0, smoothSteps)) + " steps: " + std::to_string(elapsedTimeTotal.count()) + " s\n";
		if (printStepOutput) std::cout << total_time_message;
		if (writeGenericLog) log << total_time_message;
		if (writeTimeLog) timingLog << total_time_message;
	}
}

void Evolver::clearSystem()
{
	if (SysMatrix != nullptr) {
		for (uint i = 0; i < N; i++) delete[] SysMatrix[i];
		delete[] SysMatrix;
	}
	if (sysRhsX != nullptr) delete[] sysRhsX;
	if (sysRhsY != nullptr) delete[] sysRhsY;
	if (sysRhsZ != nullptr) delete[] sysRhsZ;

	SysMatrix = nullptr;
	sysRhsX = nullptr;
	sysRhsY = nullptr;
	sysRhsZ = nullptr;
}

void Evolver::initSystem()
{
	// ===== Preparing solver matrix and rhs XYZ vectors =====
	clearSystem();

	SysMatrix = new double*[N];
	sysRhsX = new double[N];
	sysRhsY = new double[N];
	sysRhsZ = new double[N];

	// fill matrix with zeros
	for (uint i = 0; i < N; i++) {
		SysMatrix[i] = new double[N];
		for (uint j = 0; j < N; j++) SysMatrix[i][j] = 0.0;
	}
}

void Evolver::getInterpolatedSDFValuesforVertex(
	Vector3* V, float* SDF_V, Vector3* gradSDF_V, std::vector<Vector3>& positionBuffer, std::vector<float>& valueBuffer)
{
	const uint Nx = sdfGrid->Nx, Ny = sdfGrid->Ny, Nz = sdfGrid->Nz;
	float gradSDFx_V, gradSDFy_V, gradSDFz_V, norm;
	// SDF value
	positionBuffer.clear(); // for old min and max positions
	valueBuffer.clear(); // for SDF cell vertex values
	sdfGrid->getSurroundingCells(*V, Nx, Ny, Nz, sdfGrid->field, &positionBuffer, &valueBuffer);
	*SDF_V = trilinearInterpolate(*V, positionBuffer, valueBuffer);

	// -------- grad SDF --------
	// grad SDF x
	positionBuffer.clear(); // for old min and max positions
	valueBuffer.clear(); // for grad SDF x cell vertex values
	sdfGrid->getSurroundingCells(*V, Nx - 2, Ny - 2, Nz - 2, sdfGrid->gradFieldX, &positionBuffer, &valueBuffer);
	gradSDFx_V = trilinearInterpolate(*V, positionBuffer, valueBuffer);

	// grad SDF y
	positionBuffer.clear(); // for old min and max positions
	valueBuffer.clear(); // for grad SDF y cell vertex values
	sdfGrid->getSurroundingCells(*V, Nx - 2, Ny - 2, Nz - 2, sdfGrid->gradFieldY, &positionBuffer, &valueBuffer);
	gradSDFy_V = trilinearInterpolate(*V, positionBuffer, valueBuffer);

	// grad SDF z
	positionBuffer.clear(); // for old min and max positions
	valueBuffer.clear(); // for grad SDF z cell vertex values
	sdfGrid->getSurroundingCells(*V, Nx - 2, Ny - 2, Nz - 2, sdfGrid->gradFieldZ, &positionBuffer, &valueBuffer);
	gradSDFz_V = trilinearInterpolate(*V, positionBuffer, valueBuffer);

	norm = sqrt(gradSDFx_V * gradSDFx_V + gradSDFy_V * gradSDFy_V + gradSDFz_V * gradSDFz_V);

	if (fabs(norm) > FLT_EPSILON) {
		gradSDF_V->set(gradSDFx_V / norm, gradSDFy_V / norm, gradSDFz_V / norm);
	}
	else {
		gradSDF_V->set(gradSDFx_V, gradSDFy_V, gradSDFz_V);
		// std::cout << "WARNING: ZERO GRADIENT!\n";
	}	
}

float Evolver::tangentialVelocitForVertex(Vector3& V, uint i)
{
	return 0.0f;
}

void Evolver::saveFVAreaScalars()
{
	fvAreas.clear();
	fvAreas = std::vector<float>(N);
	std::vector<std::vector<Vector3>> fvVerts = {};
	std::vector<std::vector<std::vector<uint>>> adjacentPolys = {};
	evolvedSurface->getVertexFiniteVolumes(&fvVerts, &adjacentPolys);

	for (uint i = 0; i < N; i++) {
		uint m = adjacentPolys[i].size();
		// =====================================================================>

		float coVolArea = 0.0f;

		for (uint p = 0; p < m; p++) {
			// vertex
			Vector3 Fi = Vector3(
				evolvedSurface->uniqueVertices[i].x,
				evolvedSurface->uniqueVertices[i].y,
				evolvedSurface->uniqueVertices[i].z
			);

			// midpoints:
			Vector3* M0 = &fvVerts[i][(size_t)2 * p];
			Vector3* baryCenter = &fvVerts[i][(size_t)2 * p + 1];
			Vector3* M1 = &fvVerts[i][((size_t)2 * p + 2) % fvVerts[i].size()];

			// co-vol areas:
			float cV2Area0 = cross(*M0 - Fi, *baryCenter - Fi).length();
			float cV2Area1 = cross(*baryCenter - Fi, *M1 - Fi).length();
			coVolArea += (0.5f * cV2Area0 + 0.5f * cV2Area1);
		}

		fvAreas[i] = coVolArea;
	}

	evolvedSurface->setScalarData(&fvAreas, "CoVolArea");
}

void Evolver::saveInterpolatedSDFValues()
{
	vDistances.clear();
	vDistances = std::vector<float>(N);
	const uint Nx = sdfGrid->Nx, Ny = sdfGrid->Ny, Nz = sdfGrid->Nz;
	for (uint i = 0; i < N; i++) {
		// vertex
		Vector3 Fi = Vector3(
			evolvedSurface->uniqueVertices[i].x,
			evolvedSurface->uniqueVertices[i].y,
			evolvedSurface->uniqueVertices[i].z
		);

		// SDF interpolation from values surrounding V
		std::vector<Vector3> positionBuffer = {};
		std::vector<float> valueBuffer = {};
		float SDF_V;

		// SDF value
		positionBuffer.clear(); // for old min and max positions
		valueBuffer.clear(); // for SDF cell vertex values
		sdfGrid->getSurroundingCells(Fi, Nx, Ny, Nz, sdfGrid->field, &positionBuffer, &valueBuffer);
		SDF_V = trilinearInterpolate(Fi, positionBuffer, valueBuffer);

		vDistances[i] = SDF_V;
	}

	evolvedSurface->setScalarData(&vDistances, "SignedDistance");
}

void Evolver::saveInterpolatedDotValues()
{
	if (vNormals.empty()) vNormals = evolvedSurface->getAngleWeightedVertexPseudoNormals();

	vDotProducts.clear();
	vDotProducts = std::vector<float>(N);
	const uint Nx = sdfGrid->Nx, Ny = sdfGrid->Ny, Nz = sdfGrid->Nz;
	for (uint i = 0; i < N; i++) {
		// vertex
		Vector3 Fi = Vector3(
			evolvedSurface->uniqueVertices[i].x,
			evolvedSurface->uniqueVertices[i].y,
			evolvedSurface->uniqueVertices[i].z
		);

		// SDF interpolation from values surrounding V
		std::vector<Vector3> positionBuffer = {};
		std::vector<float> valueBuffer = {};
		float SDF_V; Vector3 gradSDF_V = Vector3();
		this->getInterpolatedSDFValuesforVertex(&Fi, &SDF_V, &gradSDF_V, positionBuffer, valueBuffer);

		vDotProducts[i] = dot(-1.0f * gradSDF_V, vNormals[i]);
	}

	evolvedSurface->setScalarData(&vDotProducts, "negGradDotN");
}

void Evolver::saveInterpolatedSDFGradients()
{
	vNormals.clear();
	vNormals = evolvedSurface->getAngleWeightedVertexPseudoNormals();

	// clear and allocate the gradient buffer
	vGradients.clear();
	vGradients = std::vector<Vector3>(N);

	for (int i = 0; i < N; i++) {
		// vertex
		Vector3 Fi = Vector3(
			evolvedSurface->uniqueVertices[i].x,
			evolvedSurface->uniqueVertices[i].y,
			evolvedSurface->uniqueVertices[i].z
		);

		// SDF interpolation from values surrounding V
		std::vector<Vector3> positionBuffer = {};
		std::vector<float> valueBuffer = {};
		float SDF_V; Vector3 gradSDF_V = Vector3();
		this->getInterpolatedSDFValuesforVertex(&Fi, &SDF_V, &gradSDF_V, positionBuffer, valueBuffer);

		vGradients[i] = gradSDF_V;
	}
}

void Evolver::saveMeanCurvatureScalars()
{
	evolvedSurface->setScalarData(&vCurvatures, "Mean_Curvature");
}

void Evolver::saveNormalVelocityScalars()
{
	evolvedSurface->setScalarData(&vNormalVelocities, "Normal_Velocity");
}

void Evolver::saveTangentialVelocityScalars()
{
	evolvedSurface->setScalarData(&vTangentialVelocities, "Tangential_Velocity");
}

float Evolver::laplaceBeltramiCtrlFunc(float& SDF_V)
{
	if (epsConstant) return C1;

	return C1 * (1.0f - exp(-(SDF_V * SDF_V) / C2));
}

float Evolver::laplaceBeltramiSmoothFunc(float t)
{
	return initSmoothRate * exp(-smoothDecay * t);
}

float Evolver::etaCtrlFunc(float& SDF_V, Vector3& gradSDF_V, Vector3& nV)
{
	if (etaConstant) return C;

	float gradDotN = dot(-1.0f * gradSDF_V, nV);
	return SDF_V * (gradDotN + D * sqrt(1 - gradDotN * gradDotN)) * C;
}

void Evolver::computeSurfaceNormalsAndCoVolumes()
{
	vNormals.clear();
	vNormals = evolvedSurface->getAngleWeightedVertexPseudoNormals();

	fvVerts.clear(); adjacentPolys.clear();

	fvAreas.clear();
	evolvedSurface->getVertexFiniteVolumes(&fvVerts, &adjacentPolys);
}

//
// [1] Meyer, Desbrun, Schroder, Barr - Discrete Differential-Geometry Operators for Triangulated 2-Manifolds (p. 7)
// [2] Mikula, Remesikova, Sarkoci, Sevcovic - Manifold Evolution With Tangential Redistribution of Points (SIAM Vol 36, p. A1394)
// [3] Tomek, Mikula - DISCRETE DUALITY FINITE VOLUME METHOD WITH TANGENTIAL REDISTRIBUTION OF POINTS FOR SURFACES EVOLVING BY MEAN CURVATURE (p. 1808)
//
void Evolver::getTriangleEvolutionSystem(float smoothStep, float& meanArea)
{
	if (saveDistanceStates) {
		// clear and allocate the distance buffer
		vDistances.clear();
		vDistances = std::vector<float>(N);
	}

	if (saveGradientStates) {
		// clear and allocate the gradient buffer
		vGradients.clear();
		vGradients = std::vector<Vector3>(N);
	}

	meanArea = 0.0f;

	for (uint i = 0; i < N; i++) {
		// <=== part identical to getTriangleEvolutionSystem based on [4] =======
		// vertex
		Vector3 Fi = Vector3(
			evolvedSurface->uniqueVertices[i].x,
			evolvedSurface->uniqueVertices[i].y,
			evolvedSurface->uniqueVertices[i].z
		);

		double eps = 1.0;
		// SDF interpolation from values surrounding V
		std::vector<Vector3> positionBuffer = {};
		std::vector<float> valueBuffer = {};
		float SDF_V; Vector3 gradSDF_V = Vector3();

		if (!meanCurvatureFlow) {
			this->getInterpolatedSDFValuesforVertex(&Fi, &SDF_V, &gradSDF_V, positionBuffer, valueBuffer);
			eps = laplaceBeltramiCtrlFunc(SDF_V);
		}
		else if (smoothing_flag) {
			if (sdfGrid && (saveDistanceStates || saveGradientStates)) {
				this->getInterpolatedSDFValuesforVertex(&Fi, &SDF_V, &gradSDF_V, positionBuffer, valueBuffer);
			}
			eps = laplaceBeltramiSmoothFunc(smoothStep);
		}

		if (saveGradientStates) vGradients[i] = gradSDF_V;
		if (saveDistanceStates) vDistances[i] = SDF_V;

		uint m = adjacentPolys[i].size();
		// =====================================================================>
		
		float coVolArea = 0.0f;

		for (uint p = 0; p < m; p++) {	

			// fill in values from cotan scheme based on [1], [2], [3]:

			// midpoints:
			Vector3* M0 = &fvVerts[i][ (size_t)2 * p ];
			Vector3* baryCenter = &fvVerts[i][ (size_t)2 * p + 1 ];
			Vector3* M1 = &fvVerts[i][ ((size_t)2 * p + 2) % fvVerts[i].size() ];

			// co-vol areas:
			float cV2Area0 = cross(*M0 - Fi, *baryCenter - Fi).length();
			float cV2Area1 = cross(*baryCenter - Fi, *M1 - Fi).length();
			coVolArea += (0.5f * cV2Area0 + 0.5f * cV2Area1);

			// tri vertices:
			Vector3* F1prev = &evolvedSurface->uniqueVertices[adjacentPolys[i][(p > 0 ? p - 1 : m - 1)][1]];
			Vector3* F1 = &evolvedSurface->uniqueVertices[adjacentPolys[i][p][1]];
			Vector3* F2prev = &evolvedSurface->uniqueVertices[adjacentPolys[i][(p > 0 ? p - 1 : m - 1)][2]];
			Vector3* F2 = &evolvedSurface->uniqueVertices[adjacentPolys[i][p][2]];

			// areas:
			float Area1 = cross(*F1prev - Fi, *F1prev - *F2prev).length();
			float Area2 = cross(Fi - *F2, *F1 - *F2).length();

			if (Area1 <= FLT_EPSILON || Area2 <= FLT_EPSILON || coVolArea <= FLT_EPSILON) {
				// count and skip degenerate element occurence
				// std::cout << "WARNING: DEGENERATE ELEMENT at vertex " << i << " !\n";
				degenerateCount++;
				continue;
			}

			// cotans:
			float cotan1 = dot(*F1prev - Fi, *F1prev - *F2prev) / Area1;
			float cotan2 = dot(Fi - *F2, *F1 - *F2) / Area2;
			float cotan1main = dot(Fi - *F1, *F2 - *F1) / Area2;
			float cotan2main = dot(*F1 - *F2, Fi - *F2) / Area2;

			SysMatrix[i][i] += 0.5 * ((double)cotan1main + (double)cotan2main);
			SysMatrix[i][adjacentPolys[i][p][1]] += 0.5 * ((double)cotan1 + (double)cotan2);
		}

		meanArea += coVolArea;

		// diag:
		SysMatrix[i][i] *= ((double)dt * eps) / coVolArea; 
		SysMatrix[i][i] += 1.0;
		// off-diag:
		for (uint p = 0; p < m; p++) SysMatrix[i][adjacentPolys[i][p][1]] *= -((double)dt * eps) / coVolArea;

		// tangential redist:
		float vT = this->tangentialVelocitForVertex(Fi, i);

		// grad dist func:
		float eta = ( meanCurvatureFlow ? 0.0f : this->etaCtrlFunc(SDF_V, gradSDF_V, vNormals[i]) );

		sysRhsX[i] = (double)Fi.x + (double)dt * eta * vNormals[i].x + (double)dt * vT;
		sysRhsY[i] = (double)Fi.y + (double)dt * eta * vNormals[i].y + (double)dt * vT;
		sysRhsZ[i] = (double)Fi.z + (double)dt * eta * vNormals[i].z + (double)dt * vT;

		if (sphereTest || saveAreaStates) fvAreas.push_back(coVolArea); // save areas as weights for numerical tests
	}

	if (saveAreaStates) evolvedSurface->setScalarData(&fvAreas, "CoVolArea");
	if (saveDistanceStates) evolvedSurface->setScalarData(&vDistances, "SignedDistance");

	meanArea /= N;
}

//
// [4] M. Medla - SOLVING PARTIAL DIFFERENTIAL EQUATIONS USING FINITE VOLUME METHOD ON NON-UNIFORM GRIDS (PhD Thesis, p. 29)
//
void Evolver::getQuadEvolutionSystem(float smoothStep, float& meanArea)
{
	for (uint i = 0; i < N; i++) {
		// <=== part identical to getQuadEvolutionSystem based on [1], [2], [3] =======
		// vertex
		Vector3 Fi = Vector3(
			evolvedSurface->uniqueVertices[i].x,
			evolvedSurface->uniqueVertices[i].y,
			evolvedSurface->uniqueVertices[i].z
		);

		// SDF interpolation from values surrounding V
		std::vector<Vector3> positionBuffer = {};
		std::vector<float> valueBuffer = {};
		float SDF_V, gradSDFx_V, gradSDFy_V, gradSDFz_V; Vector3 gradSDF_V = Vector3();

		this->getInterpolatedSDFValuesforVertex(&Fi, &SDF_V, &gradSDF_V, positionBuffer, valueBuffer);

		float eps = laplaceBeltramiCtrlFunc(SDF_V);

		uint m = adjacentPolys[i].size();
		// ===========================================================================>

		float coVolArea = 0.0f;

		for (uint p = 0; p < m; p++) {

			// fill in values from bilinear scheme based on [4]:


		}

		// diag:
		SysMatrix[i][i] *= -(dt * eps) / (2.0 * coVolArea); SysMatrix[i][i] += 1.0;

		// tangential redist:
		float vT = this->tangentialVelocitForVertex(Fi, i);

		// grad dist func:
		float eta = (meanCurvatureFlow ? 0.0f : this->etaCtrlFunc(SDF_V, gradSDF_V, vNormals[i]));

		sysRhsX[i] += (double)Fi.x + (double)dt * eta * vNormals[i].x + (double)dt * vT;
		sysRhsY[i] += (double)Fi.y + (double)dt * eta * vNormals[i].y + (double)dt * vT;
		sysRhsZ[i] += (double)Fi.z + (double)dt * eta * vNormals[i].z + (double)dt * vT;
	}
}

void Evolver::getCurvaturesAndNormalVelocities()
{
	vCurvatures.clear(); vNormalVelocities.clear();
	vCurvatures = std::vector<float>(N);
	vNormalVelocities = std::vector<float>(N);

	for (int i = 0; i < N; i++) {
		// vertex
		Vector3 Fi = Vector3(
			evolvedSurface->uniqueVertices[i].x,
			evolvedSurface->uniqueVertices[i].y,
			evolvedSurface->uniqueVertices[i].z
		);

		double eps = 1.0;
		// SDF interpolation from values surrounding V
		std::vector<Vector3> positionBuffer = {};
		std::vector<float> valueBuffer = {};
		float SDF_V; Vector3 gradSDF_V = Vector3();

		if (!meanCurvatureFlow) {
			this->getInterpolatedSDFValuesforVertex(&Fi, &SDF_V, &gradSDF_V, positionBuffer, valueBuffer);
			eps = laplaceBeltramiCtrlFunc(SDF_V);
		}

		uint m = adjacentPolys[i].size();
		// ===========================================================================>

		float coVolArea = 0.0f;
		Vector3 vCurvatureVect = Vector3();

		for (uint p = 0; p < m; p++) {
			// midpoints:
			Vector3* M0 = &fvVerts[i][(size_t)2 * p];
			Vector3* baryCenter = &fvVerts[i][(size_t)2 * p + 1];
			Vector3* M1 = &fvVerts[i][((size_t)2 * p + 2) % fvVerts[i].size()];

			// co-vol areas:
			float cV2Area0 = cross(*M0 - Fi, *baryCenter - Fi).length();
			float cV2Area1 = cross(*baryCenter - Fi, *M1 - Fi).length();
			coVolArea += (0.5f * cV2Area0 + 0.5f * cV2Area1);

			// tri vertices:
			Vector3* F1prev = &evolvedSurface->uniqueVertices[adjacentPolys[i][(p > 0 ? p - 1 : m - 1)][1]];
			Vector3* F1 = &evolvedSurface->uniqueVertices[adjacentPolys[i][p][1]];
			Vector3* F2prev = &evolvedSurface->uniqueVertices[adjacentPolys[i][(p > 0 ? p - 1 : m - 1)][2]];
			Vector3* F2 = &evolvedSurface->uniqueVertices[adjacentPolys[i][p][2]];

			// areas:
			float Area1 = cross(*F1prev - Fi, *F1prev - *F2prev).length();
			float Area2 = cross(Fi - *F2, *F1 - *F2).length();

			if (Area2 <= FLT_EPSILON || Area1 <= FLT_EPSILON || coVolArea <= FLT_EPSILON) {
				// skip degenerate element
				continue;
			}

			// cotans:
			float cotan1 = dot(*F1prev - Fi, *F1prev - *F2prev) / Area1;
			float cotan2 = dot(Fi - *F2, *F1 - *F2) / Area2;

			vCurvatureVect += 0.5 * ((double)cotan1 + (double)cotan2) * (Fi - *F1);
		}

		vCurvatures[i] = ((1.0f / coVolArea) * vCurvatureVect).length();

		// grad dist func:
		float eta = (meanCurvatureFlow ? 0.0f : this->etaCtrlFunc(SDF_V, gradSDF_V, vNormals[i]));

		// normal velocity magnitude (will scale surface unit normal):
		vNormalVelocities[i] = eps * vCurvatures[i] + eta;
	}
}

void Evolver::openLogs()
{
	if (writeGenericLog) {
		log = std::fstream(geomName + "_evolutionLog.txt", std::fstream::out);
		log << log_header;
	}

	if (writeMeanAreaLog) {
		meanAreaLog = std::fstream(geomName + "_meanAreaLog.txt", std::fstream::out);
	}

	if (writeErrorLog) {
		errorLog = std::fstream(geomName + "_errorLog(" + std::to_string(testId) + ").txt", std::fstream::out);
		errorLog << log_header;
	}

	if (writeTimeLog) {
		timingLog = std::fstream(geomName + "_timingLog.txt", std::fstream::out);
		timingLog << log_header;
	}
}

void Evolver::closeLogs()
{
	if (writeGenericLog) log.close();
	if (writeErrorLog) errorLog.close();
	if (writeMeanAreaLog) meanAreaLog.close();
	if (writeTimeLog) timingLog.close();
}

void Evolver::updateGeometry(double* Fx, double* Fy, double* Fz)
{
	evolvedSurface->normals.clear();
	for (uint i = 0; i < N; i++) {
		evolvedSurface->uniqueVertices[i].set(Fx[i], Fy[i], Fz[i]);
		// interrupt if solution explodes
		this->interruptEvolution = !this->sdfGrid->bbox.isInside(evolvedSurface->uniqueVertices[i]);
	}
	evolvedSurface->fillVerticesFromUniqueVertices();
}

void Evolver::exportGeometry(int step)
{
	VTKExporter e = VTKExporter();
	e.initExport(evolvedSurface, geomName + (testId == -1 ? ""  : ("(" + std::to_string(testId) + ")")) + "_" + std::to_string(step));
}

void Evolver::exportResultGeometry(std::string suffix)
{
	VTKExporter e = VTKExporter();
	e.initExport(evolvedSurface, geomName + suffix);
}

void Evolver::exportVectorStates(int step)
{
	VTKExporter e = VTKExporter();
	e.exportVectorDataOnGeometry(evolvedSurface, &vNormals, geomName + "_vNormals_" + std::to_string(step));
	e.exportVectorDataOnGeometry(evolvedSurface, &vGradients, geomName + "_vGrads_" + std::to_string(step));
}

void Evolver::exportTestGeometry(int step, float t)
{
	VTKExporter e = VTKExporter();
	float r = sqrt(r0 * r0 - 4 * t);
	Geometry sg;
	if (elType == ElementType::tri) {
		uint n = 3;
		sg = IcoSphere(n, r);
	}
	else {
		uint n = 5;
		sg = CubeSphere(n, r);
	}
	e.initExport(&sg, "exactSphere_" + (testId == -1 ? "" : ("(" + std::to_string(testId) + ")")) + "_" + std::to_string(step));
}

Evolver::Evolver()
{
}

// ========== sphere test constructor ========
//
Evolver::Evolver(EvolutionParams& eParams, SphereTestParams& stParams)
{
	this->sphereTest = true;
	this->meanCurvatureFlow = true;
	this->dt = eParams.dt; 
	this->tStop = eParams.tStop;
	this->NSteps = std::floor(eParams.tStop / eParams.dt);
	this->subdiv = eParams.subdiv;
	this->elType = eParams.elType;

	this->geomName = eParams.name;
	// console flags
	this->printHappenings = eParams.printHappenings;
	this->printStepOutput = eParams.printStepOutput;
	this->printSolution = eParams.printSolution;
	// log flags
	this->writeGenericLog = eParams.writeGenericLog;
	this->writeErrorLog = stParams.writeErrorLog;
	this->writeTimeLog = eParams.writeTimeLog;
	// save flags
	this->saveStates = eParams.saveStates;
	this->saveResult = eParams.saveResult;
	
	// test params
	this->testId = stParams.testId;
	this->r0 = stParams.r0;	

	if (eParams.sourceGeometry) this->evolvedSurface = eParams.sourceGeometry;
	// -----------------------------------------------------

	this->init();

	this->openLogs();
	this->evolve();
	this->closeLogs();
}

// ========= applied evolve constructor ============
//
Evolver::Evolver(EvolutionParams& eParams, MeanCurvatureParams& mcfParams, GradDistanceParams& sdfParams, TangentialRedistParams* tanParams)
{
	this->sphereTest = false;
	this->meanCurvatureFlow = (sdfParams.sdfGrid == nullptr); // no sdf => only MCF
	this->dt = eParams.dt;
	this->tStop = eParams.NSteps * eParams.dt;
	this->subdiv = eParams.subdiv;
	this->NSteps = eParams.NSteps;
	this->elType = eParams.elType;

	this->geomName = eParams.name;
	// console flags
	this->printHappenings = eParams.printHappenings;
	this->printStepOutput = eParams.printStepOutput;
	this->printSolution = eParams.printSolution;
	// log flags
	this->writeGenericLog = eParams.writeGenericLog;
	this->writeMeanAreaLog = mcfParams.writeMeanAreaLog;
	this->writeTimeLog = eParams.writeTimeLog;
	// save flags
	this->saveStates = eParams.saveStates;
	this->saveResult = eParams.saveResult;
	this->saveAreaStates = mcfParams.saveAreaStates;
	this->saveDistanceStates = sdfParams.saveDistanceStates;
	this->saveGradientStates = sdfParams.saveGradientStates;
	this->saveCurvatureStates = mcfParams.saveCurvatureStates;
	this->saveNormalVelocityStates = mcfParams.saveNormalVelocityStates;

	// evolution type	
	if (eParams.sourceGeometry) this->evolvedSurface = eParams.sourceGeometry;
	this->sdfGrid = sdfParams.sdfGrid;
	this->targetGeom = sdfParams.targetGeom;

	// scalar ctrl params
	this->C1 = mcfParams.C1;
	this->C2 = mcfParams.C2;
	this->C = sdfParams.C;
	this->D = sdfParams.D;

	this->epsConstant = mcfParams.constant;
	this->etaConstant = sdfParams.constant;

	// additional MCF smoothing
	this->smoothSteps = mcfParams.smoothSteps;
	this->initSmoothRate = mcfParams.initSmoothRate;
	this->smoothDecay = mcfParams.smoothDecay;

	// tangential redistribution
	this->tangential_redist = tanParams != nullptr;
	if (tangential_redist) {
		this->omega = tanParams->omega;
		this->saveTangentialVelocityStates = tanParams->saveTangentialVelocityStates;
	}
	// -----------------------------------------------------

	this->init();

	this->openLogs();
	this->evolve();

	// additional MCF smooting
	if (smoothSteps > 0) {
		this->tBegin = this->NSteps + 1;
		this->meanCurvatureFlow = true;
		this->smoothing_flag = true;
		this->evolve();
	}
	this->closeLogs();
}

Evolver::~Evolver()
{
	this->clearSystem();
	this->fvAreas.clear();
	this->vCurvatures.clear();
	this->vNormals.clear();
	this->vDistances.clear();
	this->vGradients.clear();
	this->vDotProducts.clear();

	delete solveX; delete solveY; delete solveZ;
}


float Evolver::getSphereStepError(float t)
{
	// r(t) = sqrt(r0 * r0 - 4 * t);
	// error for step t is the norm of the difference (Fi - center) - r(t), weighted by the vertex co-volume area
	// summed through all vertices, of course

	float FRadius, rt = sqrt(r0 * r0 - 4 * t);
	float result = 0.0f;

	for (uint i = 0; i < N; i++) {
		FRadius = (evolvedSurface->uniqueVertices[i] - center).length();
		result += fabs(FRadius - rt) * fvAreas[i];
	}
	return result;
}

float Evolver::getSphereStepL2Error(float t)
{
	// r(t) = sqrt(r0 * r0 - 4 * t);
	// error for step t is the square of the difference ((Fi - center) - r(t)), weighted by the vertex co-volume area
	// summed through all vertices, of course

	float FRadius, rt = sqrt(r0 * r0 - 4 * t);
	float result = 0.0f;

	for (uint i = 0; i < N; i++) {
		FRadius = (evolvedSurface->uniqueVertices[i] - center).length();
		result += (FRadius - rt) * (FRadius - rt) * fvAreas[i];
	}
	return result;
}