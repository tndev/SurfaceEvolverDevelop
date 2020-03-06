#include "SurfaceEvolutionSolver.h"

void SurfaceEvolutionSolver::clearSystem()
{
	if (SysMatrix != nullptr) {
		for (uint i = 0; i < N; i++) delete[] SysMatrix[i];
		delete[] SysMatrix;
	}
	if (sysRhsX != nullptr) delete[] sysRhsX;
	if (sysRhsY != nullptr) delete[] sysRhsY;
	if (sysRhsY != nullptr) delete[] sysRhsY;

	SysMatrix = nullptr;
	sysRhsX = nullptr;
	sysRhsY = nullptr;
	sysRhsZ = nullptr;
}

void SurfaceEvolutionSolver::initSystem()
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

void SurfaceEvolutionSolver::getInterpolatedSDFValuesforVertex(
	Vector3* V, float* SDF_V, Vector3* gradSDF_V, std::vector<Vector3>& positionBuffer, std::vector<float>& valueBuffer)
{
	const uint Nx = sdfGrid->Nx, Ny = sdfGrid->Ny, Nz = sdfGrid->Nz;
	float gradSDFx_V, gradSDFy_V, gradSDFz_V;
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

	gradSDF_V->set(gradSDFx_V, gradSDFy_V, gradSDFz_V);
}

float SurfaceEvolutionSolver::tangentialVelocitForVertex(Vector3& V, uint i)
{
	return 0.0f;
}

float SurfaceEvolutionSolver::laplaceBeltramiCtrlFunc(float& SDF_V)
{
	// C1 = 1.0f; C2 = 1.0f;
	// return C1 * (1.0f - exp(SDF_V * SDF_V / C2);
	return 1.0f;
}

float SurfaceEvolutionSolver::etaCtrlFunc(float& SDF_V, Vector3& gradSDF_V, Vector3& nV)
{
	float C = 0.0f;
	float gradDotN = dot(gradSDF_V, nV);
	return SDF_V * (fabs(gradDotN) + sqrt(1 - gradDotN * gradDotN)) * C;
}

//
// [1] Meyer, Desbrun, Schroder, Barr - Discrete Differential-Geometry Operators for Triangulated 2-Manifolds (p. 7)
// [2] Mikula, Remesikova, Sarkoci, Sevcovic - Manifold Evolution With Tangential Redistribution of Points (SIAM Vol 36, p. A1394)
// [3] Tomek, Mikula - DISCRETE DUALITY FINITE VOLUME METHOD WITH TANGENTIAL REDISTRIBUTION OF POINTS FOR SURFACES EVOLVING BY MEAN CURVATURE (p. 1808)
//
void SurfaceEvolutionSolver::getTriangleEvolutionSystem(
	std::vector<Vector3>& vNormals, std::vector<std::vector<Vector3>>& fvVerts, std::vector<std::vector<std::vector<uint>>>& adjacentPolys)
{
	for (uint i = 0; i < N; i++) {
		// <=== part identical to getTriangleEvolutionSystem based on [4] =======
		// vertex
		Vector3 Fi = Vector3(
			evolvedSurface->uniqueVertices[i].x,
			evolvedSurface->uniqueVertices[i].y,
			evolvedSurface->uniqueVertices[i].z
		);

		float eps = 1.0f;
		// SDF interpolation from values surrounding V
		std::vector<Vector3> positionBuffer = {};
		std::vector<float> valueBuffer = {};
		float SDF_V, gradSDFx_V, gradSDFy_V, gradSDFz_V; Vector3 gradSDF_V = Vector3();

		if (!meanCurvatureFlow) {
			this->getInterpolatedSDFValuesforVertex(&Fi, &SDF_V, &gradSDF_V, positionBuffer, valueBuffer);
			eps = laplaceBeltramiCtrlFunc(SDF_V);
		}		

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
			Vector3* F1prev = &evolvedSurface->uniqueVertices[adjacentPolys[i][(p - 1) % m][1]];
			Vector3* F1 = &evolvedSurface->uniqueVertices[adjacentPolys[i][p][1]]; Vector3* F2prev = F1;
			Vector3* F2 = &evolvedSurface->uniqueVertices[adjacentPolys[i][p][2]];

			// areas:
			float Area1 = cross(*F1prev - Fi, *F1prev - *F2prev).length();
			float Area2 = cross(Fi - *F2, *F1 - *F2).length();

			// cotans:
			float cotan1 = dot(*F1prev - Fi, *F1prev - *F2prev) / Area1;
			float cotan2 = dot(Fi - *F2, *F1 - *F2) / Area2;

			SysMatrix[i][i] -= 0.5 * ((double)cotan1 + (double)cotan2);
			SysMatrix[i][adjacentPolys[i][p][i]] += ((double)cotan1 + (double)cotan2) * (dt * eps) / (4.0 * coVolArea);
		}

		// diag:
		SysMatrix[i][i] *= -(dt * eps) / (2.0 * coVolArea); SysMatrix[i][i] += 1.0;

		// tangential redist:
		float vT = this->tangentialVelocitForVertex(Fi, i);

		// grad dist func:
		float eta = ( meanCurvatureFlow ? 0.0f : this->etaCtrlFunc(SDF_V, gradSDF_V, vNormals[i]) );

		sysRhsX[i] += (double)Fi.x + (double)dt * eta + (double)dt * vT;

		if (performSphereTest) fvAreas.push_back(coVolArea); // save areas as weights for numerical tests
	}
}

//
// [4] M. Medla - SOLVING PARTIAL DIFFERENTIAL EQUATIONS USING FINITE VOLUME METHOD ON NON-UNIFORM GRIDS (PhD Thesis, p. 29)
//
void SurfaceEvolutionSolver::getQuadEvolutionSystem(
	std::vector<Vector3>& vNormals, std::vector<std::vector<Vector3>>& fvVerts, std::vector<std::vector<std::vector<uint>>>& adjacentPolys)
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
		float eta = this->etaCtrlFunc(SDF_V, gradSDF_V, vNormals[i]);

		sysRhsX[i] += (double)Fi.x + (double)eta + (double)dt * vT;
	}
}

// ======= Solver helpers ============

void SurfaceEvolutionSolver::printArray1(std::string name, double* a, int printLim, bool inRow)
{
	if (inRow) {
		std::cout << name << " = " << std::endl;
		for (int i = 0; i < printLim; i++) {
			std::cout << " " << a[i];
		}
		std::cout << "  ... ";
		for (int i = N - printLim; i < N; i++) {
			std::cout << " " << a[i];
		}
		std::cout << std::endl;
	}
	else {
		std::string offset = std::string((name + " = ").length() + 1, ' ');
		std::cout << name << " = " << std::endl;
		for (int i = 0; i < printLim; i++) {
			std::cout << offset << a[i] << std::endl;
		}
		for (int i = 0; i < 3; i++) std::cout << offset << "  ." << std::endl;
		for (int i = N - printLim; i < N; i++) {
			std::cout << offset << a[i] << std::endl;
		}
		std::cout << std::endl;
	}
}

void SurfaceEvolutionSolver::printArray2(std::string name, double** A, int printLim)
{
	std::string offset = std::string((name + " = ").length() + 1, ' ');
	std::cout << name << " = " << std::endl;
	for (int i = 0; i < printLim; i++) {
		std::cout << offset;
		for (int j = 0; j < printLim; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << "  ...  ";
		for (int j = N - printLim; j < N; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << std::endl;
	}
	for (int i = 0; i < 3; i++) std::cout << offset << "  ." << std::endl;
	for (int i = N - printLim; i < N; i++) {
		std::cout << offset;
		for (int j = 0; j < printLim; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << "  ...  ";
		for (int j = N - printLim; j < N; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

double SurfaceEvolutionSolver::vectorDot(double* a, double* b)
{
	double result = 0.;
	for (int i = 0; i < N; i++) result += a[i] * b[i];
	return result;
}

double SurfaceEvolutionSolver::vectorNorm(double* a)
{
	return sqrt(vectorDot(a, a));
}

// ========== SOLVER ======================

void SurfaceEvolutionSolver::Bi_CGSTAB_solve(double** A, double* b, double* x, bool print)
{
	// ctrl. constants
	int maxIter = 100;
	double tol = 1e-6;

	// iter vectors
	double* x_curr = new double[N];
	double* x_next = new double[N];

	double* r_curr = new double[N];
	double* r_next = new double[N];

	double* rp0 = new double[N];

	double* p_curr = new double[N];
	double* p_next = new double[N];

	double* s = new double[N];

	double* tmp = new double[N];
	double* tmp1 = new double[N];

	// iter scalars
	double omega, alpha, beta, norm;

	// x0 = (1000,1000,...,1000)
#pragma omp parallel for
	for (int i = 0; i < N; i++) x_curr[i] = 1000.;
	// r0 = b - A x0
	// choose rp0 such that <r0, rp0> != 0
	// p0 = r0
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		r_curr[i] = b[i];
		for (int j = 0; j < N; j++) {
			r_curr[i] -= A[i][j] * x_curr[j];
		}
		rp0[i] = r_curr[i] + 100;
		p_curr[i] = r_curr[i];
	}
	if (print) {
		std::cout << "==================================================" << std::endl;
		std::cout << "----------- Initializing Bi-CGSTAB Method --------" << std::endl;
		printArray2("systemMatrix", A, 4);
		printArray1("systemRhs", b, 5);
		printArray1("x0", x_curr, 2);
		printArray1("r0", r_curr, 5);
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << "------------ Launching iterations ----------------" << std::endl;
	}

	// begin iterations
	for (int k = 0; k < maxIter; k++) {
#pragma omp parallel
#pragma omp single
		if (print) std::cout << "::: iter : " << k << std::endl;

		// alpha[k] = <r[k], rp0> / <Ap[k], rp0>

		double num = 0.; double den = 0.;
#pragma omp parallel for reduction (+:num, den)
		for (int i = 0; i < N; i++) {
			tmp[i] = 0.;
			for (int j = 0; j < N; j++) {
				tmp[i] += A[i][j] * p_curr[j];
			}
			num += r_curr[i] * rp0[i];
			den += tmp[i] * rp0[i];
		}
		alpha = num / den;

		// s[k] = r[k] - alpha[k] * A p[k]
#pragma omp parallel for 
		for (int i = 0; i < N; i++) {
			s[i] = r_curr[i] - alpha * tmp[i];
		}

		norm = vectorNorm(s);
		if (print) std::cout << "||s|| = " << norm << std::endl;
		if (norm < tol) {
			// x[k + 1] = x[k] + alpha[k] * p[k]
#pragma omp parallel for
			for (int i = 0; i < N; i++) {
				x_next[i] = x_curr[i] + alpha * p_curr[i];
			}

			if (print) std::cout << "||s|| < tol = " << tol << ", exiting iterations" << std::endl;
			break;
		}

		// omega[k] = <A s[k], s[k]> / <A s[k], A s[k]>

		num = 0; den = 0;
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			tmp[i] = 0;
			for (int j = 0; j < N; j++) {
				tmp[i] += A[i][j] * s[j];
			}
			num += tmp[i] * s[i];
			den += tmp[i] * tmp[i];
		}
		omega = num / den;

		// x[k + 1] = x[k] + alpha[k] * p[k] + omega[k] * s[k]
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			x_next[i] = x_curr[i] + alpha * p_curr[i] + omega * s[i];
		}

		// r[k + 1] = s[k] - omega[k] * A s[k]
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			r_next[i] = s[i] - omega * tmp[i];
		}

		norm = vectorNorm(r_next);
#pragma omp parallel
#pragma omp single
		if (print) std::cout << "||r[k + 1]|| = " << norm << std::endl;
		if (norm < tol) {
#pragma omp parallel
#pragma omp single
			if (print) std::cout << "||r[k + 1]|| < tol = " << tol << ", exiting iterations" << std::endl;
			break;
		}

		// beta[k] = (alpha[k] / omega[k]) * <r[k + 1], rp0> / <r[k], rp0>

		num = 0; den = 0;
#pragma omp parallel for reduction(+: num, den)
		for (int i = 0; i < N; i++) {
			num += r_next[i] * rp0[i];
			den += r_curr[i] * rp0[i];
		}

		beta = (alpha / omega) * num / den;

		// p[k + 1] = r[k + 1] + beta[k] * (p[k] - omega[k] * A p[k])
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			tmp[i] = 0;
			for (int j = 0; j < N; j++) {
				tmp[i] += A[i][j] * p_curr[j];
			}
			p_next[i] = r_next[i] + beta * (p_curr[i] - omega * tmp[i]);
		}

		norm = fabs(vectorDot(r_next, rp0));
#pragma omp parallel
#pragma omp single
		if (print) std::cout << "|< r[k + 1], rp0 >| = " << norm << std::endl;
		if (norm < tol) {
			// rp0 = r[k + 1]; p[k + 1] = r[k + 1]
#pragma omp parallel for
			for (int i = 0; i < N; i++) {
				rp0[i] = r_next[i]; p_next[i] = r_next[i];
			}
		}
		// current = next
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			x_curr[i] = x_next[i];
			r_curr[i] = r_next[i];
			p_curr[i] = p_next[i];
		}
#pragma omp parallel
#pragma omp single
		if (print) std::cout << "===> finishing iter " << k << std::endl;
	}

	// result: x = x_next
#pragma omp parallel for
	for (int i = 0; i < N; i++) x[i] = x_next[i];

	// clean up
	delete[] x_curr; delete[] x_next;
	delete[] r_curr; delete[] r_next;
	delete[] p_curr; delete[] p_next;
	delete[] s; delete[] tmp; delete[] tmp1;
}

void SurfaceEvolutionSolver::updateGeometry(double* Fx, double* Fy, double* Fz)
{
	evolvedSurface->normals.clear();
	for (uint i = 0; i < N; i++) {
		evolvedSurface->uniqueVertices[i].set(Fx[i], Fy[i], Fz[i]);
	}
	evolvedSurface->fillVerticesFromUniqueVertices();
}

// ====== end solver stuff ================================

SurfaceEvolutionSolver::SurfaceEvolutionSolver()
{
}

SurfaceEvolutionSolver::SurfaceEvolutionSolver(float dt, float tStop, int NSteps, ElementType type, Grid* sdfGrid, std::string name, bool sphereTest, bool saveStates)
{
	// ===== Evolution params ==============================
	this->saveStates = saveStates; 
	this->performSphereTest = sphereTest; 
	if (sphereTest || sdfGrid == nullptr) this->meanCurvatureFlow = true;
	this->geomName = name;
	
	if (NSteps < 0) {
		this->NSteps = std::floor(tStop / dt);
	}
	else {
		this->NSteps = NSteps; 
	}
	
	this->dt = dt; this->tStop = tStop;
	this->sdfGrid = sdfGrid;
	this->type = type;
	// -----------------------------------------------------

	this->init();
	this->evolve();
}

SurfaceEvolutionSolver::~SurfaceEvolutionSolver()
{
	this->clearSystem(); this->time_logs.clear();
	this->fvAreas.clear();
}

void SurfaceEvolutionSolver::init()
{
	// ===== Preparing a source geometry for immersion =====
	// -----------------------------------------------------
	// radius
	float r = (meanCurvatureFlow ? r0 : 0.499f * std::min({ sdfGrid->scale.x, sdfGrid->scale.y, sdfGrid->scale.z }));
	if (type == ElementType::tri) {
		uint n = 2;
		evolvedSurface = new IcoSphere(n, r);
	}
	else {
		uint n = 3;
		evolvedSurface = new CubeSphere(n, r);
	}
	Box3 bbox = (meanCurvatureFlow ? evolvedSurface->getBoundingBox() : sdfGrid->bbox);
	this->center = bbox.getCenter();
	Matrix4 M = Matrix4().makeTranslation(center.x, center.y, center.z);
	evolvedSurface->applyMatrix(M);
	this->N = evolvedSurface->uniqueVertices.size();

	std::cout << "============================================================= " << std::endl;
	std::cout << "----------------------- Initializing ------------------------" << std::endl;
	std::cout << "------------ S U R F A C E    E V O L U T I O N -------------" << std::endl;
	std::cout << ">> NSteps = " << NSteps << ", dt = " << dt << ", tStop = " << tStop << std::endl;
	std::cout << ">> Model: " << geomName << ", NVerts = " << N << std::endl;
	if (!meanCurvatureFlow) {
		std::cout << ">> SDF grid bounds: min = " << sdfGrid->bbox.min << ", max = " << sdfGrid->bbox.max << std::endl;
		std::cout << ">> SDF resolution: " << sdfGrid->Nx << " x " << sdfGrid->Ny << " x " << sdfGrid->Nz << std::endl;
	}
	else {
		std::cout << ">> Performing mean curvature flow test on sphere with r0 = " << r0 << std::endl;
	}

	std::cout << "-------------------------------------------------------------" << std::endl;
	std::cout << "starting evolution ..." << std::endl;

	// -------------------------------------------------------
}

void SurfaceEvolutionSolver::evolve()
{
	if (!meanCurvatureFlow) {
		if (!sdfGrid->hasGradient()) {
			sdfGrid->computeGradient();
		}
	}

	float sphereTestL2Error = 0.0f;

	for (int ti = 0; ti < NSteps; ti++) {
		
		float t = (float)ti / (float)NSteps * tStop;
		std::cout << "-------- time step ti = " << ti << ", t = " << t << " -----------" << std::endl;

		// <===== P S E U D O N O R M A L S ========
		auto startPseudoNormals = std::chrono::high_resolution_clock::now();

		if (printHappenings) std::cout << "computing surface pseudonormals ..." << std::endl;
		std::vector<Vector3> vNormals = evolvedSurface->getAngleWeightedVertexPseudoNormals();
		if (printHappenings) std::cout << "... pseudonormals computed" << std::endl;

		auto endPseudoNormals = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedPseudoNormals = (endPseudoNormals - startPseudoNormals);
		// ========================================>
		


		fvAreas.clear();
		std::vector<std::vector<Vector3>> fvVerts = {};
		std::vector<std::vector<std::vector<uint>>> adjacentPolys = {};

		// <===== F I N I T E   V O L U M E S ========
		auto startFiniteVolumes = std::chrono::high_resolution_clock::now();

		if (printHappenings) std::cout << "computing surf. finite volumes ..." << std::endl;
		evolvedSurface->getVertexFiniteVolumes(&fvVerts, &adjacentPolys);
		if (printHappenings) std::cout << "... finite vols computed" << std::endl;

		auto endFiniteVolumes = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedFiniteVolumes = (endFiniteVolumes - startFiniteVolumes);
		// ========================================>


		
		this->initSystem(); // prep linear system of dim NVerts * NVerts



		// <===== M A T R I X    F I L L ==========
		auto startFillMatrix = std::chrono::high_resolution_clock::now();

		if (printHappenings) std::cout << "filling sys matrix & rhs ..." << std::endl;
		if (this->type == ElementType::tri) {
			this->getTriangleEvolutionSystem(vNormals, fvVerts, adjacentPolys);
		} else {
			this->getQuadEvolutionSystem(vNormals, fvVerts, adjacentPolys);
		}
		if (printHappenings) std::cout << "... matrix & rhs filled" << std::endl;

		auto endFillMatrix = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedFillMatrix = (endFillMatrix - startFillMatrix);
		// ========================================>



		double* FnextX = new double[N];
		double* FnextY = new double[N];
		double* FnextZ = new double[N];



		// <===== L I N   S Y S   S O L V E ========
		auto startLinSolve = std::chrono::high_resolution_clock::now();

		if (printHappenings) std::cout << "3 x " << N << " x " << N << "-sys Bi CGSTAB solve ..." << std::endl;
		Bi_CGSTAB_solve(SysMatrix, sysRhsX, FnextX);
		Bi_CGSTAB_solve(SysMatrix, sysRhsY, FnextY);
		Bi_CGSTAB_solve(SysMatrix, sysRhsZ, FnextZ);
		if (printHappenings) std::cout << "... done" << std::endl;

		auto endLinSolve = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedLinSolve = (endLinSolve - startLinSolve);
		// ========================================>



		// <====== G E O M   U P D A T E ==========
		auto startUpdate = std::chrono::high_resolution_clock::now();

		if (printHappenings) std::cout << "updating surface geometry ..." << std::endl;
		updateGeometry(FnextX, FnextY, FnextZ);
		if (printHappenings) std::cout << "... done" << std::endl << std::endl;


		if (performSphereTest) sphereTestL2Error += this->getSphereStepL2Error(t);


		auto endUpdate = std::chrono::high_resolution_clock::now();
		std::chrono::duration<float> elapsedUpdate = (endUpdate - startUpdate);
		// ========================================>

		std::string time_log =
			"step " + std::to_string(t) + "\n" +
			"surface pseudonormals: " + std::to_string(elapsedPseudoNormals.count()) + " s, surface finite volumes: " + std::to_string(elapsedPseudoNormals.count()) + " s\n" +
			"filling sys matrix: " + std::to_string(elapsedFillMatrix.count()) + " s\n" +
			"(vector) linear sys solve: " + std::to_string(elapsedLinSolve.count()) + " s\n" +
			"geom update: " + std::to_string(elapsedUpdate.count()) + " s\n" +
			" evolution step total: " + std::to_string(elapsedPseudoNormals.count() + elapsedPseudoNormals.count() + elapsedFillMatrix.count() + elapsedLinSolve.count() + elapsedUpdate.count()) + " s\n";
		if (printStepOutput) {
			std::cout << time_log << std::endl;
		}
	}

	if (performSphereTest) {
		sphereTestL2Error *= dt;
	}
}

float SurfaceEvolutionSolver::getSphereStepL2Error(float t)
{
	// r(t) = sqrt(r0 * r0 - 4 * t);
	// error for step t is the norm of the difference (Fi - center) - r(t), weighted by the vertex co-volume area
	// summed through all vertices, of course

	float FRadius;
	float result = 0.0f;

	for (uint i = 0; i < N; i++) {
		FRadius = (evolvedSurface->uniqueVertices[i] - center).length();
		result += fabs(FRadius - sqrt(r0 * r0 - 4 * t)) * fvAreas[i];
	}
	return result;
}
