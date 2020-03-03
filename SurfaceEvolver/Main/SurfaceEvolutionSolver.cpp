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
		for (uint j = 0; j < N; j++) SysMatrix[i][j] = 0.0f;
	}
}


void SurfaceEvolutionSolver::getTriangleVertexLaplaceBeltrami(uint i, std::vector<std::vector<uint>>* adjacentPolys)
{
}

void SurfaceEvolutionSolver::getQuadVertexLaplaceBeltrami(uint i, std::vector<std::vector<uint>>* adjacentPolys)
{
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

// ====== end solver stuff ================================

SurfaceEvolutionSolver::SurfaceEvolutionSolver()
{
}

SurfaceEvolutionSolver::SurfaceEvolutionSolver(const SurfaceEvolutionSolver& other)
{
}

SurfaceEvolutionSolver::SurfaceEvolutionSolver(Grid* sdfGrid, uint NSteps, float dt, ElementType type, bool saveStates)
{
	// ===== Evolution params ==============================
	this->NSteps = NSteps; this->saveStates = saveStates;

	this->dt = dt;
	this->sdfGrid = sdfGrid;
	this->type = type;
	// -----------------------------------------------------

	this->init();

	this->evolve();
}

SurfaceEvolutionSolver::~SurfaceEvolutionSolver()
{
}

void SurfaceEvolutionSolver::init()
{
	// ===== Preparing a source geometry for immersion =====
	// -----------------------------------------------------
	// radius
	float r = 0.499f * std::min({ sdfGrid->scale.x, sdfGrid->scale.y, sdfGrid->scale.z });
	if (type == ElementType::tri) {
		uint n = 2;
		evolvedSurface = new IcoSphere(n, r);
	}
	else {
		uint n = 3;
		evolvedSurface = new CubeSphere(n, r);
	}
	Vector3 center = sdfGrid->bbox.getCenter();
	Matrix4 M = Matrix4().makeTranslation(center.x, center.y, center.z);
	evolvedSurface->applyMatrix(M);
	this->N = evolvedSurface->uniqueVertices.size();

	// -------------------------------------------------------
}

void SurfaceEvolutionSolver::evolve()
{
	const uint Nx = sdfGrid->Nx, Ny = sdfGrid->Ny, Nz = sdfGrid->Nz;
	Vector3 V = Vector3();
	// for SDF interpolation from values surrounding V
	std::vector<Vector3> positionBuffer = {};
	std::vector<float> valueBuffer = {};
	float SDF_V, gradSDFx_V, gradSDFy_V, gradSDFz_V; Vector3 gradSDF_V;

	if (!sdfGrid->hasGradient()) {
		sdfGrid->computeGradient();
	}

	for (int t = 0; t < NSteps; t++) {
		std::vector<Vector3> vNormals = evolvedSurface->getAngleWeightedVertexPseudoNormals();
		std::vector<std::vector<Vector3>> fvVerts = {};
		std::vector<std::vector<std::vector<uint>>> adjacentPolys = {};
		evolvedSurface->getVertexFiniteVolumes(&fvVerts, &adjacentPolys);
		
		this->initSystem(); // prep linear system of dim NVerts * NVerts

		for (uint i = 0; i < N; i++) {
			// ===== vertex =================================
			V.set(
				evolvedSurface->uniqueVertices[i].x,
				evolvedSurface->uniqueVertices[i].y,
				evolvedSurface->uniqueVertices[i].z
			);

			// ===== Control coefs ==========================
			// SDF value
			positionBuffer.clear(); // for old min and max positions
			valueBuffer.clear(); // for SDF cell vertex values
			sdfGrid->getSurroundingCells(V, Nx, Ny, Nz, sdfGrid->field, &positionBuffer, &valueBuffer);
			SDF_V = trilinearInterpolate(V, positionBuffer, valueBuffer);

			// -------- grad SDF --------
			// grad SDF x
			positionBuffer.clear(); // for old min and max positions
			valueBuffer.clear(); // for grad SDF x cell vertex values
			sdfGrid->getSurroundingCells(V, Nx - 2, Ny - 2, Nz - 2, sdfGrid->gradFieldX, &positionBuffer, &valueBuffer);
			gradSDFx_V = trilinearInterpolate(V, positionBuffer, valueBuffer);

			// grad SDF y
			positionBuffer.clear(); // for old min and max positions
			valueBuffer.clear(); // for grad SDF y cell vertex values
			sdfGrid->getSurroundingCells(V, Nx - 2, Ny - 2, Nz - 2, sdfGrid->gradFieldY, &positionBuffer, &valueBuffer);
			gradSDFy_V = trilinearInterpolate(V, positionBuffer, valueBuffer);

			// grad SDF z
			positionBuffer.clear(); // for old min and max positions
			valueBuffer.clear(); // for grad SDF z cell vertex values
			sdfGrid->getSurroundingCells(V, Nx - 2, Ny - 2, Nz - 2, sdfGrid->gradFieldZ, &positionBuffer, &valueBuffer);
			gradSDFz_V = trilinearInterpolate(V, positionBuffer, valueBuffer);

			gradSDF_V = Vector3(gradSDFx_V, gradSDFy_V, gradSDFz_V);

			// Laplace-Beltrami ctrl func: 
			float eps = 1.0f;			

			if (this->type == ElementType::tri) this->getTriangleVertexLaplaceBeltrami(i, &adjacentPolys[i]);
			else this->getQuadVertexLaplaceBeltrami(i, &adjacentPolys[i]);
		}
	}
}
