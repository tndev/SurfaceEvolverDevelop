#include "LinearSolver.h"

LinearSolver::LinearSolver(const unsigned int N)
{
	this->N = N;
}

LinearSolver::~LinearSolver()
{
}

void LinearSolver::Bi_CGSTAB_Solve(double** A, double* b, double* x, bool print)
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

	bool hasNaNs = false;

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
		if (std::isnan(norm)) hasNaNs = true;
		if (print) std::cout << "||s|| = " << norm << std::endl;
		if (norm < tol || hasNaNs) {
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
		if (std::isnan(norm)) hasNaNs = true;
#pragma omp parallel
#pragma omp single
		if (print) std::cout << "||r[k + 1]|| = " << norm << std::endl;
		if (norm < tol || hasNaNs) {
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
		if (std::isnan(norm)) hasNaNs = true;
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

void LinearSolver::printArray1(std::string name, double* a, int printLim, bool inRow, std::ostream& out)
{
	if (inRow) {
		out << name << " = " << std::endl;
		for (int i = 0; i < printLim; i++) {
			out << " " << a[i];
		}
		out << "  ... ";
		for (int i = N - printLim; i < N; i++) {
			out << " " << a[i];
		}
		out << std::endl;
	}
	else {
		std::string offset = std::string((name + " = ").length() + 1, ' ');
		out << name << " = " << std::endl;
		for (int i = 0; i < printLim; i++) {
			out << offset << a[i] << std::endl;
		}
		for (int i = 0; i < 3; i++) out << offset << "  ." << std::endl;
		for (int i = N - printLim; i < N; i++) {
			out << offset << a[i] << std::endl;
		}
		out << std::endl;
	}
}

void LinearSolver::printArray2(std::string name, double** A, int printLim, std::ostream& out)
{
	std::string offset = std::string((name + " = ").length() + 1, ' ');
	out << name << " = " << std::endl;
	for (int i = 0; i < printLim; i++) {
		out << offset;
		for (int j = 0; j < printLim; j++) {
			out << A[i][j] << " ";
		}
		out << "  ...  ";
		for (int j = N - printLim; j < N; j++) {
			out << A[i][j] << " ";
		}
		out << std::endl;
	}
	for (int i = 0; i < 3; i++) out << offset << "  ." << std::endl;
	for (int i = N - printLim; i < N; i++) {
		out << offset;
		for (int j = 0; j < printLim; j++) {
			out << A[i][j] << " ";
		}
		out << "  ...  ";
		for (int j = N - printLim; j < N; j++) {
			out << A[i][j] << " ";
		}
		out << std::endl;
	}
	out << std::endl;
}

double LinearSolver::vectorDot(double* a, double* b)
{
	double result = 0.;
	for (int i = 0; i < N; i++) result += a[i] * b[i];
	return result;
}

double LinearSolver::vectorNorm(double* a)
{
	return sqrt(vectorDot(a, a));
}
