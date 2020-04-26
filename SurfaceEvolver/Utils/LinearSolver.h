#ifndef LINEARSOLVER_H_
#define LINEARSOLVER_H_

#include <iostream>

class LinearSolver
{
public:
	// Bi_CGSTAB_Solve
	LinearSolver(const unsigned int N);
	~LinearSolver();

	void Bi_CGSTAB_Solve(double** A, double* b, double* x, bool print = false);

	// output
	void printArray1(std::string name, double* a, int printLim, bool inRow = true, std::ostream& out = std::cout);
	void printArray2(std::string name, double** A, int printLim, std::ostream& out = std::cout);
private:
	// system dimension
	unsigned int N;

	// solver helpers
	double vectorDot(double* a, double* b);
	double vectorNorm(double* a);
};

#endif