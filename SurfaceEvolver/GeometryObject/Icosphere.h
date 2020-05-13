#ifndef ICOSPHERE_H_
#define ICOSPHERE_H_

#include "../Geometry/Geometry.h"

class IcoSphere :
	public Geometry
{
public:
	unsigned int detail = 0;
	double radius = 50.0;
	IcoSphere();
	IcoSphere(const IcoSphere& other);
	IcoSphere(unsigned int detail, double radius, std::string name = "");
	~IcoSphere();

	void build();
};

#endif