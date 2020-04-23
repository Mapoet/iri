#include "Iri.h"


// test
Iri::Iri()
{
	this->jmag = 0;
	for (size_t i = 0; i < 30; i++)
		this->jf.push_back(true);
	this->jf[3] = false;
	this->jf[4] = false;
	this->jf[5] = false;
	this->jf[20] = false;
	this->jf[27] = false;
	this->jf[28] = false;
	this->jf[29] = false;
	this->jf[32] = false;
	this->jf[33] = false;
	this->jf[34] = false;
	this->jf[38] = false;

}

Iri::Iri(int jmag, double xlat, double xlon, int iy, double mmdd, int iut,
	double hour, double height, double htec_max, int ivar, double vbeg,
	double vend, double vstp) {
}
