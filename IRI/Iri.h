#pragma once
#include <vector>


class Iri
{
	int jmag; // geog/geom
	double xlat;
	double xlon;
	int iy; // year
	double mmdd; // mmdd or -ddd
	int iut; // local or universal time
	double hour;
	double height; // in km
	double htec_max; // =0 no TEC otherwise upper boundary for integral
	double vbeg;
	double vend;
	double vstp;
	//  =1 altitude
	//  =2,3 latitude, longitude
	//  =4,5,6 year, month, day		
	//  =7 day of year
	//  =8 hour(UT or LT)
	int ivar; 
	// options
	std::vector<bool> jf;

	// call iri_web(jmag,jf,xlat,xlon,iy,mmdd,iut,hour,
	//&hxx, htec_max, ivar, vbeg, vend, vstp, outf, oar)
public:
	Iri(int jmag, double xlat, double xlon, int iy, double mmdd, int iut,
		double hour, double height, double htec_max, int ivar, double vbeg,
		double vend, double vstp);
	Iri();
};

