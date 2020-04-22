#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <functional>
#include "constants.h"
#include "commons.h"


// irisub.for
void iri_sub(unsigned int *jf, int jmag, double alati,
	double along, int iyyyy, int mmdd, double dhour,
	double heibeg, double heiend, double heistp,
	double **outf, double *oarr);

void iri_web(int jmag, unsigned int *jf, double alati,
	double along, int iyyyy, int mmdd, int iut,
	double dhour, double height, double h_tec_max,
	int ivar, double vbeg, double vend, double vstp,
	double **a, double **b);


// iritec.for
void iri_tec(double hstart, double hend, int istep,
	double &tectot, double &tectop, double &tecbot);


// irifun.for
void read_ig_rz();

void readapf107();

void moda(int in, int iyear, int &month, int &iday,
	int &idoy, int &nrdaymo);

void fieldg(double dlat, double dlong, double alt,
	double &x, double &y, double &z, double &f, double &dip,
	double &dec, double &smodip);

void sun(int iyear, int iday, int ihour, int min, int isec,
	double &gst, double &slong, double &srasn, double &sdec);

void tcon(int yr, int mm, int day, int idn, double *rz,
	double *ig, double &rsn, int &nmonth);

void apf_only(int iyyyy, int imn, int id, double &f107d, double &f107pd,
	double &f107_81, double &f107_365, int &iapda, int &isdate);

void soco(int ld, double t, double flat, double elon, double height,
	double &declin, double &zenith, double &sunrse, double &sunset);

double foeedi(double cov, double xhi, double xhim, double xlati);

void storm(int *ap, double rga, double rgo, int coor, double &rgma, int ut,
	int doy, double &cf);

double conver(double rga, double rgo);

void regfa1(double x11, double x22, double fx11, double fx22,
	double eps, double fw, double(*f)(double), bool schalt, double &x);

void shamdhmf2(double rlat, double flon, double t, double rz, double &hmf2);

void schnevpdh(double rz, double flat, double flon, double r, double t,
	int l, double &bn, double &be, double &bv);

void schnevpd(double rz, double flat, double flon, double r, double t,
	int l, double &bn, double &be, double &bv);

void tbfit(double t1, double t2, int ibf,
	double &thint, double &tzero);

void legfun(int m, double fn, double constp, double colat, 
	double &p, double &dp, double &pms);

double model_hmf2(int day, int month, double ut, double xmodip, double xlong,
	double f107_81);

double sdmf2(double ut, int monthut, double f107a, double xmodip, double xlong);

double hmf2_med_sd(int iut, int monthut, double f107a, double xmodip, double xlong);

void read_data_sd(int month, double **coeff_month);

double fun_hmf2_sd(double teta, double xlong, double *kf);

void fun_gk(double teta, double xlong, double *gk);

void legendre(int mm, int nn, double **p, double teta);

double fun_hmf2ut(double t);

void koeff_ut(int mm, int mk, double *kf_ut);

void fun_akp_ut(int mm, int mk, double **akp_ut, double *dk_ut);

void fun_fk_ut(int mk, double *gk_ut, double **akp_ut, double *fk_ut);

void fun_gk_ut(int mm, int mk, double t, double *gk_ut);

double b0_98(double hour, double sax, double sux, int nseasn, double r,
	double zlo, double zmodip);

double eptr(double x, double sc, double hx);

//ROGUL(SEADAY,XHI,SEAX,GRAT)

void shamdb0d(double rlat, double flon, double t, double rz, double &b);

void shab1d(double flat, double flon, double t, double rz, double &b);

void schnevpdb1(double rz, double flat, double flon, double r, double t,
	int l, double &bn, double &be, double &bv);

double fof1ed(double ylati, double r, double chi);

double f1_c1(double xmodip, double hour, double suxnon, double saxnon);

void f1_prob(double sza, double glat, double rz12,
	double &f1prob, double &f1probl);

double hpol(double hour, double tw, double xnw, double sa, double su,
	double dsa, double dsu);

double epst(double x, double sc, double hx);

void tal(double shabr, double sdelta, double shbr, double sdtdh0,
	bool aus6, double *spt);

double xmded(double xhi, double r, double yw);

void apfmsis(int isdate, double hour, double *iapo);


// igrf.for
void sphcar(double &r, double &teta, double &phi,
	double &x, double &y, double &z, int j);

void geodip(int &iyr, double &sla, double &slo,
	double &dla, double &dlo, int j);

void geomag(double &xgeo, double &ygeo, double &zgeo, double &xmag,
	double &ymag, double &zmag, int j, int iyr);

void recalc(int iyr, int iday, int ihour, int min, int isec);

void getshc(std::string fspec, int &nmax, double &erad, double *gh);

void intershc(double date, double dte1, int nmax1, double *gh1, double dte2, int nmax2,
	double *gh2, int &nmax, double *gh);

void extrashc(double date, double dte1, int nmax1, double *gh1, int nmax2,
	double *gh2, int &nmax, double *gh);

void feldcof(double year);

void igrf_dip(double xlat, double xlong, double year, double height,
	double &dec, double &dip, double &dipl, double &ymodip);

void feldg(double glat, double glon, double alt,
	double &bnorth, double &beast, double &bdown, double &babs);

void feldi();

void igrf_sub(double xlat, double xlong, double year, double height,
	double &xl, int &icode, double &dipl, double &babs);

void shellg(double glat, double glon, double alt,
	double &fl, int &icode, double &b0);

void stoer(double **p, int idx, double &bq, double &r);

void clcmlt(int iyyyy, int ddd, double uthr, double glat, double glon, double &mlt);

void dpmtrx(int iyyyy, int ddd, double *xm, double *ym, double *zm);

double fmodip(double xlat);


// cira.for

void tselec(double *sv);

void gtd7(int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
	double f107, double *ap, int mass, double *d, double *t);

double vtst7(int iyd, double sec, double glat, double glong, double stl, double f107a,
	double f107, double *ap, int ic);

void glatf(double lat, double &gv, double &reff);

void gts7(int iyd, double sec, double alt, double glat, double glong, double stl,
	double f107a, double f107, double *ap, int mass, double *d, double *t);

double globe7(double yrd, double sec, double xlat, double xlong, 
	double tloc, double f107a, double f107, double *ap, double p[]);

double glob7s(double *p);


// TODO: invdpc,