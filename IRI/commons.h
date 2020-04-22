#pragma once

#ifndef COMMONS_H
#define COMMONS_H
struct model_ {
	int nmax;
	double time;
	double g[196];
	std::string name;
} model;

struct c1_ {
	double st0;
	double ct0;
	double sl0;
	double cl0;
	double ctcl; 
	double stcl;
	double ctsl;
	double stsl; 
	double ab[19];
	int k;
	int iy;
	double bb[8];
} c1;

struct dipol_ {
	double ghi1, ghi2, ghi3;
} dipol;

struct igrf2_ {
	double xi[3];
	double x[3];
	double h[196];
} igrf2;

struct fidb0_ {
	double sp[3];
} fidb0;

struct igrz_ {
	double aig[806], arz[806];
	int iymst, iymend;
} igrz;

struct apfa_ {
	int aap[23000][9];
	double af107[23000][3];
	int n;
} apfa;

struct findrlat_ {
	double flon, ryear;
} findrlat;

struct amtb_ {
	double bint[9][9][5],
		bext[9][9][5],
		re,
		tzero;
	int ifit, ib, kint, lint, kext,
		lext, kmax;
	double fn[9][9];
} amtb;

struct hmf2ut_ {
	double hmf2_ut[24];
} hmf2ut;

struct csw_ {
	double sw[25];
	double swc[25];
	int isw;
} csw;

struct gts3c_ {
	double tlb, s, db04, db16, db28, db32, db40, db48, db01,
		za, t0, z0, g0, rl, dd, db14, tr12;
} gts3c;

struct meso7_ {
	double tn1[5], tn2[4], tn3[5], tgn1[2], tgn2[2], tgn3[2];
} meso7;

struct lower7_ {
	double ptm[10], pdm[10][8];
} lower7;

struct parm7_ {
	double pt[150], pda1[150], pda2[150], pda3[150],
		pda4[150], pda5[150], pda6[150], pda7[150], pda8[150], pda9[150],
		ps[150], pdl[25][2], ptl[100][4], pma[100][10], sam[100];
} parm7;

struct datim7_ {
	std::string isd[3], ist[2], nam[2];
} datim7;

struct datime_ {
	std::string isdate[3], istime[2], name[2];
} datime;

struct mavg7_ {
	double pavgm[10];
} mavg7;

struct dmix_ {
	double dm04, dm16, dm28, dm32, dm40, dm01, dm14;
} dmix;

struct parmb_ {
	double gsurf, re;
} paramb;

struct metsel_ {
	int imr;
} metsel;

struct ttest_ {
	double tinf, gb, rout, t[15];
} ttest;


struct lpoly_ {
	double plg[9][4], ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc,
		iyr, day, df, dfa, apd, apdf, apt[4], xlong;
} lpoly;
#endif