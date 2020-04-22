#include "subroutines.h"

const std::string ITEXT[] = { "H", "LATI", "LONG", "YEAR", "MONTH", "DAY", "DAYOF", "HOUR" };

int main() {
	int jm, jmag, iy, imd, mmdd, iut, ivar,
		iymst, iymend, n, numstp,
		jchoice, piktab, nummax = 1000;
	double xlat, xlon, hour, phour, hx, hxx, htec_max,
		vbeg, vend, vstp, var, f107d, pf107d;
	unsigned int *jf = new unsigned int[50];
	std::string pname[7];
	std::string popt, map, hopt, bopt, iopt, dopt, sopt, seopt, topt, tsopt, f1opt, rzopt, igopt, fdopt, f8opt;

	double **oar = new double*[100];
	double **outf = new double*[20];
	for (size_t i = 0; i < 100; i++) {
		oar[i] = new double[1000];
		oar[i][0] = -1.0;
		if (i < 20)
			outf[i] = new double[1000];
	}

	read_ig_rz();
	readapf107();

	// user input of IRI input parameters
	std::cout << "Coord. Type (=0/1, geof/geom), "
		<< "Latitude(-90. - 90. deg.), "
		<< "Longitude(0. - 360. deg.):" << std::endl;
	std::cin >> jm >> xlat >> xlon;

	std::cout << "Year(yyyy), mmdd(or -ddd), iut(=0/1, LT/UT), hour:" 
		<< std::endl;
	std::cin >> iy >> imd >> iut >> hour;

	std::cout << "Height (km, from 60. to 2000.):" << std::endl;
	std::cin >> hx;

	std::cout << "Output-Option:" << std::endl;
	std::cout << "\tPlasma frequencies, B0, M3000, valley, width and depth"
		<< std::endl;
	piktab = 2;

	std::cout << "Upper height [km] for TEC integration (0 for no TEC):"
		<< std::endl;
	std::cin >> htec_max;

	std::cout << "Variable? (1/2/../8 for height/lat/long/year/month/day/day of year/hour):"
		<< std::endl;
	std::cin >> ivar;
	ivar--;

	std::cout << "Begin, end, and stepsize for the selected variable:"
		<< std::endl;
	std::cin >> vbeg >> vend >> vstp;

	std::cout << "Options: 1 (true) or 0 (false)" << std::endl;
	std::cout << "\tEnter 0 to use standard or 1 to enter your own"
		<< '\t' << std::endl;
	std::cin >> jchoice;

	for (size_t i = 0; i < 50; i++)
		jf[i] = 1;

	if (jchoice == 0) {
		jf[3] = 0; // t=B0table f=other models (0)
		jf[4] = 0; // t=CCIR  f=URSI foF2 model (0)
		jf[5] = 0; // t=DS95+DY85   f=RBV10+TTS03 (0)
		jf[20] = 0; // f=ion drift not computed (0)
		jf[22] = 0; // t=AEROS/ISIS f=TTS Te with PF10.7 (0)
		jf[27] = 0; // f=spread-F not computed (0)
		jf[28] = 0; // t=old  f=New Topside options (0)
		jf[29] = 0; // t=corr f=NeQuick topside (0)
		jf[32] = 0; // f=auroral boundary off (0)
		jf[33] = 0; // t=messages on f= off (0)
		jf[34] = 0; // f=auroral E-storm model off (0)
		jf[38] = 0; // t=M3000F2 model f=new hmF2 models (0)
	}
	else {
		std::cout << "\tCompute Ne, T, Ni? (enter: 1 1 1  if you want all):"
			<< '\t' << std::endl;
		std::cin >> jf[0] >> jf[1] >> jf[2];
		if (jf[0]) {
			std::cout << "\tLAY version: 1=standard ver., 0=LAY version. {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[10];

			std::cout << "\tNe Topside: 1=IRI-2001/h0.5, 0=new options {0}:"
				<< '\t' << std::endl;
			std::cin >> jf[28];

			std::cout << "\tNe Topside: 1=IRI01_corrt, 0=NeQuick/h0.5 {0}:"
				<< '\t' << std::endl;
			std::cin >> jf[29];

			std::cout << "\tNe Topside: 1=F10.7<188, 0=unlimited {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[6];

			std::cout << "\tNe topside: 1=w/o foF2 storm model, 0=with {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[36];

			std::cout << "\tfoF2 model: 1=CCIR, 0=URSI-88 {0}:"
				<< '\t' << std::endl;
			std::cin >> jf[4];

			std::cout << "\tfoF2: 1=with storm model, 0=without {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[25];

			std::cout << "\thmF2: 1=f(M3000F2), 0=new models {0}:"
				<< '\t' << std::endl;
			std::cin >> jf[38];

			std::cout << "\thmF2: 1=AMTB-model, 0=Shubin-COSMIC model {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[39];

			if (jf[39]) {
				std::cout << "\thmF2: 1=w/o foF2 storm model, 0=with {1}:"
					<< '\t' << std::endl;
				std::cin >> jf[35];
			}

			std::cout << "\tF2 peak density or foF2: 1=model, 0=user input {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[7];

			std::cout << "\tF2 peak height or M3000F2: 1=model, 0=user input {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[8];

			std::cout << "\tAuroral boundary model: 1=on, 0=off {0}:"
				<< '\t' << std::endl;
			std::cin >> jf[32];

			std::cout << "\tBottomside thickness B0: 1=Bil-2000, 0=other options {0}:"
				<< '\t' << std::endl;
			std::cin >> jf[3];

			std::cout << "\tBottomside thickness B0: 1=ABT-2009, 0=Gul-1987 {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[30];

			std::cout << "\tBottomside thickness B0: 1=model, 0=user input {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[42];

			std::cout << "\tF1 peak density or foF1: 1=model, 0=user input {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[12];

			if (!jf[10])
				std::cout << "\tF1 peak height: 1=model, 0=user input {1}:"
					<< '\t' << std::endl;
				std::cin >> jf[13];

			std::cout << "\tF1: 1=with probability model, 0=without {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[18];

			std::cout << "\tF1: 1=standard probability, 0=with L condition {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[19];

			std::cout << "\tE peak density or foE: 1=model, 0=user input {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[14];

			std::cout << "\tE peak height: 1=model, 0=user input {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[15];

			std::cout << "\tE peak auroral storm model: 1=on, 0=off {0}:"
				<< '\t' << std::endl;
			std::cin >> jf[34];

			std::cout << "\tD: 1=IRI-1990, 0=FT-2001 {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[23];
		}
		if (jf[1]) {
			std::cout << "\tTe(Ne) model: 1=not used, 0=correlation is used. {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[9];

			std::cout << "\tTe: 1=Bil-1985, 0=TBT-2012 {0}:"
				<< '\t' << std::endl;
			std::cin >> jf[22];

			std::cout << "\tTe: 1=TBT-2012 with PF107 dep., 0=w/o {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[41];
		}
		if (jf[2]) {
			std::cout << "\tIon comp. model: 1=DS95/DY85, 0=RBV10/TTS05 {0}:"
				<< '\t' << std::endl;
			std::cin >> jf[5];

			if (!jf[5]) {
				std::cout << "\tIRIFLIP: 1=messages off, 0=on {1}:"
					<< '\t' << std::endl;
				std::cin >> jf[37];
			}

			std::cout << "\tNi: 1=ion composition in %, 0=ion densities in cm-3 {1}:"
				<< '\t' << std::endl;
			std::cin >> jf[21];
		}
		std::cout << "\tEquat. Vert. Ion Drift: 1=computed, 0=not computed {1}:"
			<< '\t' << std::endl;
		std::cin >> jf[20];

		std::cout << "\tSpread-F probability: 1=computed, 0=not computed {1}:"
			<< '\t' << std::endl;
		std::cin >> jf[27];

		std::cout << "\tCOV: 1: COV=F10.7_365, 0: COV=func(IG12). {1}:"
			<< '\t' << std::endl;
		std::cin >> jf[40];

		std::cout << "\tSunspot index: 1=from file, 0=user input. {1}:"
			<< '\t' << std::endl;
		std::cin >> jf[16];

		std::cout << "\tIonospheric index: 1=from file, 0=user input. {1}:"
			<< '\t' << std::endl;
		std::cin >> jf[26];

		std::cout << "\tF10.7D Index: 1=from file, 0=user input {1}:"
			<< '\t' << std::endl;
		std::cin >> jf[24];

		std::cout << "\tF10.7_81 Index: 1=from file, 0=user input {1}:"
			<< '\t' << std::endl;
		std::cin >> jf[31];

		std::cout << "\tdip, magbr, modip: t=IGRF, f=old FIELDG using "
			<< "POGO68/10 for 1973 {1}:" << '\t' << std::endl;
		std::cin >> jf[17];

		std::cout << "\tMessages 1=on, 0=off {1}:"
			<< '\t' << std::endl;
		std::cin >> jf[33];

		std::cout << "\tMessage output unit: 1=(UNIT=6), 0=(UNIT=11). {1}:"
			<< '\t' << std::endl;
		std::cin >> jf[11];
	}
	
	// option to enter measured values for NmF2, hmF2, NmF1, hmF1, NmE, hmE,
    // B0, N(300), N(400), N(600) if available;
	std::cout << std::endl << std::endl << std::endl;
	numstp = (vend - vbeg) / vstp + 1;
	if (ivar == 0)
		numstp = 1;
	if (jf[0]) {
		if (!jf[7] || !jf[8] || !jf[12] || !jf[13] || !jf[14] || !jf[15] || !jf[42]) {
			var = vbeg;
			int i = 0;
			do {
				if (!jf[7]) {
					jf[25] = 0;
					std::cout << "foF2/Mhz or NmF2/m-3 for " << ITEXT[ivar] << "="
						<< var << std::endl;
					std::cin >> oar[0][i];
					pname[0] = "foF2/MHz";
					if (oar[0][i] > 30)
						pname[0] = "NmF2/m-3";
				}
				if (!jf[8]) {
					std::cout << "hmF2/km or M3000F2 for " << ITEXT[ivar] << "="
						<< var << std::endl;
					std::cin >> oar[1][i];
					pname[1] = "M(3000)F2";
					if (oar[1][i] > 50)
						pname[1] = "hmF2/km";
				}
				if (!jf[12]) {
					std::cout << "foF1/MHz or NmF1/m-3 for " << ITEXT[ivar] << "="
						<< var << std::endl;
					std::cin >> oar[2][i];
					pname[2] = "foF1/MHz";
					if (oar[2][i] > 30)
						pname[2] = "NmF1/m-3";
				}
				if (!jf[13]) {
					std::cout << "hmF1/km for " << ITEXT[ivar] << "="
						<< var << std::endl;
					std::cin >> oar[3][i];
					pname[3] = "hmF1/km";
				}
				if (!jf[14]) {
					std::cout << "foE/MHz or NmE/m-3 for " << ITEXT[ivar] << "="
						<< var << std::endl;
					std::cin >> oar[4][i];
					pname[4] = "foE/MHz";
					if (oar[4][i] > 30)
						pname[4] = "NmE/m-3";
				}
				if (!jf[15]) {
					std::cout << "hmE/km for " << ITEXT[ivar] << "="
						<< var << std::endl;
					std::cin >> oar[5][i];
					pname[5] = "hmE/km";
				}
				if (!jf[42]) {
					std::cout << "B0/km for " << ITEXT[ivar] << "="
						<< var << std::endl;
					std::cin >> oar[9][i];
					pname[6] = "B0/km ";
				}
				if (!jf[43]) {
					std::cout << "B1 for " << ITEXT[ivar] << "="
						<< var << std::endl;
					std::cin >> oar[86][i];
					pname[6] = "B1";
				}
				i = i + 1;
				var = var + vstp;
			} while (ivar > 1 && var <= vend);
		}
	}
	// option to enter Ne for Te-Ne relationship
	if (jf[1] && !jf[9]) {
		var = vbeg;
		for (size_t i = 0; i < numstp; i++) {
			std::cout << "Ne(300km),Ne(400km)/m-3 for " << ITEXT[ivar]
				<< '=' << var <<  " [-1 if not]" << std::endl;
			std::cin >> oar[14][i] >> oar[15][i];
			var = var + vstp;
		}
	}
	// option to enter F107D and/or PF107
	if (!jf[24]) {
		std::cout << "User input for F107D:" << std::endl;
		std::cin >> f107d;
		for (size_t i = 0; i < 100; i++)
			oar[40][i] = f107d;
	}
	if (!jf[31]) {
		std::cout << "User input for PF107:" << std::endl;
		std::cin >> pf107d;
		for (size_t i = 0; i < 100; i++)
			oar[45][i] = pf107d;
	}
	// option to enter Rz12 and/or IG12
	if (!jf[16]) {
		std::cout << "User input for Rz12:" << std::endl;
		std::cin >> oar[32][0];
		for (size_t i = 1; i < 100; i++)
			oar[32][i] = oar[32][0];
	}
	if (!jf[26]) {
		std::cout << "User input for IG12:" << std::endl;
		std::cin >> oar[38][0];
		for (size_t i = 1; i < 100; i++)
			oar[38][i] = oar[38][0];
	}
	if (!jf[31]) {
		std::cout << "User input for F10.7_81d:" << std::endl;
		std::cin >> oar[45][0];
		for (size_t i = 1; i < 100; i++)
			oar[45][i] = oar[45][0];
	}
	// end of user input

	numstp = abs((vend - vbeg) / vstp + 1);
	if (numstp > nummax)
		numstp = nummax;

	if (jf[28])
		popt = "IRI2001";
	else
		if (jf[29])
			popt = "IRIcorr";
		else
			popt = "NeQuick";

	map = "URSI";
	if (jf[4])
		map = "CCIR";

	if (jf[38])
		hopt = "CCIR-M3000";
	else
		if (jf[39])
			hopt = "AMTB-2013";
		else
			hopt = "Shubin2015";

	if (jf[3])
		bopt = "BIl-2000";
	else
		if (jf[30])
			bopt = "ABT-2009";
		else
			bopt = "Gul-1987";

	iopt = "RBV10+TTS05";
	if (jf[5])
		iopt = "DS95 + DY85";

	dopt = "FT01+DRS95";
	if (jf[23])
		dopt = "IRI-1990";

	sopt = "off";
	if (jf[25])
		sopt = "on ";

	seopt = "off";
	if (jf[34])
		seopt = "on ";

	topt = "TBT-2012";
	tsopt = " with";
	if (jf[22]) {
		topt = "BIl-1985";
		tsopt = "  w/o";
	}
	
	if (jf[18]) {
		f1opt = "Scotto-97 no L";
		if (!jf[19])
			f1opt = "Scotto-97 with L";
		else
			f1opt = "IRI-95";
	}

	rzopt = " user input";
	if (jf[16])
		rzopt = " ";
	igopt = " user input";
	if (jf[26])
		igopt = " ";
	fdopt = " user input";
	if (jf[24])
		fdopt = " ";
	f8opt = " user input";
	if (jf[31])
		f8opt = " ";

	hxx = hx;
	jmag = jm;
	mmdd = imd;

	//calling IRI subroutine
	iri_web(jmag, jf, xlat, xlon, iy, mmdd, iut, hour, 
		hxx, htec_max, ivar, vbeg, vend, vstp, outf, oar);

	//preparation of results page
	// TODO: !

	return 0;
}