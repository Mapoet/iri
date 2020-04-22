#include "subroutines.h"

void read_ig_rz() {
//----------------------------------------------------------------
// Reads the Rz12 and IG12 indices file IG_RZ.DAT from I / O UNIT = 12
// and stores the indices in COMMON :
//		common / igrz / aig, arziyst, iyed   with aig(806), arz(806),
//											start year(iyst)
//                                           end year(iyed)
//
//	The indices file IG_RZ.DAT is structured as follows(values are
//	separated by comma) :
//  day, month, year of the last update of this file,
//  a blank line
//  start month, start year, end month, end year,
//  a blank line
//  the IG index for December of start year minus 1 (this value is
//		needed for interpolating from 1st to 15th of first year)
//  the 12 IG indices(13 - months running mean) for start year,
//  the 12 IG indices for the second year
//      .. and so on until the last year,
//  the 12 IG indices for the last year
//  the IG index for January of end year plus 1 (needed for interpolation)
//  a blank line
//  the Rz index for December of start year minus 1
//  the 12 Rz indices(13 - months running mean) for the start year,
//  the 12 Rz indices for the second year
//      .. and so on until the last year.
//  the 12 Rz indices for the last year
//  the Rz index for January of end year plus 1
//
//A negative Rz index means that the given index is the 13 - months -
//running mean of the solar radio flux(F10.7).The close correlation
//between(Rz)12 and (F10.7)12 is used to compute the(Rz)12 indices.
//
//An IG index of - 111 indicates that no IG values are available for the
//time period.In this case a correlation function between(IG)12 and
//(Rz)12 is used to obtain(IG)12.
//
//The computation of the 13 - month - running mean for month M requires the
//indices for the six months preceeding M and the six months following
//M(month: M - 6, ..., M + 6).To calculate the current running mean one
//therefore requires predictions of the indix for the next six months.
//Starting from six months before the UPDATE DATE(listed at the top of
//the file) and onward the indices are therefore based on indices
//predictions.
//----------------------------------------------------------------
	int iyst, iyend, iupd,
		iupm, iupy, imst,
		imend, inum_vals, inum_chan;

	std::string value;
	std::fstream myfile("data/ig_rz.dat", std::ios_base::in);

	if (myfile.is_open())
	{
		std::getline(myfile, value, ',');
		iupd = stoi(value);

		std::getline(myfile, value, ',');
		iupm = stoi(value);

		std::getline(myfile, value, ',');
		iupy = stoi(value);

		std::getline(myfile, value, ',');
		imst = stoi(value);

		std::getline(myfile, value, ',');
		iyst = stoi(value);

		std::getline(myfile, value, ',');
		imend = stoi(value);

		std::getline(myfile, value, ',');
		iyend = stoi(value);

		igrz.iymst = iyst * 100 + imst;
		igrz.iymend = iyend * 100 + imend;

		inum_vals = 3 - imst + (iyend - iyst) * 12 + imend;

		for (int i = 1; i <= inum_vals; i++) {
			std::getline(myfile, value, ',');
			igrz.aig[i] = stod(value);
		}
		for (int i = 1; i <= inum_vals; i++) {
			std::getline(myfile, value, ',');
			igrz.arz[i] = stod(value);
		}

		if ((iupy * 100 + iupm) > 201609)
		{
			inum_chan = 3 - imst + (2014 - iyst) * 12;
			for (size_t i = inum_chan; i <= inum_vals; i++)
				igrz.arz[i] = igrz.arz[i] * 0.7;
		}

		myfile.close();
	}
	else
		std::cout << "ERROR. FILE WASN'T OPENED!" << std::endl;
}

void readapf107() {
//------------------------------------------------------------------------ -
// Reads APF107.DAT file(on UNIT = 13) and stores contents in COMMON block :
// 	COMMON / AAP, AF107, N / with  AAP(23000, 9) and AF107(23000, 3)
//		AAP(*, 1)	3 - hour Ap indices for the UT interval )0 - 3)
//		AAP(*, 2)	3 - hour Ap indices for the UT interval )3 - 6)
//         ....                       ....
//		AAP(*, 8)	3 - hour Ap indices for the UT interval )21 - 6)
//		AAP(*, 9)	daily Ap
//		AF107(*, 1)	F10.7 radio flux for the day
//		AF107(*, 2)	81 - day average of F10.7 radio flux
//		AF107(*, 3)	365 - day average of F10.7
//      N           total number of records
//
//APF107.DAT is structured as follows :
// 		JY(I3), JMN(I3), JD(I3)	year, month, day
//		IIAP(8)	(8I3)			3 - hour Ap indices for the UT intervals
//(0 - 3(, (3 - 6(, (6 - 9(, .., (18 - 21(, (21 - 24(
//		IAPD(I3)				daily Ap
//		IR(I3)					sunspot number for the day(empty)
//		F107(F5.1)				F10.7 radio flux for the day
//		F107_81(F5.1)			81 - day average of F10.7 radio flux
//      F107_365(F5.1)         365 - day average of F10.7 centered on
//                              the date of interest.At start and end
//								of index file it takes all available
//                              indices, e.g. for the first date the
//                              average is only over 40 F10.7 values
//                              and over 41 values on the 2nd date.
//
//If date is outside the range of the Ap indices file then IAP(1) = -5
//------------------------------------------------------------------------ -
	std::fstream myfile("data/apf107.dat", std::ios::in);

	int /*JY, JMN, JD,*/ iapda, IR;
	int *iiap = new int[8];
	double F107D, F107_81, F107_365;
	std::string value;

	if (myfile.is_open()) {

		for (size_t i = 0; (i < 23000) && myfile.good(); i++) {
			std::getline(myfile, value);
			/*JY = stoi(value.substr(0, 3));
			JMN = stoi(value.substr(3, 3));
			JD = stoi(value.substr(6, 3));*/

			for (size_t k = 0; k < 8; k++)
				iiap[k] = stoi(value.substr(9 + k * 3, 3));

			iapda = stoi(value.substr(33, 3));
			
			IR = stoi(value.substr(36, 3));
			F107D = stod(value.substr(39, 5));
			F107_81 = stod(value.substr(44, 5));
			F107_365 = stod(value.substr(49, 5));

			for (size_t j = 0; j < 8; j++) {
				apfa.aap[i][j] = iiap[j];
			}
			apfa.aap[i][8] = iapda;
			if (F107_81 < -4) F107_81 = F107D;
			if (F107_365 < -4) F107_365 = F107D;
			apfa.af107[i][0] = F107D;
			apfa.af107[i][1] = F107_81;
			apfa.af107[i][2] = F107_365;
			apfa.n = i;
		}
		myfile.close();
	}

}

void iri_web(int jmag, unsigned int *jf, double alati,
	double along, int iyyyy, int mmdd, int iut,
	double dhour, double height, double h_tec_max,
	int ivar, double vbeg, double vend, double vstp,
	double **a, double **b) {
//------------------------------------------------------------------------
// input : jmag, alati, along, iyyyy, mmdd, dhour  see IRI_SUB
//          height  height in km
//          h_tec_max = 0 no TEC otherwise upper boundary for integral
//          iut = 1 for UT = 0 for LT
//          ivar = 1      altitude
// = 2, 3    latitude, longitude
// = 4, 5, 6  year, month, day
// = 7      day of year
// = 8      hour(UT or LT)
//          vbeg, vend, vstp  variable range(begin, end, step)
// output : a       similar to outf in IRI_SUB
//          b       similar to oarr in IRI_SUB
//
//          numstp  number of steps; maximal 1000
//------------------------------------------------------------------------
	double xhour, tec, tecb, tect;
	double *xvar = new double[8];
	double **outf = new double*[20];
	for (size_t i = 0; i < 20; i++)
		outf[i] = new double[1000];
	int nummax = 1000;
	int numstp = (vend - vbeg) / vstp + 1;

	if (numstp > nummax)
		numstp = nummax;

	double *oar = new double[100];
	for (size_t i = 0; i < 100; i++)
		oar[i] = b[i][0];

	double *oarr = new double[100];
	if (ivar == 0) {
		for (size_t i = 0; i < 100; i++)
			oarr[i] = oar[i];
		xhour = dhour + iut * 25;

		iri_sub(jf, jmag, alati, along, iyyyy, mmdd, xhour, vbeg, vend, vstp, a, oarr);
		
		if (h_tec_max > 50) {
			iri_tec(50.0, h_tec_max, 2, tec, tect, tecb);
			oarr[36] = tec;
			oarr[37] = tect;
		}
		for (size_t i = 0; i < 100; i++)
			b[i][0] = oarr[i];

		return;
	}

	if (height < 0.0)
		height = 100;

	xvar[1] = alati;
	xvar[2] = along;
	xvar[3] = iyyyy;
	xvar[4] = mmdd / 100;
	xvar[5] = mmdd % 100;
	xvar[6] = abs(mmdd);
	xvar[7] = dhour;

	xvar[ivar] = vbeg;

	alati = xvar[1];
	along = xvar[2];
	iyyyy = int(xvar[3]);

	if (ivar == 6)
		mmdd = -int(vbeg);
	else
		mmdd = int(xvar[4] * 100 + xvar[5]);

	dhour = xvar[7] + iut * 25.0;

	for (size_t i = 0; i < numstp; i++)
	{
		for (size_t j = 0; j < 100; j++)
			oarr[j] = b[j][i];

		iri_sub(jf, jmag, alati, along, iyyyy, mmdd, dhour, height, height, 1.0, outf, oarr);

		if (h_tec_max > 50) {
			iri_tec(50.0, h_tec_max, 2, tec, tect, tecb);
			oarr[36] = tec;
			oarr[37] = tect;
		}

		for (size_t j = 0; j < 100; j++) {
			b[j][i] = oarr[j];
			if (j < 20)
				a[j][i] = outf[j][0];
		}

		xvar[ivar] += vstp;

		alati = xvar[1];
		along = xvar[2];
		iyyyy = int(xvar[3]);

		if (ivar == 6)
			mmdd = -int(vbeg);
		else
			mmdd = int(xvar[4] * 100 + xvar[5]);

		dhour = xvar[7] + iut * 25.0;
	}
	return;
}

void iri_sub(unsigned int *jf, int jmag, double alati,
	double along, int iyyyy, int mmdd, double dhour,
	double heibeg, double heiend, double heistp,
	double **outf, double *oarr) {
/*/---------------------------------------------------------------- -
//
// INPUT : JF(1:50)      true / false switches for several options
//         JMAG = 0 geographic = 1 geomagnetic coordinates
//         ALATI, ALONG   LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
//         IYYYY         Year as YYYY, e.g. 1985
//         MMDD(-DDD)   DATE(OR DAY OF YEAR AS A NEGATIVE NUMBER)
//         DHOUR         LOCAL TIME(OR UNIVERSAL TIME + 25) IN DECIMAL
//                          HOURS
//         HEIBEG, HEIGHT RANGE IN KM; maximal 100 heights, i.e.
//          HEIEND, HEISTP        int((heiend - heibeg) / heistp) + 1.le.100
//
//    JF switches to turn off / on(.true. / .false.) several options
//
//    i       .true.                  .false.standard version
//---------------------------------------------------------------- -
//    1    Ne computed            Ne not computed                     t
//    2    Te, Ti computed        Te, Ti not computed                 t
//    3    Ne & Ni computed       Ni not computed                     t
//    4    B0, B1 - Bil - 2000       B0, B1 - other models jf(31)     false
//    5    foF2 - CCIR            foF2 - URSI                     false
//    6    Ni - DS - 1995 & DY - 1985 Ni - RBV - 2010 & TTS - 2005        false
//    7    Ne - Tops: f10.7 < 188   f10.7 unlimited                     t
//    8    foF2 from model        foF2 or NmF2 - user input           t
//    9    hmF2 from model        hmF2 or M3000F2 - user input        t
//   10    Te - Standard          Te - Using Te / Ne correlation        t
//   11    Ne - Standard Profile  Ne - Lay - function formalism         t
//   12    Messages to unit 6     to messages.txt on unit 11          t
//   13    foF1 from model        foF1 or NmF1 - user input           t
//   14    hmF1 from model        hmF1 - user input(only Lay version)t
//   15    foE  from model        foE or NmE - user input             t
//   16    hmE  from model        hmE - user input                    t
//   17    Rz12 from file         Rz12 - user input                   t
//   18    IGRF dip, magbr, modip old feldg using POGO68 / 10 for 1973 t
//   19    F1 probability model   critical solar zenith angle(old)   t
//   20    standard F1            standard F1 plus L condition        t
//   21    ion drift computed     ion drift not computed          false
//   22    ion densities in %     ion densities in m - 3                t
//   23    Te_tops(Bil - 1985)     Te_topside(TBT - 2012)           false
//   24    D - region: IRI - 1990     FT - 2001 and DRS - 1995                t
//   25    F107D from APF107.DAT  F107D user input(oarr(41))         t
//   26    foF2 storm model       no storm updating                   t
//   27    IG12 from file         IG12 - user                         t
//   28    spread - F probability 	 not computed                    false
//   29    IRI01 - topside          new options as def.by JF(30)   false
//   30    IRI01 - topside corr.NeQuick topside model   	     false
//(29, 30) = (t, t) IRIold, (f, t) IRIcor, (f, f) NeQuick
//   31    B0, B1 ABT - 2009	     B0 Gulyaeva - 1987 h0.5               t
//(4, 31) = (t, t) Bil - 00, (f, t) ABT - 09, (f, f) Gul - 87, (t, f) not used
//   32    F10.7_81 from file     F10.7_81 - user input(oarr(46))    t
//   33    Auroral boundary model on / off  true / false	             false
//   34    Messages on            Messages off                        t
//   35    foE storm model        no foE storm updating           false
//   36    hmF2 w / out foF2_storm  with foF2 - storm                     t
//   37    topside w / out foF2 - storm  with foF2 - storm                  t
//   38    turn WRITEs off in IRIFLIP   turn WRITEs on                t
//   39    hmF2(M3000F2)         new models                      false
//   40    hmF2 AMTB - model        Shubin - COSMIC model                 t
//   41    Use COV = F10.7_365      COV = f(IG12) (IRI before Oct 2015)   t
//   42    Te with PF10.7 dep.w / o PF10.7 dependance               t
//   43    B0 from model          B0 user input                       t
//   44    B1 from model          B1 user input                       t
//      ....
//   50
//------------------------------------------------------------------
//
//  Depending on the jf() settings additional INPUT parameters may
//  be required :
//
//       Setting              INPUT parameter
//---------------------------------------------------------------- -
//    jf(8) = .false.OARR(1) = user input for foF2 / MHz or NmF2 / m - 3
//    jf(9) = .false.OARR(2) = user input for hmF2 / km or M(3000)F2
//    jf(10) = .false.OARR(15), OARR(16) = user input for Ne(300km),
//       Ne(400km) / m - 3. Use OARR() = -1 if one of these values is not
//       available.If jf(23) = .false.then Ne(300km), Ne(550km) / m - 3.
//    jf(13) = .false.OARR(3) = user input for foF1 / MHz or NmF1 / m - 3
//    jf(14) = .false.OARR(4) = user input for hmF1 / km
//    jf(15) = .false.OARR(5) = user input for foE / MHz or NmE / m - 3
//    jf(16) = .false.OARR(6) = user input for hmE / km
//    jf(17) = .flase.OARR(33) = user input for Rz12
//    jf(25) = .false.OARR(41) = user input for daily F10.7 index
//    jf(27) = .false.OARR(39) = user input for IG12
//    jf(43) = .false.OARR(10) = user input for B0
//    jf(44) = .false.OARR(87) = user input for B1
//
//
//  OUTPUT : OUTF(1:20, 1 : 1000)
//               OUTF(1, *)  ELECTRON DENSITY / M - 3
//               OUTF(2, *)  NEUTRAL TEMPERATURE / K
//               OUTF(3, *)  ION TEMPERATURE / K
//               OUTF(4, *)  ELECTRON TEMPERATURE / K
//               OUTF(5, *)  O + ION DENSITY / % or /M - 3 if jf(22) = f
//               OUTF(6, *)  H + ION DENSITY / % or /M - 3 if jf(22) = f
//               OUTF(7, *)  HE + ION DENSITY / % or /M - 3 if jf(22) = f
//               OUTF(8, *)  O2 + ION DENSITY / % or /M - 3 if jf(22) = f
//               OUTF(9, *)  NO + ION DENSITY / % or /M - 3 if jf(22) = f
//                 AND, IF JF(6) = .FALSE.:
//               OUTF(10, *)  CLUSTER IONS DEN / % or /M - 3 if jf(22) = f
//               OUTF(11, *)  N + ION DENSITY / % or /M - 3 if jf(22) = f
//               OUTF(12, *)
//               OUTF(13, *)
//  if (jf(24)    OUTF(14, 1:11) standard IRI - Ne for 60, 65, .., 110km
	// = .false.)        12:22) Friedrich(FIRI) model at these heights
//                      23:33) standard Danilov(SW = 0, WA = 0)
//                      34:44) for minor Stratospheric Warming(SW = 0.5)
//                      45:55) for major Stratospheric Warming(SW = 1)
//                      56:66) weak Winter Anomaly(WA = 0.5) conditions
//                      67:77) strong Winter Anomaly(WA = 1) conditions
//               OUTF(15 - 20, *)  free
//
//            OARR(1:100)   ADDITIONAL OUTPUT PARAMETERS
//
//      #OARR(1) = NMF2 / M - 3           #OARR(2) = HMF2 / KM
//      #OARR(3) = NMF1 / M - 3           #OARR(4) = HMF1 / KM
//      #OARR(5) = NME / M - 3            #OARR(6) = HME / KM
//       OARR(7) = NMD / M - 3             OARR(8) = HMD / KM
//       OARR(9) = HHALF / KM           #OARR(10) = B0 / KM
//       OARR(11) = VALLEY - BASE / M - 3     OARR(12) = VALLEY - TOP / KM
//       OARR(13) = TE - PEAK / K          OARR(14) = TE - PEAK HEIGHT / KM
//      #OARR(15) = TE - MOD(300KM)     #OARR(16) = TE - MOD(400KM) / K
//       OARR(17) = TE - MOD(600KM)      OARR(18) = TE - MOD(1400KM) / K
//       OARR(19) = TE - MOD(3000KM)     OARR(20) = TE(120KM) = TN = TI / K
//       OARR(21) = TI - MOD(430KM)      OARR(22) = X / KM, WHERE TE = TI
//       OARR(23) = SOL ZENITH ANG / DEG OARR(24) = SUN DECLINATION / DEG
//       OARR(25) = DIP / deg            OARR(26) = DIP LATITUDE / deg
//       OARR(27) = MODIFIED DIP LAT.OARR(28) = Geographic latitude
//       OARR(29) = sunrise / dec.hours OARR(30) = sunset / dec.hours
//       OARR(31) = ISEASON(1 = spring) OARR(32) = Geographic longitude
//      #OARR(33) = Rz12               OARR(34) = Covington Index
//      #OARR(35) = B1                 OARR(36) = M(3000)F2
//      $OARR(37) = TEC / m - 2           $OARR(38) = TEC_top / TEC * 100.
//      #OARR(39) = gind(IG12)        OARR(40) = F1 probability
//      #OARR(41) = F10.7 daily        OARR(42) = c1(F1 shape)
//       OARR(43) = daynr              OARR(44) = equatorial vertical
//       OARR(45) = foF2_storm / foF2_quiet         ion drift in m / s
//      #OARR(46) = F10.7_81           OARR(47) = foE_storm / foE_quiet
//       OARR(48) = spread - F probability
//       OARR(49) = Geomag.latitude   OARR(50) = Geomag.longitude
//       OARR(51) = ap at current time OARR(52) = daily ap
//       OARR(53) = invdip / degree      OARR(54) = MLT - Te
//       OARR(55) = CGM - latitude       OARR(56) = CGM - longitude
//       OARR(57) = CGM - MLT            OARR(58) = CGM lat eq.aurl bodry
//       OARR(59) = CGM - lati(MLT = 0)    OARR(60) = CGM - lati for MLT = 1
//       OARR(61) = CGM - lati(MLT = 2)    OARR(62) = CGM - lati for MLT = 3
//       OARR(63) = CGM - lati(MLT = 4)    OARR(64) = CGM - lati for MLT = 5
//       OARR(65) = CGM - lati(MLT = 6)    OARR(66) = CGM - lati for MLT = 7
//       OARR(67) = CGM - lati(MLT = 8)    OARR(68) = CGM - lati for MLT = 9
//       OARR(69) = CGM - lati(MLT = 10)   OARR(70) = CGM - lati for MLT = 11
//       OARR(71) = CGM - lati(MLT = 12)   OARR(72) = CGM - lati for MLT = 13
//       OARR(73) = CGM - lati(MLT = 14)   OARR(74) = CGM - lati for MLT = 15
//       OARR(75) = CGM - lati(MLT = 16)   OARR(76) = CGM - lati for MLT = 17
//       OARR(77) = CGM - lati(MLT = 18)   OARR(78) = CGM - lati for MLT = 19
//       OARR(79) = CGM - lati(MLT = 20)   OARR(80) = CGM - lati for MLT = 21
//       OARR(81) = CGM - lati(MLT = 22)   OARR(82) = CGM - lati for MLT = 23
//       OARR(83) = Kp at current time OARR(84) = magnetic declination
//       OARR(85) = L - value            OARR(86) = dipole moment
//                # INPUT as well as OUTPUT parameter
//                $ special for IRIWeb(only place - holders)
//---------------------------------------------------------------------- -
// * ****************************************************************
//*** THE ALTITUDE LIMITS ARE : LOWER(DAY / NIGHT)  UPPER        ***
//***     ELECTRON DENSITY         60 / 80 KM       1500 KM       ***
//***     TEMPERATURES               60 KM        2500 / 3000 KM  ***
//***     ION DENSITIES             100 KM        1500 KM       ***
//*****************************************************************
//*****************************************************************
//*********            INTERNALLY                    **************
//*********       ALL ANGLES ARE IN DEGREE           **************
//*********       ALL DENSITIES ARE IN M - 3 * *************
//*********       ALL ALTITUDES ARE IN KM            **************
//*********     ALL TEMPERATURES ARE IN KELVIN       **************
//*********     ALL TIMES ARE IN DECIMAL HOURS       **************
//*****************************************************************
//*****************************************************************
//*****************************************************************/
	unsigned int mess;
	//int icalls = 0;
	mess = jf[33];

	// set switches for NRLMSIS00
	int isw = 0;

	double *swmi = new double[25];
	for (size_t i = 0; i < 25; i++)
		swmi[i] = 1.0;

	int nummax = 1000;
	for (size_t i = 0; i < 20; i++)
		for (size_t j = 0; j < nummax; j++)
			outf[i][j] = -1.0;

	// oarr(1:6,15,16,33,39:41) is used for inputs
	for (size_t i = 6; i < 9; i++)
		oarr[i] = -1.0;
	for (size_t i = 10; i < 14; i++)
		oarr[i] = -1.0;
	for (size_t i = 16; i < 32; i++)
		oarr[i] = -1.0;
	oarr[33] = -1.0;
	for (size_t i = 35; i < 38; i++)
		oarr[i] = -1.0;
	oarr[39] = -1.0;
	for (size_t i = 41; i < 45; i++)
		oarr[i] = -1.0;
	for (size_t i = 46; i < 100; i++)
		oarr[i] = -1.0;

	int montho = -1,
		nmono = -1,
		iyearo = -1,
		idaynro = -1,
		rzino = 1,
		igino = 1,
		ut0 = -1,
		ursifo = 1;

	int numhei = int(abs(heiend - heibeg) / abs(heistp)) + 1;
	if (numhei > nummax)
		numhei = nummax;

	/*C NEW-GUL------------------------------
	c         Y05=.6931473
	c         QF=1.
	c	     h05top=0.
	C NEW-GUL------------------------------

	C
	C Code inserted to aleviate block data problem for PC version.
	C Thus avoiding DATA statement with parameters from COMMON block */

	double xdels[] = { 5., 5., 5., 10. },
		dnds[] = { 16e-3, 1e-2, 16e-3, 16e-3 },
		xnar[] = { 0., 0. },
		dte[] = { 5., 5., 10., 20., 20. },
		dti[] = { 10., 10., 20., 20. };

	int ddo[] = { 9, 5, 5, 25 },
		do2[] = { 5, 5 },
		tecon[] = { 0,0 };

	/*C
	C FIRST SPECIFY YOUR COMPUTERS CHANNEL NUMBERS ....................
	C AGNR = OUTPUT(OUTPUT IS DISPLAYED OR STORED IN FILE OUTPUT.IRI)...
	C IUCCIR = UNIT NUMBER FOR CCIR COEFFICIENTS ........................
	c*/
	const int IUCCIR = 10;
	/*c - web - special for web version
	c - web - messages should be turned off with mess = jf(34) = .false */
	const int KONSOL = 6;

	// selection of density, temperature and ion composition options
	int noden = !jf[0],
		notem = !jf[1],
		noion = !jf[2];
	if (!noion)
		noden = 0;
	int rbtt = !jf[5],
		layver = !jf[10],
		old79 = !jf[6],
		f1_ocpro = jf[18],
		f1_l_cond = !jf[19],
		dreg = jf[23];

	// rz12, IG12, F10.7D, PF10.7 input option
	int rzin = !jf[16];
	double arzin, aigin, f107din, f10781in;
	if (rzin)
		arzin = oarr[32];
	else
		oarr[32] = -1;

	int igin = !jf[26];
	if (igin)
		aigin = oarr[38];
	else
		oarr[38] = -1;

	if (!jf[24])
		f107din = oarr[40];
	else
		oarr[40] = -1;

	if (!jf[31])
		f10781in = oarr[45];
	else
		oarr[45] = -1;

	// Topside density
	int itopn;
	if (jf[28])
		itopn = 0;
	else
		if (jf[29])
			itopn = 1;
		else
			itopn = 2;

	double afof2, anmf2;
	// F2 peak density
	int fof2in = !jf[8];
	if (fof2in) {
		afof2 = oarr[0];
		anmf2 = oarr[0];
		if (oarr[0] < 100.)
			anmf2 = 1.24e10 * afof2 * afof2;
		else
			afof2 = sqrt(anmf2 / 1.24e10);
	}
	else
		oarr[0] = -1;
	int ursif2 = !jf[4];

	double ahmf2, afof1, anmf1;
	// F2 peak altitude
	int hmf2in = !jf[8];
	if (hmf2in)
		ahmf2 = oarr[1];
	else
		oarr[1] = -1;

	// F1 peak density
	int fof1in = !jf[12];
	if (fof1in) {
		afof1 = oarr[2];
		anmf1 = oarr[2];
		if (oarr[2] < 100)
			anmf1 = 1.24e10 * afof1 * afof1;
		else
			afof1 = sqrt(anmf1 / 1.24e10);
	}
	else
		oarr[2] = -1;

	// F1 peak altitude
	int hmf1in = !jf[13];
	if (hmf1in)
		double ahmf1 = oarr[3];
	else
		oarr[3] = -1;

	// E peak density
	int foein = !jf[14];
	double afoe, anme, ahme,
		b0_us, b1_us;
	if (foein) {
		afoe = oarr[4];
		anme = oarr[4];
		if (oarr[4] < 100)
			anme = 1.24e10 * afoe * afoe;
		else
			afoe = sqrt(anme / 1.24e10);
	}
	else
		oarr[4] = -1;

	// E peak altitude
	int hmein = !jf[15];
	if (hmein)
		ahme = oarr[5];
	else
		oarr[5] = -1;

	// B0 bottomside thickness
	int b0in = !jf[42];
	if (b0in)
		b0_us = oarr[9];
	else
		oarr[9] = -1;

	// B1 bottomside profile shape
	int b1in = !jf[43];
	if (jf[42])
		b1in = 0;
	if (b1in)
		b1_us = oarr[34];
	else
		oarr[34] = -1;

	// TE-NE MODEL OPTION
	int teneop = !jf[9];
	if (teneop)
		for (size_t i = 0; i < 2; i++) {
			xnar[i] = oarr[i + 14];
			tecon[i] = 0;
			if (xnar[i] > 0)
				tecon[i] = 1;
		}
	else {
		oarr[14] = -1;
		oarr[15] = -1;
	}

	/*C CALCULATION OF DAY OF YEAR OR MONTH/DAY AND DECIMAL YEAR 
	c NRDAYM is the number of days in the current month 
	c IDAYY is the number of days in the current year
	c
	c  leap year rule: years evenly divisible by 4 are leap years, except
	c  years also evenly divisible by 100 are not leap years, except years 
	c  also evenly divisible by 400 are leap years. The year 2000 is a 100 
	c  and 400 year exception and therefore it is a normal leap year. 
	c  The next 100 year exception will be in the year 2100!*/
	int iyear = iyyyy;
	if (iyear < 100)
		iyear += 1900;
	if (iyear < 30)
		iyear += 2000;
	int idayy = 365;
	if (iyear / 4 * 4 == iyear)
		idayy = 366; // leap year

	int month, iday, nrdaym, daynr, iyd;
	double amx, radj, height_center;
	if (mmdd < 0) {
		daynr = -mmdd;
		moda(1, iyear, month, iday, daynr, nrdaym);
	}
	else {
		month = mmdd / 100;
		iday = mmdd % 100;
		moda(0, iyear, month, iday, daynr, nrdaym);
	}
	findrlat.ryear = iyear + (daynr - 1.0) / idayy;
	iyd = iyear * 1000 + daynr;
	amx = constants::pi * (daynr - 3.) / 182.6;
	radj = 1. - constants::eexc * (cos(amx) + constants::eexc * (cos(2 * amx) - 1.) / 2.);

	// calculate center height for CGM computation
	height_center = (heibeg + heiend) / 2;

	/*C CALCULATION OF GEODETIC/GEOMAGNETIC COORDINATES (LATI, LONGI AND 
	C MLAT, MLONG), MAGNETIC INCLINATION (DIP), DIP LATITUDE (MAGBR) 
	C AND MODIFIED DIP (MODIP), ALL IN DEGREES*/
	double mlat, mlong, lati, longi,
		dec, dip, magbr, modip,
		xma, yma, zma, bet;
	if (along < 0)
		along = along + 360; // -180/180 to 0-360
	if (jmag > 0) {
		mlat = alati;
		mlong = along;
	}
	else {
		lati = alati;
		longi = along;
	}

	geodip(iyear, lati, longi, mlat, mlong, jmag);
	feldcof(findrlat.ryear);

	if (jf[17])
		igrf_dip(lati, longi, findrlat.ryear, 300.0, dec, dip, magbr, modip);
	else {
		fieldg(lati, longi, 300.0, xma, yma, zma, bet, dip, dec, modip);
		magbr = atan(0.5*tan(dip*constants::umr)) / constants::umr;
	}
	// calculate L - value, dip lati, and B_abs needed for invdip computation
	// calculating invdip at 600 km
	double invdip = -100.0,
		fl, dipl, babs, dimo,
		abslat, absmlt, absmdp,
		absmbr;
	int icode;
	if ((jf[1] && !jf[22]) || (jf[2] && !jf[5])) {
		igrf_sub(lati, longi, findrlat.ryear, 600.0, fl, icode, dipl, babs);
		if (fl > 10)
			fl = 10;
		invdip = invdpc(fl, dimo, babs, dipl);
	}

	abslat = abs(lati);
	absmlt = abs(mlat);
	absmdp = abs(modip);
	absmbr = abs(magbr);

	// CALCULATION OF UT/LT and XMLT
	double hourut, hour;
	if (dhour < 24) {
		hour = dhour; // dhour =< 24 is LT
		hourut = hour - longi / 15.;
		if (hourut < 0)
			hourut += 24;
	}
	else {
		hourut = dhour - 25;
		hour = hourut + longi / 15;
		if (hour > 24)
			hour -= 24;
	}

	double xmlt;
	clcmlt(iyear, daynr, hourut, lati, longi, xmlt);

	/*c SEASON assumes equal length seasons (92 days) with spring 
	c (SEASON=1) starting at day-of-year=45; for lati < 0 adjustment 
	c for southern hemisphere is made. Some models require the
	c seasonal month (ISEAMON) or the seasonal day-of year (SEADAY)
	c ZMONTH is decimal month (Jan 1 = 1.0 and Dec 31 = 12.97)
	c SDAY is the day number reduced to a 360 day year (TOPH05)
	c NRDAYM is the number of days in the current month 
	c IDAYY is the number of days in the current year*/
	int season = (daynr + 45.0) / 92.0;
	if (season < 1)
		season = 4;
	int nseasn = season;
	double zmonth = month + (iday - 1)*1. / nrdaym;
	double sday = daynr / idayy * 360;
	int seaday = daynr,
		iseamon = month;
	if (lati < 0.0) {
		season -= 2;
		if (season < 1)
			season += 4;
		iseamon = month + 6;
		if (iseamon > 12)
			iseamon -= 12;
		seaday = daynr + idayy / 2.;
		if (seaday > idayy)
			seaday -= idayy;
		sday += 180;
		if (sday > 360)
			sday -= 360;
	}
	/*C 12-month running mean sunspot number (rssn) and Ionospheric Global 
	C index (gind), daily F10.7 cm solar radio flux (f107d) and monthly 
	C F10.7 (cov) index*/
	unsigned int sam_mon = int(month == montho),
		sam_yea = int(iyear == iyearo),
		sam_doy = int(daynr == idaynro),
		sam_date = int(sam_yea && sam_doy),
		sam_ut = int(hourut == ut0);

	double *rzar = new double[3],
		*arig = new double[3],
		ttt, zi, rrr, rssn, gind,
		cov, f107d, f107y, f10781,
		f107365, f107_daily, f107pd,
		f107_81, f107_365,
		pf107, f_adj, f107yo,
		f10781o, covsat;
	int nmonth, iap_daily, isdate;
	if (!sam_date || rzino || rzin || igin || igino) {
		tcon(iyear, month, iday, daynr, rzar, arig, ttt, nmonth);
		if (nmonth < 0)
			return;
		if (rzin) {
			rrr = arzin;
			rzar[0] = rrr;
			rzar[1] = rrr;
			rzar[2] = rrr;
			if (!igin) {
				zi = (-0.0031*arzin + 1.5332)*arzin - 11.5634;
				arig[0] = zi;
				arig[1] = zi;
				arig[2] = zi;
			}
		}
		if (igin) {
			zi = aigin;
			arig[0] = zi;
			arig[1] = zi;
			arig[2] = zi;
			if (!rzin) {
				double xigin = aigin;
				if (xigin > 178.0066)
					xigin = 178.0066;
				rrr = 247.29 - 17.96*sqrt(178.0066 - xigin);
				rzar[0] = rrr;
				rzar[1] = rrr;
				rzar[2] = rrr;
			}
		}
		rssn = rzar[2];
		gind = arig[2];
		cov = 63.75 + rssn * (0.728 + rssn * 0.00089);
		/*C Getting F10.7 index: daily (f107d), previous day (f107y; 
		C required by MSIS), 81-day average (f10781), 365-day average
		C (f107365), and PF10.7=(F10.7_daily + F10.7_81_day)/2. 
		C F10.7 should be adjusted (to top of atmosphere) value not 
		C observed (at the ground) value */
		f107d = cov;
		f107y = cov;
		f10781 = cov;
		f107365 = cov;
		if (jf[31] || jf[24]) {
			apf_only(iyear, month, iday, f107_daily, f107pd, f107_81,
				f107_365, iap_daily, isdate);
			if (f107_daily > -11.1) {
				f107d = f107_daily;
				f107y = f107pd;
				f10781 = f107_81;
				f107365 = f107_365;
			}
		}
		if (!jf[24]) {
			f107d = f107din;
			f107y = f107din;
		}
		if (jf[31]) {
			f10781 = f10781in;
			f107365 = f10781in;
		}
		pf107 = (f107d + f10781) / 2.;
		f_adj = radj * radj;
		f107yo = f107y / f_adj;
		f10781o = f10781 / f_adj;

		if (jf[40])
			cov = f107365;
		covsat = cov;
		if (covsat > 188)
			covsat = 188;
	}
	double sundec, xhi1, sax80, sux80,
		sud1, xhi2, sax110, sux110,
		xhi, sax200, sux200, xhi3,
		sax300, sux300, sunde1, xhinon,
		sax1, sux1;
	unsigned int dnight = 0,
		enight = 0,
		fnight = 0;
	/*C CALCULATION OF SOLAR ZENITH ANGLE (XHI/DEG), SUN DECLINATION ANGLE 
	C (SUNDEC),SOLAR ZENITH ANGLE AT NOON (XHINON) AND TIME OF LOCAL 
	C SUNRISE/SUNSET (SAX, SUX; dec. hours) AT 70 KM (D-REGION), 110 KM
	C (E-REGION), 200 KM (F1-REGION), AND 500 KM (TE, TI) */
	soco(daynr, hour, lati, longi, 80., sundec, xhi1, sax80, sux80);
	soco(daynr, hour, lati, longi, 110., sud1, xhi2, sax110, sux110);
	soco(daynr, hour, lati, longi, 200., sud1, xhi, sax200, sux200);
	soco(daynr, hour, lati, longi, 300., sud1, xhi3, sax300, sux300);
	soco(daynr, 12.0, lati, longi, 110., sunde1, xhinon, sax1, sux1);

	if (abs(sax80) > 25) {
		if (sax80 < 0)
			dnight = 1;
	}
	else if (sax80 <= sux80) {
		if ((hour > sux80) || (hour < sax80))
			dnight = 1;
	}
	else if ((hour > sux80) && (hour < sax80)) {
		dnight = 1;
	}

	if (abs(sax110) > 25) {
		if (sax110 < 0)
			enight = 1;
	}
	else if (sax110 <= sux110) {
		if ((hour > sux110) || (hour < sax110))
			enight = 1;
	}
	else if ((hour > sux110) && (hour < sax110)) {
		enight = 1;
	}

	if (abs(sax200) > 25) {
		if (sax200 < 0)
			fnight = 1;
	}
	else if (sax200 <= sux200) {
		if ((hour > sux200) || (hour < sax200))
			fnight = 1;
	}
	else if ((hour > sux200) && (hour < sax200)) {
		fnight = 1;
	}

	/*C CALCULATION OF ELECTRON DENSITY PARAMETERS................
	C lower height boundary (HNEA), upper boundary (HNEE)*/
	double hnea = 65.,
		hnee = 2000.,
		dela = 4.32,
		dell, foe, nme, hme;
	if (dnight)
		hnea = 80;
	std::string filename;
	if (!noden) {
		if (absmdp >= 18)
			dela = 1.0 + exp(-(absmdp - 30.0) / 10.0);
		dell = 1 + exp(-(abslat - 20.) / 10.);
		// E peak critical frequency (foE), density (NmE), and height (hmE)
		if (foein) {
			foe = afoe;
			nme = anme;
		}
		else {
			foe = foeedi(cov, xhi, xhinon, abslat);
			nme = 1.24e10 * foe * foe;
		}
		if (hmein)
			hme = ahme;
		else
			hme = 110.;

		/*c F2 peak critical frequency foF2, density NmF2, and height hmF2
		c
		C READ CCIR AND URSI COEFFICIENT SET FOR CHOSEN MONTH*/
		bool f4292 = false, f4293 = false;
		double ***f2 = new double**[13];
		for (size_t i = 0; i < 13; i++)	{
			f2[i] = new double*[76];
			for (size_t j = 0; j < 76; j++)
				f2[i][j] = new double[2];
		}
		double ***fm3 = new double**[9];
		for (size_t i = 0; i < 9; i++) {
			f2[i] = new double*[49];
			for (size_t j = 0; j < 49; j++)
				f2[i][j] = new double[2];
		}
		if (!(fof2in && hmf2in && (itopn != 2))) {
			if (ursif2 == ursifo) {
				if (rzin || rzino || igin || igino) {
					if (sam_mon && (nmonth == nmono) && sam_yea)
						f4292 = true;
					else if (sam_mon)
						f4293 = true;
				}
			} // 7797
			if (!f4292) {
				ursifo = ursif2;
				filename = "ccir" + std::to_string(month + 10) + ".asc";

				std::fstream myfile("data/" + filename, std::ios_base::in);

				if (myfile.is_open())
				{
					for (size_t k = 0; k < 2; k++)
						for (size_t j = 0; j < 76; j++)
							for (size_t i = 0; i < 13; i++)
								myfile >> f2[i][j][k];

					for (size_t k = 0; k < 2; k++)
						for (size_t j = 0; j < 49; j++)
							for (size_t i = 0; i < 9; i++)
								myfile >> fm3[i][j][k];

				
					myfile.close();
				}
				else
					std::cout << "ERROR. FILE WASN'T OPENED!" << std::endl;
			}

		} // 501
	} // 4933
}
// TODO:
void iri_tec(double hstart, double hend, int istep,
	double &tectot, double &tectop, double &tecbot) {
/*C subroutine to compute the total ionospheric content
C INPUT:      
C   hstart  altitude (in km) where integration should start
C   hend    altitude (in km) where integration should end
C   istep   =0 [fast, but higher uncertainty <5%]
C           =1 [standard, recommended]
C           =2 [stepsize of 1 km; best TEC, longest CPU time]
C OUTPUT:
C   tectot  total ionospheric content in tec-units (10^16 m^-2)
C   tectop  topside content (in %)
C   tecbot  bottomside content (in %)
C
C The different stepsizes for the numerical integration are 
c defined as follows (h1=100km, h2=hmF2-10km, h3=hmF2+10km, 
c h4=hmF2+150km, h5=hmF2+250km):
C       istep   h1-h2   h2-h3   h3-h4   h4-h5   h5-hend
C       0       2.0km   1.0km   2.5km   exponential approximation
C       1       2.0km   1.0km   2.5km   10.0km  30.0km
C       2       1.0km   0.5km   1.0km   1.0km   1.0km*/
	bool expo = false;
	int numstep = 5;
	double xnorm = xnmf2 / 1000.;

}

void moda(int in, int iyear, int &month, int &iday, int &idoy, int &nrdaymo) {
/*C CALCULATES DAY OF YEAR (IDOY, ddd) FROM YEAR (IYEAR, yy or yyyy), 
C MONTH (MONTH, mm) AND DAY OF MONTH (IDAY, dd) IF IN=0, OR MONTH 
C AND DAY FROM YEAR AND DAY OF YEAR IF IN=1. NRDAYMO is an output 
C parameter providing the number of days in the specific month */
	int mm[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	int imo = 0, mosum,
		mobe = 0, moold;
/*c  leap year rule: years evenly divisible by 4 are leap years, except
c  years also evenly divisible by 100 are not leap years, except years 
c  also evenly divisible by 400 are leap years. The year 2000 therefore 
C  is a leap year. The 100 and 400 year exception rule
c     if((iyear/4*4.eq.iyear).and.(iyear/100*100.ne.iyear)) mm(2)=29
c  will become important again in the year 2100 which is not a leap 
C  year */
	mm[1] = 28;
	if (iyear / 4 * 4 == iyear) mm[1] = 29;

	if (in == 0) {
		mosum = 0;
		if (month > 1)
			for (size_t i = 0; i < month - 1; i++)
				mosum += mm[i];
		idoy = mosum + iday;
		nrdaymo = mm[month-1];
	}
	else {
		do {
			imo += 1;
			if (imo > 12)
				break;
			moold = mobe;
			nrdaymo = mm[imo - 1];
			mobe += nrdaymo;
		} while (mobe < idoy);
		month = imo;
		iday = idoy - moold;
	}
}

void geodip(int &iyr, double &sla, double &slo, double &dla, double &dlo, int j) {
/*C  Calculates dipole geomagnetic coordinates from geocentric coordinates
C  or vice versa.

C                     J=0           J=1
C		INPUT:     J,SLA,SLO     J,DLA,DLO
C		OUTPUT:     DLA,DLO       SLA,SLO

C  Last revision: November 2005 (Vladimir Papitashvili)
C  The code is modifed from GEOCOR written by V.Popov and V.Papitashvili
C  in mid-1980s */
	// Earth's radius (km) RE = 6371.2
	// The radius of the sphere to compute the coordinates (in Re)
	//        RH = (RE + HI) / RE
	double r = 1.;
	double x, y, z, xm, ym, zm, rm, th, pf, col, rlo;
	if (j == 0) {
		col = (90. - sla) * constants::umr;
		rlo = slo * constants::umr;
		sphcar(r, col, rlo, x, y, z, 1);
		geomag(x, y, z, xm, ym, zm, 1, iyr);
		sphcar(rm, th, pf, xm, ym, zm, -1);
		dlo = pf / constants::umr;
		dla = 90. - th / constants::umr;
	}
	else {
		col = (90. - dla) * constants::umr;
		rlo = dlo * constants::umr;
		sphcar(r, col, rlo, xm, ym, zm, 1);
		geomag(x, y, z, xm, ym, zm, -1, iyr);
		sphcar(rm, th, pf, x, y, z, -1);
		slo = pf / constants::umr;
		sla = 90. - th / constants::umr;
	}
}

void sphcar(double &r, double &teta, double &phi,
	double &x, double &y, double &z, int j) {
/*C  *********************************************************************
C   CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICA VERSA
C    (TETA AND PHI IN RADIANS).
C                  J>0            J<0
C-----INPUT:   J,R,TETA,PHI     J,X,Y,Z
C----OUTPUT:      X,Y,Z        R,TETA,PHI
C  AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
C      STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
C      (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)
C  *********************************************************************/
	double sq;
	if (j == 0) {
		sq = x * x + y * y;
		r = sqrt(sq + z * z);
		if (sq != 0.) {
			sq = sqrt(sq);
			phi = atan2(y, x);
			teta = atan2(sq, z);
			if (phi < 0)
				phi += 2 * constants::pi;
		}
		else {
			phi = 0.;
			if (z < 0)
				teta = constants::pi;
			else
				teta = 0;
		}
	}
	else {
		sq = r * sin(teta);
		x = sq * cos(phi);
		y = sq * sin(phi);
		z = r * cos(teta);
	}
}

void geomag(double &xgeo, double &ygeo, double &zgeo, double &xmag, 
	double &ymag, double &zmag, int j, int iyr) {
/*C  *********************************************************************
C CONVERTS GEOCENTRIC (GEO) TO DIPOLE (MAG) COORDINATES OR VICA VERSA.
C IYR IS YEAR NUMBER (FOUR DIGITS).

C                           J>0                J<0
C-----INPUT:  J,XGEO,YGEO,ZGEO,IYR   J,XMAG,YMAG,ZMAG,IYR
C-----OUTPUT:    XMAG,YMAG,ZMAG        XGEO,YGEO,ZGEO

C  AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
C      STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
C      (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)
C  **********************************************************************/
	int ii = 1;

	if (iyr != 1) {
		ii = iyr;
		recalc(ii, 0, 25, 0, 0);
	}
	if (j >= 0) {
		xmag = xgeo * c1.ctcl + ygeo * c1.ctsl - zgeo * c1.st0;
		ymag = ygeo * c1.cl0 - xgeo * c1.sl0;
		zmag = xgeo * c1.stcl + ygeo * c1.stsl + zgeo * c1.ct0;
	}
	else {
		xgeo = xmag * c1.ctcl - ymag * c1.sl0 + zmag * c1.stcl;
		ygeo = xmag * c1.ctsl + ymag * c1.cl0 + zmag * c1.stsl;
		zgeo = zmag * c1.ct0 - xmag * c1.st0;
	}
}

void recalc(int iyr, int iday, int ihour, int min, int isec) {
	/*C  *********************************************************************
	C  If only IYR is given then CALL RECALC(IYR,0,25,0,0)
	C  THIS IS A MODIFIED VERSION OF THE SUBROUTINE RECOMP WRITTEN BY
	C  N. A. TSYGANENKO. SINCE I WANT TO USE IT IN PLACE OF SUBROUTINE
	C  RECALC, I HAVE RENAMED THIS ROUTINE RECALC AND ELIMINATED THE
	C  ORIGINAL RECALC FROM THIS VERSION OF THE <GEOPACK.FOR> PACKAGE.
	C  THIS WAY ALL ORIGINAL CALLS TO RECALC WILL CONTINUE TO WORK WITHOUT
	C  HAVING TO CHANGE THEM TO CALLS TO RECOMP.
	C
	C  AN ALTERNATIVE VERSION OF THE SUBROUTINE RECALC FROM THE GEOPACK
	C  PACKAGE BASED ON A DIFFERENT APPROACH TO DERIVATION OF ROTATION
	C  MATRIX ELEMENTS
	C
	C  THIS SUBROUTINE WORKS BY 20% FASTER THAN RECALC AND IS EASIER TO
	C  UNDERSTAND
	C  #####################################################
	C  #  WRITTEN BY  N.A. TSYGANENKO ON DECEMBER 1, 1991  #
	C  #####################################################
	C  Modified by Mauricio Peredo, Hughes STX at NASA/GSFC Code 695,
	C  September 1992
	C
	C  Modified to accept years up to year 2000 and updated IGRF coeficients
	C     from 1945 (updated by V. Papitashvili, February 1995)
	C
	C  Modified to accept years up to 2005 (V. Papitashvili, January 2001)
	C
	C  Modified to accept years from 1900 through 2010 using the DGRF &
	C     IGRF-10 coeficients (updated by V. Papitashvili, November 2005)
	C
	C  Modified to accept years up to 2015 (V. Papitashvili, January 2011)
	C
	C  Modified to accept years up to 2020 (D. Bilitza, October 2015)
	C
	C   OTHER SUBROUTINES CALLED BY THIS ONE: SUN
	C
	C     IYR = YEAR NUMBER (FOUR DIGITS)
	C     IDAY = DAY OF YEAR (DAY 1 = JAN 1)
	C     IHOUR = HOUR OF DAY (00 TO 23)
	C     MIN = MINUTE OF HOUR (00 TO 59)
	C     ISEC = SECONDS OF DAY(00 TO 59)
	C  **********************************************************************/
	int iye = 0, ide = 0, iy;
	double f1, f2, g10, g11, h11, dt;
	if ((iyr != iye) || (iday != ide)) {
		ide = iday;

		if (iyr < 1900)
			iy = 1900;
		else if (iyr > 2020)
			iy = 2020;
		else
			iy = iyr;

		/*C  LINEAR INTERPOLATION OF THE GEODIPOLE MOMENT COMPONENTS BETWEEN THE
		C  VALUES FOR THE NEAREST EPOCHS:*/

		if (iy < 1905) {
			f2 = (iy + 1.0 * iday / 365 - 1900) / 5;
			f1 = 1 - f2;
			g10 = 31543 * f1 + 31464 * f2;
			g11 = -2298 * f1 - 2298 * f2;
			h11 = 5922 * f1 + 5909 * f2;
		}
		else if (iy < 1910) {
			f2 = (iy + 1.0 * iday / 365 - 1905) / 5;
			f1 = 1 - f2;
			g10 = 31464 * f1 + 31354 * f2;
			g11 = -2298 * f1 - 2297 * f2;
			h11 = 5909 * f1 + 5898 * f2;
		}
		else if (iy < 1915) {
			f2 = (iy + 1.0 * iday / 365 - 1910) / 5;
			f1 = 1 - f2;
			g10 = 31354 * f1 + 31212 * f2;
			g11 = -2297 * f1 - 2306 * f2;
			h11 = 5898 * f1 + 5875 * f2;
		}
		else if (iy < 1920) {
			f2 = (iy + 1.0 * iday / 365 - 1915) / 5;
			f1 = 1 - f2;
			g10 = 31212 * f1 + 31060 * f2;
			g11 = -2306 * f1 - 2317 * f2;
			h11 = 5875 * f1 + 5845 * f2;
		}
		else if (iy < 1925) {
			f2 = (iy + 1.0 * iday / 365 - 1920) / 5;
			f1 = 1 - f2;
			g10 = 31060 * f1 + 30926 * f2;
			g11 = -2317 * f1 - 2318 * f2;
			h11 = 5845 * f1 + 5817 * f2;
		}
		else if (iy < 1930) {
			f2 = (iy + 1.0 * iday / 365 - 1925) / 5;
			f1 = 1 - f2;
			g10 = 30926 * f1 + 30805 * f2;
			g11 = -2318 * f1 - 2316 * f2;
			h11 = 5817 * f1 + 5808 * f2;
		}
		else if (iy < 1935) {
			f2 = (iy + 1.0 * iday / 365 - 1930) / 5;
			f1 = 1 - f2;
			g10 = 30805 * f1 + 30715 * f2;
			g11 = -2316 * f1 - 2306 * f2;
			h11 = 5808 * f1 + 5812 * f2;
		}
		else if (iy < 1940) {
			f2 = (iy + 1.0 * iday / 365 - 1935) / 5;
			f1 = 1 - f2;
			g10 = 30715 * f1 + 30654 * f2;
			g11 = -2306 * f1 - 2292 * f2;
			h11 = 5812 * f1 + 5821 * f2;
		}
		else if (iy < 1945) {
			f2 = (iy + 1.0 * iday / 365 - 1940) / 5;
			f1 = 1 - f2;
			g10 = 30654 * f1 + 30594 * f2;
			g11 = -2292 * f1 - 2285 * f2;
			h11 = 5821 * f1 + 5810 * f2;
		}
		else if (iy < 1950) {
			f2 = (iy + 1.0 * iday / 365 - 1945) / 5;
			f1 = 1 - f2;
			g10 = 30594 * f1 + 30554 * f2;
			g11 = -2285 * f1 - 2250 * f2;
			h11 = 5810 * f1 + 5815 * f2;
		}
		else if (iy < 1955) {
			f2 = (iy + 1.0 * iday / 365 - 1950) / 5;
			f1 = 1 - f2;
			g10 = 30554 * f1 + 30500 * f2;
			g11 = -2250 * f1 - 2215 * f2;
			h11 = 5815 * f1 + 5820 * f2;
		}
		else if (iy < 1960) {
			f2 = (iy + 1.0 * iday / 365 - 1955) / 5;
			f1 = 1 - f2;
			g10 = 30500 * f1 + 30421 * f2;
			g11 = -2215 * f1 - 2169 * f2;
			h11 = 5820 * f1 + 5791 * f2;
		}
		else if (iy < 1965) {
			f2 = (iy + 1.0 * iday / 365 - 1960) / 5;
			f1 = 1 - f2;
			g10 = 30421 * f1 + 30334 * f2;
			g11 = -2169 * f1 - 2119 * f2;
			h11 = 5791 * f1 + 5776 * f2;
		}
		else if (iy < 1970) {
			f2 = (iy + 1.0 * iday / 365 - 1965) / 5;
			f1 = 1 - f2;
			g10 = 30334 * f1 + 30220 * f2;
			g11 = -2119 * f1 - 2068 * f2;
			h11 = 5776 * f1 + 5737 * f2;
		}
		else if (iy < 1975) {
			f2 = (iy + 1.0 * iday / 365 - 1970) / 5;
			f1 = 1 - f2;
			g10 = 30220 * f1 + 30100 * f2;
			g11 = -2068 * f1 - 2013 * f2;
			h11 = 5737 * f1 + 5675 * f2;
		}
		else if (iy < 1980) {
			f2 = (iy + 1.0 * iday / 365 - 1975) / 5;
			f1 = 1 - f2;
			g10 = 30100 * f1 + 29992 * f2;
			g11 = -2013 * f1 - 1956 * f2;
			h11 = 5675 * f1 + 5604 * f2;
		}
		else if (iy < 1985) {
			f2 = (iy + 1.0 * iday / 365 - 1980) / 5;
			f1 = 1 - f2;
			g10 = 29992 * f1 + 29873 * f2;
			g11 = -1956 * f1 - 1905 * f2;
			h11 = 5604 * f1 + 5500 * f2;
		}
		else if (iy < 1990) {
			f2 = (iy + 1.0 * iday / 365 - 1985) / 5;
			f1 = 1 - f2;
			g10 = 29873 * f1 + 29775 * f2;
			g11 = -1905 * f1 - 1848 * f2;
			h11 = 5500 * f1 + 5406 * f2;
		}
		else if (iy < 1995) {
			f2 = (iy + 1.0 * iday / 365 - 1990) / 5;
			f1 = 1 - f2;
			g10 = 29775 * f1 + 29692 * f2;
			g11 = -1848 * f1 - 1784 * f2;
			h11 = 5406 * f1 + 5306 * f2;
		}
		else if (iy < 2000) {
			f2 = (iy + 1.0 * iday / 365 - 1995) / 5;
			f1 = 1 - f2;
			g10 = 29692 * f1 + 29619.4*f2;
			g11 = -1784 * f1 - 1728.2*f2;
			h11 = 5306 * f1 + 5186.1*f2;
		}
		else if (iy < 2005) {
			f2 = (iy + 1.0 * iday / 365 - 2000) / 5;
			f1 = 1 - f2;
			g10 = 29619.4 * f1 + 29554.63 * f2;
			g11 = -1728.2 * f1 - 1669.05 * f2;
			h11 = 5186.1 * f1 + 5077.99 * f2;
		}
		else if (iy < 2010) {
			f2 = (iy + 1.0 * iday / 365 - 2005) / 5;
			f1 = 1 - f2;
			g10 = 29554.63 * f1 + 29496.57 * f2;
			g11 = -1669.05 * f1 - 1586.42 * f2;
			h11 = 5077.99 * f1 + 4944.26 * f2;
		}
		else if (iy < 2015) {
			f2 = (iy + 1.0 * iday / 365 - 2010) / 5;
			f1 = 1 - f2;
			g10 = 29496.57 * f1 + 29442.0 * f2;
			g11 = -1586.42 * f1 - 1501.0 * f2;
			h11 = 4944.26 * f1 + 4797.1 * f2;
		}
		else {
			dt = iy + 1.0 * iday / 365 - 2015;
			g10 = 29442.0 - 10.3 * dt;
			g11 = -1501.0 + 18.1 * dt;
			h11 = 4797.1 - 26.6 * dt;
		}
		/*C  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EzMAG IN GEO COORD
		C  SYSTEM:
		C  SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
		C         ST0 * CL0                ST0 * SL0                CT0*/
		double sq, sqq, sqr;
		sq = g11 * g11 + h11 * h11;
		sqq = sqrt(sq);
		sqr = sqrt(g10 * g10 + sq);
		c1.sl0 = -h11 / sqq;
		c1.cl0 = -g11 / sqq;
		c1.st0 = sqq / sqr;
		c1.ct0 = g10 / sqr;
		c1.stcl = c1.st0 * c1.cl0;
		c1.stsl = c1.st0 * c1.sl0;
		c1.ctsl = c1.ct0 * c1.sl0;
		c1.ctcl = c1.ct0 * c1.cl0;
	}
	/*C  THE CALCULATIONS ARE TERMINATED IF ONLY GEO-MAG TRANSFORMATION
	C  IS TO BE DONE  (IHOUR>24 IS THE AGREED CONDITION FOR THIS CASE):*/
	if (ihour < 25)
		std::cout << "\nError! recalc must be continued" << std::endl;
}

void feldcof(double year) {
/*c-----------------------------------------------------------------------        
C  DETERMINES COEFFICIENTS AND DIPOL MOMENT FROM IGRF MODELS
C
C       INPUT:  YEAR    DECIMAL YEAR FOR WHICH GEOMAGNETIC FIELD IS TO
C                       BE CALCULATED
C        				COMMON/IGRF1/ERAD,AQUAD,BQUAD,DIMO /CONST/UMR,PI
C       OUTPUT:         COMMON/MODEL/NMAX,TIME,GH1,FIL1
C        				COMMON/DIPOL/GHI1,GHI2,GHI3
C
C THE GEOMAGNETIC DIPOL MOMENT (DIMO) IN GAUSS (NORMALIZED TO EARTH'S 
C RADIUS) AT THE TIME (YEAR) IS COMPUTED BUT NOT USED.
C
C 05/31/2000 updated to IGRF-2000 version (###) 
C 03/24/2000 updated to IGRF-2005 version (###) 
C 07/22/2009 NMAX=13 for DGRF00 and IGRF05; H/G-arrays(195)
C 02/26/2010 update to IGRF-11 (2010) (###)  
C 10/05/2011 added COMMON/DIPOL/ for MLT computation in DPMTRX (IRIFUN)
C 02/10/2015 update to IGRF-12 (2015) (###)
c-----------------------------------------------------------------------*/
	// C ### numye is number of IGRF coefficient files minus 1
	int numye = 15;
	/*C  IS=0 FOR SCHMIDT NORMALIZATION   IS=1 GAUSS NORMALIZATION
	C  IU  IS INPUT UNIT NUMBER FOR IGRF COEFFICIENT SETS*/
	int iu = 14, is = 0;
	//C-- DETERMINE IGRF-YEARS FOR INPUT-YEAR
	double time = year;
	int iyea = int(year / 5) * 5;
	int l = (iyea - 1945) / 5 + 1;
	if (l < 1) l = 1;
	else if (l > numye) l = numye;
	
	std::string filmod[] = { "dgrf1945.dat", "dgrf1950.dat", "dgrf1955.dat",
		"dgrf1960.dat", "dgrf1965.dat", "dgrf1970.dat", "dgrf1975.dat",
		"dgrf1980.dat", "dgrf1985.dat", "dgrf1990.dat", "dgrf1995.dat",
		"dgrf2000.dat", "dgrf2005.dat", "dgrf2010.dat", "igrf2015.dat",
		"igrf2015s.dat" };
	double dtemod[] = { 945., 1950., 1955., 1960., 1965.,
		1970., 1975., 1980., 1985., 1990., 1995., 2000.,2005.,
		2010., 2015., 2020. };

	double dte1 = dtemod[l - 1];
	std::string fil1 = "data/" + filmod[l - 1];
	double dte2 = dtemod[l];
	std::string fil2 = "data/" + filmod[l];

	// C-- GET IGRF COEFFICIENTS FOR THE BOUNDARY YEARS
	int nmax1, nmax2, ier;
	double erad, *gh1 = new double[196],
		*gh2 = new double[196],
		*gha = new double[196];
	getshc(fil1, nmax1, erad, gh1);
	getshc(fil2, nmax2, erad, gh2);
	// C-- DETERMINE IGRF COEFFICIENTS FOR YEAR
	if (l <= numye - 1)
		intershc(year, dte1, nmax1, gh1, dte2, nmax2, gh2, model.nmax, gha);
	else
		extrashc(year, dte1, nmax1, gh1, nmax2, gh2, model.nmax, gha);

	// C-- DETERMINE MAGNETIC DIPOL MOMENT AND COEFFIECIENTS G
	double f0 = 0, f;
	for (size_t i = 0; i < 3; i++)	{
		f = gha[i] * 1e-5;
		f0 += f * f;
	}
	double dimo = sqrt(f0);
	dipol.ghi1 = gha[0];
	dipol.ghi2 = gha[1];
	dipol.ghi3 = gha[2];
	gh1[0] = 0;
	int i = 1;
	f0 = 1e-5;
	if (is == 0) f0 = -f0;

	for (size_t n = 1; n <= model.nmax; n++) {
		f0 = f0 * n*n / (4 * n - 2);
		if (is == 0) f0 *= (2 * n - 1) / n;
		f = f0 * 0.5;
		if (is == 0) f *= sqrt(2);
		gh1[i] = gha[i - 1] * f0;
		i++;
		for (size_t m = 1; m <= n; m++) {
			f *= (n + m) / (n - m + 1);
			if (is == 0) f *= sqrt((n - m + 1) / (n + m));
			gh1[i] = gha[i - 1] * f;
			gh1[i + 1] = gha[i] * f;
		}
	}
}

void getshc(std::string fspec, int &nmax, double &erad, double *gh) {
/*C ===============================================================               
C       Reads spherical harmonic coefficients from the specified     
C       file into an array.                                          
C       Input:                                                       
C           IU    - Logical unit number                              
C           FSPEC - File specification                               
C       Output:                                                      
C           NMAX  - Maximum degree and order of model                
C           ERAD  - Earth's radius associated with the spherical     
C                   harmonic coefficients, in the same units as      
C                   elevation                                        
C           GH    - Schmidt quasi-normal internal spherical          
C                   harmonic coefficients                            
C           IER   - Error number: =  0, no error                     
C                                 = -2, records out of order         
C                                 = FORTRAN run-time error number    
C ===============================================================*/
	for (size_t i = 0; i < 196; i++)
		gh[i] = 0;

	std::string fout;
	double xmyear;
	std::fstream myfile(fspec, std::ios_base::in);

	if (myfile.is_open())
	{
		myfile >> fout >> nmax >> erad >> xmyear;
		int nm = nmax * (nmax + 2);

		for (size_t i = 0; i < nm; i++)
			myfile >> gh[i];

		myfile.close();
	}
	else
		std::cout << "ERROR. FILE " << fspec << " WASN'T OPENED!" << std::endl;
}

void intershc(double date, double dte1, int nmax1, double *gh1, double dte2, int nmax2,
	double *gh2, int &nmax, double *gh) {
/*C ===============================================================               
C                                                                               
C       Version 1.01                                                 
C                                                                               
C       Interpolates linearly, in time, between two spherical        
C       harmonic models.                                             
C                                                                               
C       Input:                                                       
C           DATE  - Date of resulting model (in decimal year)        
C           DTE1  - Date of earlier model                            
C           NMAX1 - Maximum degree and order of earlier model        
C           GH1   - Schmidt quasi-normal internal spherical          
C                   harmonic coefficients of earlier model           
C           DTE2  - Date of later model                              
C           NMAX2 - Maximum degree and order of later model          
C           GH2   - Schmidt quasi-normal internal spherical          
C                   harmonic coefficients of later model             
C                                                                               
C       Output:                                                      
C           GH    - Coefficients of resulting model                  
C           NMAX  - Maximum degree and order of resulting model      
C                                                                               
C       A. Zunde                                                     
C       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225    
C                                                                               
C ===============================================================  */
/*C ---------------------------------------------------------------               
C       The coefficients (GH) of the resulting model, at date        
C       DATE, are computed by linearly interpolating between the     
C       coefficients of the earlier model (GH1), at date DTE1,       
C       and those of the later model (GH2), at date DTE2. If one     
C       model is smaller than the other, the interpolation is        
C       performed with the missing coefficients assumed to be 0.     
C --------------------------------------------------------------- */
	double factor = (date - dte1) / (dte2 - dte1);
	int k, l;
	if (nmax1 == nmax2) {
		k = nmax1 * (nmax1 + 2);
		nmax = nmax1;
	}
	else if (nmax1 > nmax2) {
		k = nmax2 * (nmax2 + 2);
		l = nmax1 * (nmax1 + 2);
		for (size_t i = k; i < l; i++)
			gh[i] = gh1[i] - factor * gh1[i];
		nmax = nmax1;
	}
	else {
		k = nmax1 * (nmax1 + 2);
		l = nmax2 * (nmax2 + 2);
		for (size_t i = k; i < l; i++)
			gh[i] = factor * gh2[i];
		nmax = nmax2;
	}

	for (size_t i = 0; i < k; i++)
		gh[i] = gh1[i] + factor * (gh2[i] - gh1[i]);
}

void extrashc(double date, double dte1, int nmax1, double *gh1, int nmax2,
	double *gh2, int &nmax, double *gh) {
/*C ===============================================================               
C                                                                               
C       Version 1.01                                                   
C                                                                               
C       Extrapolates linearly a spherical harmonic model with a        
C       rate-of-change model.                                          
C                                                                               
C       Input:                                                         
C           DATE  - Date of resulting model (in decimal year)          
C           DTE1  - Date of base model                                 
C           NMAX1 - Maximum degree and order of base model             
C           GH1   - Schmidt quasi-normal internal spherical            
C                   harmonic coefficients of base model                
C           NMAX2 - Maximum degree and order of rate-of-change         
C                   model                                              
C           GH2   - Schmidt quasi-normal internal spherical            
C                   harmonic coefficients of rate-of-change model      
C                                                                               
C       Output:                                                        
C           GH    - Coefficients of resulting model                    
C           NMAX  - Maximum degree and order of resulting model        
C                                                                               
C       A. Zunde                                                       
C       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225      
C                                                                               
C ===============================================================               
                                                                                
        DIMENSION       GH1(*), GH2(*), GH(*)                        
                                                                                
C ---------------------------------------------------------------               
C       The coefficients (GH) of the resulting model, at date          
C       DATE, are computed by linearly extrapolating the coef-         
C       ficients of the base model (GH1), at date DTE1, using          
C       those of the rate-of-change model (GH2), at date DTE2. If      
C       one model is smaller than the other, the extrapolation is      
C       performed with the missing coefficients assumed to be 0.       
C ---------------------------------------------------------------*/
	double factor = (date - dte1);
	int k, l;
	if (nmax1 == nmax2) {
		k = nmax1 * (nmax1 + 2);
		nmax = nmax1;
	}
	else if (nmax1 > nmax2) {
		k = nmax2 * (nmax2 + 2);
		l = nmax1 * (nmax1 + 2);
		for (size_t i = k; i < l; i++)
			gh[i] = gh1[i];
		nmax = nmax1;
	}
	else {
		k = nmax1 * (nmax1 + 2);
		l = nmax2 * (nmax2 + 2);
		for (size_t i = k; i < l; i++)
			gh[i] = factor * gh2[i];
		nmax = nmax2;
	}
	for (size_t i = 0; i < k; i++)
		gh[i] = gh1[i] + factor * gh2[i];
}

void igrf_dip(double xlat, double xlong, double year, double height,
	double &dec, double &dip, double &dipl, double &ymodip) {
/*c-----------------------------------------------------------------------        
c INPUT:
c    xlat      geodatic latitude in degrees
c    xlong     geodatic longitude in degrees
c    year      decimal year (year+month/12.0-0.5 or 
c                  year+day-of-year/365 or ../366 if leap year) 
c    height    height in km
c OUTPUT:
c    dec       magnetic declination in degrees
c    dip       magnetic inclination (dip) in degrees
c    dipl      dip latitude in degrees
c    ymodip    modified dip latitude = asin{dip/sqrt[dip^2+cos(LATI)]} 
c-----------------------------------------------------------------------*/
	double xlati = xlat,
		xlongi = xlong,
		h = height,
		bnorth, beast, bdown, babs;

	feldg(xlati, xlongi, h, bnorth, beast, bdown, babs);
	double decarg = beast / sqrt(beast * beast + bnorth * bnorth);
	if (abs(decarg) > 1.) decarg = decarg / abs(decarg);
	dec = asin(decarg);
	double bdba = bdown / babs;
	if (abs(bdba) > 1.) bdba = bdba / abs(bdba);
	dip = asin(bdba);
	double dipdiv = dip / sqrt(dip * dip + cos(xlati * constants::umr));
	if (abs(dipdiv) > 1.) dipdiv = dipdiv / abs(dipdiv);
	double smodip = asin(dipdiv);
	dipl = atan(bdown / 2.0 / sqrt(bnorth * bnorth + beast * beast)) /
		constants::umr;
	ymodip = smodip / constants::umr;
	dec /= constants::umr;
	dip /= constants::umr;
}

void feldg(double glat, double glon, double alt, 
	double &bnorth, double &beast, double &bdown, double &babs) {
/*c-----------------------------------------------------------------------        
C CALCULATES EARTH MAGNETIC FIELD FROM SPHERICAL HARMONICS MODEL
C REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTRE, INTERNAL NOTE 61, 
C      1970.
c-----------------------------------------------------------------------        
C CHANGES (D. BILITZA, NOV 87):
C   - FIELD COEFFICIENTS IN BINARY DATA FILES INSTEAD OF BLOCK DATA
C   - CALCULATES DIPOL MOMENT
C 09/07/22 NMAX=13 for DGRF00 and IGRF05; H/G-arrays(195)
c-----------------------------------------------------------------------        
C  INPUT:  ENTRY POINT FELDG
C               GLAT  GEODETIC LATITUDE IN DEGREES (NORTH)
C               GLON  GEODETIC LONGITUDE IN DEGREES (EAST)
C               ALT   ALTITUDE IN KM ABOVE SEA LEVEL
C
C          ENTRY POINT FELDC
C               V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM)
C                       X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE
C                       Y-AXIS POINTING TO EQUATOR AT 90 LONG.
C                       Z-AXIS POINTING TO NORTH POLE
C
C          COMMON BLANK AND ENTRY POINT FELDI ARE NEEDED WHEN USED
C            IN CONNECTION WITH L-CALCULATION PROGRAM SHELLG.
C       
C          COMMON /MODEL/ AND /IGRF1/
C               UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT>
C               ERA     EARTH RADIUS FOR NORMALIZATION OF CARTESIAN 
C                       COORDINATES (6371.2 KM)
C               AQUAD, BQUAD   SQUARE OF MAJOR AND MINOR HALF AXIS OF 
C                       EARTH ELLIPSOID AS RECOMMENDED BY INTERNAT. 
C                       ASTRONOMICAL UNION (6378.160, 6356.775 KM).
C               NMAX    MAXIMUM ORDER OF SPHERICAL HARMONICS
C               TIME    YEAR (DECIMAL: 1973.5) FOR WHICH MAGNETIC 
C                       FIELD IS TO BE CALCULATED
C               G(M)    NORMALIZED FIELD COEFFICIENTS (SEE FELDCOF)
C                       M=NMAX*(NMAX+2)
c-----------------------------------------------------------------------        
C  OUTPUT: BABS   MAGNETIC FIELD STRENGTH IN GAUSS
C          BNORTH, BEAST, BDOWN   COMPONENTS OF THE FIELD WITH RESPECT
C                 TO THE LOCAL GEODETIC COORDINATE SYSTEM, WITH AXIS
C                 POINTING IN THE TANGENTIAL PLANE TO THE NORTH, EAST
C                 AND DOWNWARD.   
C-----------------------------------------------------------------------
C-- IS RECORDS ENTRY POINT
C
C*****ENTRY POINT  FELDG  TO BE USED WITH GEODETIC CO-ORDINATES*/
	int is = 1;
	double rlat = glat * constants::umr,
		ct = sin(rlat),
		st = cos(rlat),
		d = sqrt(constants::aquad - (constants::aquad - constants::bquad)*ct*ct),
		rlon = glon * constants::umr,
		cp = cos(rlon),
		sp = sin(rlon),
		zzz = (alt + constants::bquad / d) * ct / constants::era,
		rho = (alt + constants::aquad / d) * st / constants::era,
		xxx = rho * cp,
		yyy = rho * sp,
		rq = 1. / (xxx * xxx + yyy * yyy + zzz * zzz);

	igrf2.xi[0] = xxx * rq;
	igrf2.xi[1] = yyy * rq;
	igrf2.xi[2] = zzz * rq;

	int ihmax = model.nmax * model.nmax + 1,
		last = ihmax + 2 * model.nmax,
		imax = 2 * model.nmax - 1;

	for (size_t i = ihmax - 1; i < last; i++)
		igrf2.h[i] = model.g[i];
	int i, ih, il;
	double f, x, y, z;
	for (size_t k = 1; k <= 3; k+=2) {
		i = imax;
		ih = ihmax;
		do {
			il = ih - i;
			f = 2. / (i - k + 2);
			x = igrf2.xi[0] * f;
			y = igrf2.xi[1] * f;
			z = igrf2.xi[2] * 2 * f;
			i -= 2;

			for (size_t m = 3; m <= i; m+=2) {
				igrf2.h[il + m] = model.g[il + m] + 
					z * igrf2.h[ih + m] +
					x * (igrf2.h[ih + m + 2] - igrf2.h[ih + m - 2]) -
					y * (igrf2.h[ih + m + 1] + igrf2.h[ih + m - 3]);
				igrf2.h[il + m - 1] = model.g[il + m - 1] +
					z * igrf2.h[ih + m - 1] +
					x * (igrf2.h[ih + m + 1] - igrf2.h[ih + m - 3]) +
					y * (igrf2.h[ih + m + 2] + igrf2.h[ih + m - 2]);
			}

			igrf2.h[il + 1] = model.g[il + 1] +
				z * igrf2.h[ih + 1] +
				x * igrf2.h[ih + 3] -
				y * (igrf2.h[ih + 2] + igrf2.h[ih - 1]);
			igrf2.h[il] = model.g[il] +
				z * igrf2.h[ih] +
				x * (igrf2.h[ih + 2] - igrf2.h[ih - 1]) +
				y * igrf2.h[ih + 3];
			igrf2.h[il - 1] = model.g[il - 1] +
				z * igrf2.h[ih - 1] +
				x * igrf2.h[ih] * 2 +
				y * igrf2.h[ih + 1] * 2;
			ih = il;
		} while (i >= k);
	}
	double s = 5 * igrf2.h[0] + 2 * (igrf2.h[1] * igrf2.xi[2] +
		igrf2.h[2] * igrf2.xi[0] + igrf2.h[3] * igrf2.xi[1]),
		t = 2 * rq * sqrt(rq),
		bxxx = t * (igrf2.h[2] - s * xxx),
		byyy = t * (igrf2.h[3] - s * yyy),
		bzzz = t * (igrf2.h[1] - s * zzz);

	babs = sqrt(bxxx * bxxx + byyy * byyy + bzzz * bzzz);
	beast = byyy * cp - bxxx * sp;
	double brho = byyy * sp + bxxx * cp;
	bnorth = bzzz * st - brho * ct;
	bdown = -bzzz * ct - brho * st;
};

void fieldg(double dlat, double dlong, double alt,
	double &x, double &y, double &z, double &f, double &dip,
	double &dec, double &smodip) {
/*C THIS IS A SPECIAL VERSION OF THE POGO 68/10 MAGNETIC FIELD                    
C LEGENDRE MODEL. TRANSFORMATION COEFF. G(144) VALID FOR 1973.                  
C INPUT: DLAT, DLONG=GEOGRAPHIC COORDINATES/DEG.(-90/90,0/360),                 
C        ALT=ALTITUDE/KM.                          
C OUTPUT: F TOTAL FIELD (GAUSS), Z DOWNWARD VERTICAL COMPONENT                  
C        X,Y COMPONENTS IN THE EQUATORIAL PLANE (X TO ZERO LONGITUDE).          
C        DIP INCLINATION ANGLE(DEGREE). SMODIP RAWER'S MODFIED DIP.             
C SHEIK,1977.*/
	double fel1[] = { 0.0, 0.1506723, 0.0101742, -0.0286519, 0.0092606,
		-0.0130846, 0.0089594, -0.0136808,-0.0001508, -0.0093977,
		0.0130650, 0.0020520, -0.0121956, -0.0023451, -0.0208555,
		0.0068416,-0.0142659, -0.0093322, -0.0021364, -0.0078910,
		0.0045586,  0.0128904, -0.0002951, -0.0237245,0.0289493,
		0.0074605, -0.0105741, -0.0005116, -0.0105732, -0.0058542,
		0.0033268, 0.0078164,0.0211234, 0.0099309, 0.0362792,
		-0.0201070,-0.0046350,-0.0058722,0.0011147,-0.0013949,
		-0.0108838,  0.0322263, -0.0147390,  0.0031247, 0.0111986,
		-0.0109394,0.0058112,  0.2739046, -0.0155682, -0.0253272,
		0.0163782, 0.0205730,  0.0022081, 0.0112749,-0.0098427,
		0.0072705, 0.0195189, -0.0081132, -0.0071889, -0.0579970,
		-0.0856642, 0.1884260,-0.7391512, 0.1210288, -0.0241888,
		-0.0052464, -0.0096312, -0.0044834, 0.0201764,  0.0258343,
		0.0083033,  0.0077187 },

		fel2[] = { 0.0586055,0.0102236,-0.0396107, -0.0167860,
		-0.2019911, -0.5810815,0.0379916,  3.7508268,
		1.8133030, -0.0564250, -0.0557352, 0.1335347, -0.0142641,
		-0.1024618,0.0970994, -0.0751830,-0.1274948, 0.0402073,
		0.0386290, 0.1883088,  0.1838960, -0.7848989,0.7591817,
		-0.9302389,-0.8560960, 0.6633250, -4.6363869, -13.2599277,
		0.1002136,  0.0855714,-0.0991981, -0.0765378,-0.0455264,
		0.1169326, -0.2604067, 0.1800076, -0.2223685, -0.6347679,
		0.5334222, -0.3459502,-0.1573697,  0.8589464, 1.7815990,
		-6.3347645, -3.1513653, -9.9927750,13.3327637, -35.4897308,
		37.3466339, -0.5257398,  0.0571474, -0.5421217,  0.2404770,
		-0.1747774,-0.3433644, 0.4829708,0.3935944, 0.4885033,
		0.8488121, -0.7640999, -1.8884945, 3.2930784,-7.3497229,
		0.1672821,-0.2306652, 10.5782146, 12.6031065, 8.6579742,
		215.5209961, -27.1419220,22.3405762,1108.6394043 };

	int k = 0;
	double *h = new double[144],
		*g = new double[144];
	for (size_t i = 1; i <= 72; i++)	{
		k++;
		g[k - 1] = fel1[i - 1];
		g[71 + k] = fel2[i - 1];
	}
	double rlat = dlat * constants::umr,
		ct = sin(rlat),
		st = cos(rlat),
		d = sqrt(40680925.0 - 272336.0 * ct * ct),
		rlong = dlong * constants::umr,
		cp = cos(rlong),
		sp = sin(rlong),
		zzz = (alt + 40408589.0 / d) * ct / 6371.2,
		rho = (alt + 40680925.0 / d) * st / 6371.2,
		xxx = rho * cp,
		yyy = rho * sp,
		rq = 1.0 / (xxx * xxx + yyy * yyy + zzz * zzz),
		xi[] = { xxx * rq, yyy * rq, zzz * rq };

	int nmax = 11,
		ihmax = nmax * nmax + 1,
		last = ihmax + nmax + nmax,
		imax = nmax + nmax - 1;
	for (size_t i = ihmax; i <= last; i++)
		h[i - 1] = g[i - 1];
	int i, ih, il;
	double f1, x1, y1, z1;
	for (size_t k = 1; k <= 3; k += 2) {
		i = imax;
		ih = ihmax;
		do {
			il = ih - i;
			f1 = 2. / (i - k + 2);
			x1 = xi[0] * f1;
			y1 = xi[1] * f1;
			z1 = igrf2.xi[2] * 2 * f1;
			i -= 2;
			if (i >= 1) {
				if (i == 1) {
					for (size_t m = 3; m <= i; m += 2) {
						h[il + m] = g[il + m] +
							z1 * h[ih + m] +
							x1 * (h[ih + m + 2] - h[ih + m - 2]) -
							y1 * (h[ih + m + 1] + h[ih + m - 3]);
						h[il + m - 1] = g[il + m - 1] +
							z1 * h[ih + m - 1] +
							x1 * (h[ih + m + 1] - h[ih + m - 3]) +
							y1 * (igrf2.h[ih + m + 2] + h[ih + m - 2]);
					}
				}
				h[il + 1] = g[il + 1] +
					z1 * h[ih + 1] +
					x1 * h[ih + 3] -
					y1 * (h[ih + 2] + h[ih - 1]);
				h[il] = g[il] +
					z1 * h[ih] +
					x1 * (h[ih + 2] - h[ih - 1]) +
					y1 * igrf2.h[ih + 3];
			}
			h[il - 1] = g[il - 1] +
				z1 * h[ih - 1] +
				x1 * h[ih] * 2 +
				y1 * h[ih + 1] * 2;
			ih = il;
		} while (i >= k);
	}
	double s = 5 * h[0] + 2 * (h[1] * xi[2] +
		h[2] * xi[0] + h[3] * xi[1]),
		xt = 2 * rq * sqrt(rq);
	x = xt * (h[2] - s * xxx);
	y = xt * (h[3] - s * yyy);
	z = xt * (h[1] - s * zzz);
	f = sqrt(x*x + y * y + z * z);
	double brh0 = y * sp + x * cp;
	y = y * cp - x * sp;
	x = z * st - brh0 * ct;
	z = -z * ct - brh0 * st;
	double zdivf = z / f;
	if (abs(zdivf) > 1.) zdivf /= abs(zdivf);
	dip = asin(zdivf);
	double ydivs = y / sqrt(x*x + y * y);
	if (abs(ydivs) > 1.) ydivs /= abs(ydivs);
	dec = asin(ydivs);
	double dipdiv = dip / sqrt(dip*dip + st);
	if (abs(dipdiv) > 1.) dipdiv /= abs(dipdiv);
	smodip = asin(dipdiv);
	dip /= constants::umr;
	dec /= constants::umr;
	smodip /= constants::umr;
};

void igrf_sub(double xlat, double xlong, double year, double height,
	double &xl, int &icode, double &dipl, double &babs) {
/*c-----------------------------------------------------------------------        
c INPUT:
c    xlat      geodatic latitude in degrees
c    xlong     geodatic longitude in degrees
c    year      decimal year (year+(month-0.5)/12.0-0.5 or 
c                  year+day-of-year/365 or ../366 if leap year) 
c    height    height in km
c OUTPUT:
c    xl        L value
c    icode      =1  L is correct; =2  L is not correct;
c               =3  an approximation is used
c    dipl      dip latitude in degrees
c    babs      magnetic field strength in Gauss
c-----------------------------------------------------------------------    */
	double lati = xlat,
		longi = xlong,
		bnorth, beast, bdown, bab1;
	feldg(lati, longi, height, bnorth, beast, bdown, babs);
	shellg(lati, longi, height, xl, icode, bab1);
	dipl = atan(bdown / 2.0 / sqrt(bnorth*bnorth + beast * beast)) / constants::umr;
};

void shellg(double glat, double glon, double alt,
	double &fl, int &icode, double &b0) {
/*c-----------------------------------------------------------------------        
C CALCULATES L-VALUE FOR SPECIFIED GEODAETIC COORDINATES, ALTITUDE
C AND GEMAGNETIC FIELD MODEL.
C REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTER, INTERNAL NOTE 
C      NO. 67, 1970.
C      G. KLUGE, COMPUTER PHYSICS COMMUNICATIONS 3, 31-35, 1972
c-----------------------------------------------------------------------        
C CHANGES (D. BILITZA, NOV 87):
C   - USING CORRECT DIPOL MOMENT I.E.,DIFFERENT COMMON/MODEL/
C   - USING IGRF EARTH MAGNETIC FIELD MODELS FROM 1945 TO 1990
C 09/07/22 NMAX=13 for DGRF00 and IGRF05; H/G-arrays(195)
c-----------------------------------------------------------------------        
C  INPUT:  ENTRY POINT SHELLG
C             GLAT  GEODETIC LATITUDE IN DEGREES (NORTH)
C             GLON  GEODETIC LONGITUDE IN DEGREES (EAST)
C             ALT   ALTITUDE IN KM ABOVE SEA LEVEL
C
C          ENTRY POINT SHELLC
C             V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM)
C                     X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE
C                     Y-AXIS POINTING TO EQUATOR AT 90 LONG.
C                     Z-AXIS POINTING TO NORTH POLE
C
C          DIMO     DIPOL MOMENT IN GAUSS (NORMALIZED TO EARTH RADIUS) 
C
C          COMMON 
C             X(3)    NOT USED
C             H(144)  FIELD MODEL COEFFICIENTS ADJUSTED FOR SHELLG
c-----------------------------------------------------------------------        
C  OUTPUT: FL           L-VALUE
C          ICODE        =1 NORMAL COMPLETION
C                       =2 UNPHYSICAL CONJUGATE POINT (FL MEANINGLESS)
C                       =3 SHELL PARAMETER GREATER THAN LIMIT UP TO
C                          WHICH ACCURATE CALCULATION IS REQUIRED;
C                          APPROXIMATION IS USED.
C          B0           MAGNETIC FIELD STRENGTH IN GAUSS
c-----------------------------------------------------------------------       
C-- RMIN, RMAX ARE BOUNDARIES FOR IDENTIFICATION OF ICODE=2 AND 3
C-- STEP IS STEP SIZE FOR FIELD LINE TRACING
C-- STEQ IS STEP SIZE FOR INTEGRATION*/
	double rmin = 0.05,
		rmax = 1.01,
		step = 0.2,
		steq = 0.03,
		bequ = 1e10;
	// C*****ENTRY POINT  SHELLG  TO BE USED WITH GEODETIC CO-ORDINATES
	double rlat = glat * constants::umr,
		ct = sin(rlat),
		st = cos(rlat),
		d = sqrt(constants::aquad - (constants::aquad - constants::bquad)*ct*ct);
	igrf2.x[0] = (alt + constants::aquad / d)*st / constants::era;
	igrf2.x[2] = (alt + constants::bquad / d)*ct / constants::era;
	double rlon = glon * constants::umr;
	igrf2.x[1] = igrf2.x[0] * sin(rlon);
	igrf2.x[0] *= cos(rlon);
	double **p = new double*[8];
	for (size_t i = 0; i < 8; i++)
		p[i] = new double[100];
	double u[3][3] = {
		{+0.3511737,-0.9148385,-0.1993679},
		{+0.9335804,+0.3583680,+0.0000000},
		{+0.0714471,-0.1861260,+0.9799247}
	};
	double rq = 1. / (igrf2.x[0] * igrf2.x[0] + igrf2.x[1] * igrf2.x[1] + igrf2.x[2] * igrf2.x[2]),
		r3h = sqrt(rq * sqrt(rq)), bq1, r1, bq2, bq3, r2, r3, zz;
	p[0][1] = (igrf2.x[0] * u[0][0] + igrf2.x[1] * u[1][0] + igrf2.x[2] * u[2][0]) * r3h;
	p[1][1] = (igrf2.x[0] * u[0][1] + igrf2.x[1] * u[1][1]) * r3h;
	p[2][1] = (igrf2.x[0] * u[0][2] + igrf2.x[1] * u[1][2] + igrf2.x[2] * u[2][2]) * rq;

	// C*****FIRST THREE POINTS OF FIELD LINE
	step = -copysign(step, p[2][1]);
	stoer(p, 1, bq2, r2);
	b0 = sqrt(bq2);
	p[0][2] = p[0][1] + 0.5 * step * p[3][2];
	p[1][2] = p[1][1] + 0.5 * step * p[4][1];
	p[2][2] = p[2][1] + 0.5 * step;
	stoer(p, 2, bq3, r3);
	p[0][0] = p[0][1] - step * (2. * p[3][1] - p[3][2]);
	p[1][0] = p[1][1] - step * (2. * p[4][1] - p[4][2]);
	p[2][0] = p[2][1] - step;
	stoer(p, 0, bq1, r1);
	p[0][2] = p[0][1] + step * (20. * p[3][2] - 3. * p[3][1] + p[3][0]) / 18.;
	p[1][2] = p[1][1] + step * (20. * p[4][2] - 3. * p[4][1] + p[4][0]) / 18.;
	p[2][2] = p[2][1] + step;
	stoer(p, 2, bq3, r3);
	// C*****INVERT SENSE IF REQUIRED
	if (bq3 > bq1) {
		step = -step;
		r3 = r1;
		bq3 = bq1;
		for (size_t i = 1; i <= 7; i++)	{
			zz = p[i - 1][0];
			p[i - 1][0] = p[i - 1][2];
			p[i - 1][2] = zz;
		}
	}
	// C*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
	int iequ = 0;
	if (bq1 < bequ) {
		bequ = bq1;
		iequ = 1;
	}
	if (bq2 < bequ) {
		bequ = bq2;
		iequ = 2;
	}
	if (bq3 < bequ) {
		bequ = bq3;
		iequ = 1;
	}
	// C*****INITIALIZATION OF INTEGRATION LOOPS 
	double step12 = step / 12,
		step2 = step + step,
		fi = 0,
		oradik = 0,
		oterm = 0,
		stp = r2 * steq,
		z = p[2][1] + stp;
	steq = copysign(steq, step);
	stp /= 0.75;
	p[7][0] = step2 * (p[0][0] * p[3][0] + p[1][0] * p[4][0]);
	p[7][1] = step2 * (p[0][1] * p[3][1] + p[1][1] * p[4][1]);
	// C*****MAIN LOOP (FIELD LINE TRACING)
	double t, hli, zq, r, ff, radik, term,
		c0, c1, c2, c3, d0, d1, d2, e0, e1, e2;
	bool goto30 = false,
		goto10 = false;
	int num;
	for (size_t n = 3; n <= 3333; n++) {
		// C*****CORRECTOR (FIELD LINE TRACING)
		p[0][n - 1] = p[0][n - 2] + step12 * (5. * p[3][n - 1] + 8. * p[3][n - 2] - p[3][n - 3]);
		p[1][n - 1] = p[1][n - 2] + step12 * (5. * p[4][n - 1] + 8. * p[4][n - 2] - p[4][n - 3]);
		// C*****PREPARE EXPANSION COEFFICIENTS FOR INTERPOLATION           
		// C*****OF SLOWLY VARYING QUANTITIES
		p[7][n - 1] = step2 * (p[0][n - 1] * p[3][n - 1] + p[1][n - 1] * p[4][n - 1]);
		c0 = p[0][n - 2] * p[0][n - 2] + p[1][n - 2] * p[1][n - 2];
		c1 = p[7][n - 2];
		c2 = (p[7][n - 1] - p[7][n - 3]) * 0.25;
		c3 = (p[7][n - 1] + p[7][n - 3] - c1 - c1) / 6.0;
		d0 = p[5][n - 2];
		d1 = (p[5][n - 1] - p[5][n - 3]) * 0.5;
		d2 = (p[5][n - 1] + p[5][n - 3] - d0 - d0) * 0.5;
		e0 = p[6][n - 2];
		e1 = (p[6][n - 1] - p[6][n - 3]) * 0.5;
		e2 = (p[6][n - 1] + p[6][n - 3] - e0 - e0) * 0.5;
		// C*****INNER LOOP (FOR QUADRATURE)
		while (true) {
			t = (z - p[2][n - 2]) / step;
			if (t > 1.) break;
			hli = 0.5 * (((c3 * t + c2) * t + c1) * t + c0);
			zq = z * z;
			r = hli + sqrt(hli * hli + zq);
			if (r <= rmin) {
				goto30 = true;
				num = n;
				break;
			}
			rq = r * r;
			ff = sqrt(1. + 3. * zq / rq);
			radik = b0 - ((d2 * t + d1) * t + d0) * r * rq * ff;
			if (r - rmax == 0) {
				icode = 2;
				radik = radik - 12. * (r - rmax) * (r - rmax);
			}
			if (radik + radik <= oradik) {
				goto10 = true;
				break;
			}
			term = sqrt(radik) * ff * ((e2*t + e1) * t + e0) / (rq + zq);
			fi = fi + stp * (oterm + term);
			oradik = radik;
			oterm = term;
			stp = r * steq;
			z = z + stp;
		}
		if (goto10 || goto30) break;
		p[0][n] = p[0][n - 1] + step12 * (23. * p[3][n - 1] - 16. * p[3][n - 2] + 5. * p[3][n - 3]);
		p[1][n] = p[1][n - 1] + step12 * (23. * p[4][n - 1] - 16. * p[4][n - 2] + 5. * p[4][n - 3]);
		p[2][n] = p[2][n - 1] + step;
		stoer(p, n, bq3, r3);
		// C*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
		if (bq3 < bequ) {
			bequ = bq3;
			iequ = n + 1;
		}
	}
	if (!goto30) {
		if (iequ < 2) iequ = 2;
		fidb0.sp[0] = p[0][iequ - 2];
		fidb0.sp[1] = p[1][iequ - 2];
		fidb0.sp[2] = p[2][iequ - 2];
		if (oradik >= 1e-15)
			fi += stp / 0.75 * oterm * oradik / (oradik - radik);
		/*C-- The minimal allowable value of FI was changed from 1E-15 to 1E-12,
		C-- because 1E-38 is the minimal allowable arg. for ALOG in our envir.
		C-- D. Bilitza, Nov 87.*/
		fi = 0.5 * abs(fi) / sqrt(b0) + 1e-12;
		/*C*****COMPUTE L FROM B AND I.  SAME AS CARMEL IN INVAR.  
		C
		C-- Correct dipole moment is used here. D. Bilitza, Nov 87.*/
		double dimob0 = constants::dimo / b0,
			arg1 = log(fi),
			arg2 = log(dimob0),
			xx = 3 * arg1 - arg2, gg;

		if (xx > 23)
			gg = xx - 3.0460681e0;
		else if (xx > 11.7)
			gg = (((((2.8212095e-8 * xx - 3.8049276e-6) * xx + 2.170224e-4) * xx -
				6.7310339e-3) * xx + 1.2038224e-1) * xx - 1.8461796e-1) * xx + 2.0007187e0;
		else if (xx > 3)
			gg = ((((((((6.3271665e-10 * xx - 3.958306e-8) * xx + 9.9766148e-07) * xx -
				1.2531932e-5) * xx + 7.9451313e-5) * xx - 3.2077032e-4) * xx +
				2.1680398e-3) * xx + 1.2817956e-2) * xx + 4.3510529e-1) * xx + 6.222355e-1;
		else if (xx > -3)
			gg = ((((((((2.6047023e-10 * xx + 2.3028767e-9) * xx - 2.1997983e-8) * xx -
				5.3977642e-7) * xx - 3.3408822e-6) * xx + 3.8379917e-5) * xx +
				1.1784234e-3) * xx + 1.4492441e-2) * xx + 4.3352788e-1) * xx + 6.228644e-1;
		else if (xx > -22)
			gg = ((((((((-8.1537735e-14 * xx + 8.3232531e-13) * xx + 1.0066362e-9) * xx +
				8.1048663e-8) * xx + 3.2916354e-6) * xx + 8.2711096e-5) * xx +
				1.3714667e-3) * xx + 1.5017245e-2) * xx + 4.3432642e-1) * xx + 6.2337691e-1;
		else
			gg = 3.33338e-1 * xx + 3.0062102e-1;

		fl = exp(log((1. + exp(gg)) * dimob0) / 3.0);
	}
	else {
		icode = 3;
		t = -p[2][num - 2] / step;
		fl = 1. / (abs(((c3 * t + c2) * t + c1) * t + c0) + 1e-15);
	}
};

void stoer(double **p, int idx, double &bq, double &r) {
/*C*******************************************************************
C* SUBROUTINE USED FOR FIELD LINE TRACING IN SHELLG                *
C* CALLS ENTRY POINT FELDI IN GEOMAGNETIC FIELD SUBROUTINE FELDG   *
C
C 09/07/22 NMAX=13 for DGRF00 and IGRF05; H/G-arrays(195)
C********************************************************************/
	// C*****XM,YM,ZM  ARE GEOMAGNETIC CARTESIAN INVERSE CO-ORDINATES
	double zm = p[2][idx],
		fli = p[0][idx] * p[0][idx] + p[1][idx] * p[1][idx] + 1e-15;
	r = 0.5 * (fli + sqrt(fli * fli + 4 * zm * zm));
	double rq = r * r,
		wr = sqrt(r),
		xm = p[0][idx] * wr,
		ym = p[1][idx] * wr;
	// C*****TRANSFORM TO GEOGRAPHIC CO-ORDINATE SYSTEM
	double u[3][3] = {
		{+0.3511737,-0.9148385,-0.1993679},
		{+0.9335804,+0.3583680,+0.0000000},
		{+0.0714471,-0.1861260,+0.9799247}
	};
	igrf2.xi[0] = xm * u[0][0] + ym * u[0][1] + zm * u[0][2];
	igrf2.xi[1] = xm * u[2][0] + ym * u[1][1] + zm * u[1][2];
	igrf2.xi[2] = xm * u[2][0] + zm * u[2][2];
	// C*****COMPUTE DERIVATIVES
	feldi();
	double q = igrf2.h[0] / rq,
		dx = igrf2.h[2] + igrf2.h[2] + q * igrf2.xi[0],
		dy = igrf2.h[3] + igrf2.h[3] + q * igrf2.xi[1],
		dz = igrf2.h[1] + igrf2.h[1] + q * igrf2.xi[2];
	// C*****TRANSFORM BACK TO GEOMAGNETIC CO - ORDINATE SYSTEM
	double dxm = u[0][0] * dx + u[1][0] * dy + u[2][0] * dz,
		dym = u[0][1] * dx + u[1][1] * dy,
		dzm = u[0][2] * dx + u[1][2] * dy + u[2][2] * dz,
		dr = (xm * dxm + ym * dym + zm * dzm) / r;
	// C*****FORM SLOWLY VARYING EXPRESSIONS
	p[3][idx] = (wr * dxm - 0.5 * p[0][idx] * dr) / (r * dzm);
	p[4][idx] = (wr * dym - 0.5 * p[1][idx] * dr) / (r * dzm);
	double dsq = rq * (dxm * dxm + dym * dym + dzm * dzm);
	bq = dsq * rq * rq;
	p[5][idx] = sqrt(dsq / (rq + 3. * zm * zm));
	p[6][idx] = p[5][idx] * (rq + zm * zm) / (rq * dzm);
};

void feldi() {
	int ihmax = model.nmax * model.nmax + 1,
		last = ihmax + 2 * model.nmax,
		imax = 2 * model.nmax - 1;

	for (size_t i = ihmax - 1; i < last; i++)
		igrf2.h[i] = model.g[i];
	int i, ih, il;
	double f, x, y, z;
	for (size_t k = 1; k <= 3; k += 2) {
		i = imax;
		ih = ihmax;
		do {
			il = ih - i;
			f = 2. / (i - k + 2);
			x = igrf2.xi[0] * f;
			y = igrf2.xi[1] * f;
			z = igrf2.xi[2] * 2 * f;
			i -= 2;

			for (size_t m = 3; m <= i; m += 2) {
				igrf2.h[il + m] = model.g[il + m] +
					z * igrf2.h[ih + m] +
					x * (igrf2.h[ih + m + 2] - igrf2.h[ih + m - 2]) -
					y * (igrf2.h[ih + m + 1] + igrf2.h[ih + m - 3]);
				igrf2.h[il + m - 1] = model.g[il + m - 1] +
					z * igrf2.h[ih + m - 1] +
					x * (igrf2.h[ih + m + 1] - igrf2.h[ih + m - 3]) +
					y * (igrf2.h[ih + m + 2] + igrf2.h[ih + m - 2]);
			}

			igrf2.h[il + 1] = model.g[il + 1] +
				z * igrf2.h[ih + 1] +
				x * igrf2.h[ih + 3] -
				y * (igrf2.h[ih + 2] + igrf2.h[ih - 1]);
			igrf2.h[il] = model.g[il] +
				z * igrf2.h[ih] +
				x * (igrf2.h[ih + 2] - igrf2.h[ih - 1]) +
				y * igrf2.h[ih + 3];
			igrf2.h[il - 1] = model.g[il - 1] +
				z * igrf2.h[ih - 1] +
				x * igrf2.h[ih] * 2 +
				y * igrf2.h[ih + 1] * 2;
			ih = il;
		} while (i >= k);
	}
}

void clcmlt(int iyyyy, int ddd, double uthr, double glat, double glon, double &mlt) {
/*C--------------------------------------------------------------------
C      calculates magnetic local time
C      Inputs:
C             IYYYY..Year as YYYY, e.g. 1998
C             DDD..day of year (1.1. = 0)
C             UTHR..universal time in decimal hours
C             GLAT,GLON..latitude north and longitude east in degrees
C      Output:
C             MLT..magnetic local time in decimal hours
C      Required subroutines: DPMTRX
C--------------------------------------------------------------------*/
	double xg = cos(glat * constants::umr) * cos(glon * constants::umr),
		yg = cos(glat * constants::umr) * sin(glon * constants::umr),
		zg = sin(glat * constants::umr);

	double *xxm = new double[3],
		*yym = new double[3],
		*zzm = new double[3];
	dpmtrx(iyyyy, ddd, xxm, yym, zzm);
	double xm = xxm[0] * xg + xxm[1] * yg + xxm[2] * zg,
		ym = yym[0] * xg + yym[1] * yg + yym[2] * zg,
		zm = zzm[0] * xg + zzm[1] * yg + zzm[2] * zg,
		gst, slong, srasn, sdec;

	int ihour = int(uthr),
		min = int((uthr - ihour) * 60),
		isec = int((uthr - ihour - min / 60.0) * 3600);

		sun(iyyyy, ddd + 1, ihour, min, isec, gst, slong, srasn, sdec);
		double be = gst,
			cal = cos(srasn),
			*sa = new double[3],
			*sg = new double[3],
			*sm = new double[3];
		sa[2] = sin(sdec);
		sa[0] = cos(sdec);
		sa[1] = sa[0] * sin(srasn);
		sa[0] = sa[0] * cal;
		double s = sin(be),
			c = cos(be);

		sg[0] = c * sa[0] + s * sa[1];
		sg[1] = c * sa[1] - s * sa[0];
		sg[2] = sa[2];

		// transform
		sm[0] = xxm[0] * sg[0] + xxm[1] * sg[1] + xxm[2] * sg[2];
		sm[1] = yym[0] * sg[0] + yym[1] * sg[1] + yym[2] * sg[2];
		sm[2] = zzm[0] * sg[0] + zzm[1] * sg[1] + zzm[2] * sg[2];

		double lam = atan2(ym, xm),
			lams = atan2(sm[1], sm[0]),
			dellam = lam - lams;
		if (dellam < 0.) dellam += 2 * constants::pi;
		mlt = remainder(dellam / constants::pi * 12. + 12., 24.);
}

void dpmtrx(int iyyyy, int ddd, double *xm, double *ym, double *zm) {
/*C--------------------------------------------------------------------------
C      calculates othonormal matrix (columns XM,YM,ZM) for transformation 
C      from geographic to magnetic coordinates
C      Inputs:
C             IYYYY..year
C               DDD..day of year (1.1 = 0)
C      Outputs:
C               XM,YM,ZM..colums of the matrix
C      Notes:
C      MX(N),MY(N),MZ(N)..coordinates of the B vector in geographic system 
C                for years stored in YR(N)
C      N..number of elements of arrays MX,MY,MZ and YR
C--------------------------------------------------------------------------*/
	double mxi = -dipol.ghi2,
		myi = -dipol.ghi3,
		mzi = -dipol.ghi1,
		m = sqrt(mxi * mxi + myi * myi + mzi * mzi),
		myz = sqrt(myi * myi + mzi * mzi);

	zm[0] = mxi / m;
	zm[1] = myi / m;
	zm[2] = mzi / m;
	double zm12 = sqrt(zm[0] * zm[0] + zm[1] * zm[1]);
	ym[0] = -zm[1] / zm12;
	ym[1] = zm[0] / zm12;
	ym[2] = 0.;
	xm[0] = ym[1] * zm[2] - ym[2] * zm[1];
	xm[1] = ym[2] * zm[0] - ym[0] * zm[2];
	xm[2] = ym[0] * zm[1] - ym[1] * zm[0];
}

void sun(int iyear, int iday, int ihour, int min, int isec,
	double &gst, double &slong, double &srasn, double &sdec) {
/*C-----------------------------------------------------------------------------
C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
C
C-------  INPUT PARAMETERS:
C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
C
C-------  OUTPUT PARAMETERS:
C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
C  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
C  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
C
C  LAST MODIFICATION:  MARCH 31, 2003 (ONLY SOME NOTATION CHANGES)
C
C     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
C-----------------------------------------------------------------------------*/
	if (iyear >= 1901 && iyear <= 2099) return;

	double fday = 1. * (ihour * 3600 + min * 60 + isec) / 86400,
		dj = 365 * (iyear - 1900) + 1. * (iyear - 1901) / 4 + iday - 0.5 + fday,
		t = dj / 36525.,
		vl = remainder(279.696678 + 0.9856473354 * dj, 360.),
		g = remainder(358.475845 + 0.985600267 * dj, 360.) * constants::umr;
	gst = remainder(279.690983 + .9856473354 * dj + 360. * fday + 180., 360.) * constants::umr;
	slong = (vl + (1.91946 - 0.004789*t) * sin(g) + 0.020094 * sin(2. * g)) * constants::umr;
	if (slong > 6.2831853) slong -= 2 * constants::pi;
	else if (slong < 0) slong += 2 * constants::pi;
	double obliq = (23.45229 - 0.0130125*t) * constants::umr,
		sob = sin(obliq),
		slp = slong - 9.924e-5;

	/*C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  DUE TO
	C   THE ORBITAL MOTION OF THE EARTH*/
	double sin1 = sob * sin(slp),
		cos1 = sqrt(1. - sin1 * sin1),
		sc = sin1 / cos1;
	sdec = atan(sc);
	srasn = constants::pi - atan2(cos(obliq) / sob * sc, -cos(slp) / cos1);
}

void tcon(int yr, int mm, int day, int idn, double *rz,
	double *ig, double &rsn, int &nmonth) {
/*c----------------------------------------------------------------
c input:        yr,mm,day       year(yyyy),month(mm),day(dd)
c               idn             day of year(ddd)
c output:       rz(3)           12-month-smoothed solar sunspot number
c               ig(3)           12-month-smoothed IG index
c               rsn             interpolation parameter
c               nmonth          previous or following month depending
c                               on day
c
c Uses read_ig_rz and common/igrz/ to get indices 
c 
c rz(1) & ig(1) contain the indices for the month mm and rz(2) & ig(2)
c for the previous month (if day less than 15) or for the following
c month (if day greater than 15). These indices are for the mid of the 
c month. The indices for the given day are obtained by linear 
c interpolation and are stored in rz(3) and ig(3).
c----------------------------------------------------------------*/
	double iytmp = yr * 100 + mm;
	if (iytmp < igrz.iymst || iytmp > igrz.iymend) {
		nmonth = -1;
		return;
	}
	
	int iyst = igrz.iymst / 100,
		imst = igrz.iymst - iyst * 100,
		// num = 12 - imst + 1 + (yr - iyst - 1) * 12 + mm + 1
		num = 2 - imst + (yr - iyst) * 12 + mm;

	rz[0] = igrz.arz[num - 1];
	ig[0] = igrz.aig[num - 1];
	int midm = 15;
	if (mm == 2) midm = 14;
	int idd1, nrdaym;
	moda(0, yr, mm, midm, idd1, nrdaym);
	int imm2, idd2, iyy2;
	if (day < midm) {
		// day is at or after mid of month
		imm2 = mm + 1;
		if (imm2 > 12) {
			imm2 = 1;
			iyy2 = yr + 1;
			idd2 = 380;
			if (yr / 4 * 4 == yr) idd2 = 381;
		}
		else {
			iyy2 = yr;
			midm = 15;
			if (imm2 == 2) midm = 14;
			moda(0, iyy2, imm2, midm, idd2, nrdaym);
		}
		rz[1] = igrz.arz[num];
		ig[1] = igrz.aig[num];
		rsn = (idn - idd1) * 1. / (idd2 - idd1);
		rz[2] = rz[0] + (rz[1] - rz[0]) * rsn;
		ig[2] = ig[0] + (ig[1] - ig[0]) * rsn;
	}
	else {
		imm2 = mm - 1;
		if (imm2 < 1) {
			imm2 = 12;
			idd2 = -16;
			iyy2 = yr - 1;
		}
		else {
			iyy2 = yr;
			midm = 15;
			if (imm2 == 2) midm = 14;
			moda(0, iyy2, imm2, midm, idd2, nrdaym);
		}
		rz[1] = igrz.arz[num - 2];
		ig[1] = igrz.aig[num - 2];
		rsn = (idn - idd2) * 1. / (idd1 - idd2);
		rz[2] = rz[1] + (rz[0] - rz[1]) * rsn;
		ig[2] = ig[1] + (ig[0] - ig[1]) * rsn;
	}
	nmonth = imm2;
}

void apf_only(int iyyyy, int imn, int id, double &f107d, double &f107pd,
	double &f107_81, double &f107_365, int &iapda, int &isdate) {
/*c-----------------------------------------------------------------------
c Finds daily F10.7, daily Ap, and 81-day and 365-day F10.7 index: 
c
c    INPUTS:   IYYYY (yyyy)	year 
c              IMN (mm)		month 
c	           ID (dd)		day 
c    OUTPUT:   F107D		F10.7 index for the day (adjusted 
c								to 1AU)
C              F107PD  		F10.7 index for one day prior (used in MSIS)
c              F107_81		F10.7 average over 3 solar rotations
c                               (81 days, centered on the current day) 
c              F107_365 	F10.7 12-month running mean
c              IAPDA		Daily Ap
c              ISDATE		Array-index for the specified date (for
c                                use in APF subroutine.
c 
c Using COMMON/apfa/ for indices
c
c Is used for vdrift and foeedi.
c
c If date is outside the range of indices file than F107D=F107_81=-11.1  
c-----------------------------------------------------------------------*/
	int lm[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 },
		iybeg = 1958, is, nyd;
	f107d = -11.1;
	f107_81 = -11.1;
	f107_365 = -11.1;
	iapda = -11;
	if (iyyyy >= iybeg) {
		is = 0;
		for (size_t i = iybeg; i < iyyyy; i++) {
			nyd = 365;
			if (i / 4 * 4 == i) nyd = 366;
			is += nyd;
		}
		lm[1] = 28;
		if (iyyyy / 4 * 4 == iyyyy) lm[1] = 29;
		for (size_t i = 1; i < imn; i++)
			is += lm[i - 1];
		is += id;
		isdate = is;
		if (is <= apfa.n) {
			f107d = apfa.af107[is - 1][0];
			f107pd = f107d;
			if (is > 1) f107pd = apfa.af107[is - 2][0];
			f107_81 = apfa.af107[is - 1][1];
			if (f107_81 < -4) f107_81 = f107d;
			f107_365 = apfa.af107[is - 1][2];
			if (f107_365 < -4) f107_365 = f107d;
			iapda = apfa.aap[is - 1][8];
		}
	}
}

void soco(int ld, double t, double flat, double elon, double height,
	double &declin, double &zenith, double &sunrse, double &sunset) {
/*c--------------------------------------------------------------------
c       s/r to calculate the solar declination, zenith angle, and
c       sunrise & sunset times  - based on Newbern Smith's algorithm
c       [leo mcnamara, 1-sep-86, last modified 16-jun-87]
c       {dieter bilitza, 30-oct-89, modified for IRI application}
c
c in:   ld      local day of year
c       t       local hour (decimal)
c       flat    northern latitude in degrees
c       elon    east longitude in degrees
c		height	height in km
c
c out:  declin      declination of the sun in degrees
c       zenith      zenith angle of the sun in degrees
c       sunrse      local time of sunrise in hours 
c       sunset      local time of sunset in hours 
c-------------------------------------------------------------------*/
	//  amplitudes of Fourier coefficients  --  1955 epoch.................
	double  p[] = { 0.017203534, 0.034407068, 0.051610602, 0.068814136, 0.103221204 };
	// s/r is formulated in terms of WEST longitude.......................
	double wlon = 360 - elon;
	// time of equinox for 1980...........................................
	double td = ld + (t + wlon / 15.) / 24.,
		te = td + 0.9369;
	// declination of the sun..............................................
	double dcl = 23.256 * sin(p[0] * (te - 82.242)) + 0.381 * sin(p[1] * (te - 44.855))
		+ 0.167 * sin(p[2] * (te - 23.355)) - 0.013 * sin(p[3] * (te + 11.97))
		+ 0.011 * sin(p[5] * (te - 10.41)) + 0.339137;
	declin = dcl;
	double dc = dcl * constants::umr;
	// the equation of time................................................
	double tf = te - 0.5,
		eqt = -7.38 * sin(p[0] * (tf - 4.)) - 9.87 * sin(p[1] * (tf + 9.))
		+ 0.27 * sin(p[2] * (tf - 53.)) - 0.2 * cos(p[3] * (tf - 17.)),
		et = eqt * constants::umr / 4.,
		fa = flat * constants::umr,
		phi = constants::humr * (t - 12.) + et,
		a = sin(fa) * sin(dc),
		b = cos(fa) * cos(dc),
		cosx = a + b * cos(phi);
	if (abs(cosx) > 1.) cosx = copysign(1., cosx);
	zenith = acos(cosx) / constants::umr;
	/*c calculate sunrise and sunset times --  at the ground...........
	c see Explanatory Supplement to the Ephemeris (1961) pg 401......
	c sunrise at height h metres is at...............................*/
	double h = height * 1000.,
		chih = 90.83 + 0.0347 * sqrt(h);
	/*c this includes corrections for horizontal refraction and........
	c semi-diameter of the solar disk................................*/
	double ch = cos(chih * constants::umr),
		cosphi = (ch - a) / b;
	/*c if abs(secphi) > 1., sun does not rise/set.....................
	c allow for sun never setting - high latitude summer.............*/
	double secphi = 999999.;
	if (cosphi != 0) secphi = 1. / cosphi;
	
	if ((secphi > -1) && (secphi <= 0)) {
		sunset = 99.;
		sunrse = 99.;
	}
	// c allow for sun never rising - high latitude winter..............
	else if ((secphi > 0) && (secphi < 1)) {
		sunset = -99.;
		sunrse = -99.;
	}
	else {
		cosx = cosphi;
		if (abs(cosx) > 1) cosx = copysign(1., cosx);
		phi = acos(cosx);
		et = et / constants::humr;
		phi = phi / constants::humr;
		sunrse = 12. - phi - et;
		sunset = 12. + phi - et;
		if (sunrse < 0) sunrse += 24;
		if (sunset >= 24.) sunset -= 24;
	}
}

double foeedi(double cov, double xhi, double xhim, double xlati) {
/*C-------------------------------------------------------
C CALCULATES FOE/MHZ BY THE EDINBURGH-METHOD.      
C INPUT: 
C 	COV		MONTHLY MEAN 10.7CM SOLAR RADIO FLUX measured at 
C           ground level  
C   XHI		SOLAR ZENITH ANGLE IN DEGREE 
C   XHIM	SOLAR ZENITH ANGLE AT NOON IN DEGREE
C   XLATI 	ABSOLUTE VALUE OF GEOGRAPHIC LATITUDE IN DEGREE, 
C REFERENCE: 
C       KOURIS-MUGGELETON, CCIR DOC. 6/3/07, 1973
C       TROST, J. GEOPHYS. RES. 84, 2736, 1979 (was used
C               to improve the nighttime varition)
C       RAWER AND BILITZA, Adv. Space Res. 10(8), 5-14, 1990
C D.BILITZA--------------------------------- AUGUST 1986.*/
	// variation with solar activity (factor A) ...............
	double a = 1.0 + 0.0094 * (cov - 66.0);
	// variation with noon solar zenith angle (B) and with latitude (C)
	double sl = cos(xlati * constants::umr),
		sm, c, b, sp;
	if (xlati < 32) {
		sm = -1.93 + 1.92 * sl;
		c = 23.0 + 116.0 * sl;
	}
	else {
		sm = 0.11 - 0.49 * sl;
		c = 92.0 + 35.0 * sl;
	}
	if (xhim >= 90) xhim = 89.999;
	b = pow(cos(xhim * constants::umr), sm);
	// variation with solar zenith angle (D) ..........................        
	if (xlati > 12) sp = 1.2;
	else sp = 1.31;
	// adjusted solar zenith angle during nighttime (XHIC) .............
	double xhic = xhi - 3 * log(1. + exp((xhi - 89.98) / 3.)),
		d = pow(cos(xhic * constants::umr), sp);
	// determine foE**4 ................................................
	double r4foe = a * b * c * d;
	// minimum allowable foE (foe_min=sqrt[SMIN])...............................
	double smin = 0.121 + 0.0015*(cov - 60.);
	smin *= smin;
	if (r4foe < smin) r4foe = smin;
	return pow(r4foe, 0.25);
}

void storm(int *ap, double rga, double rgo, int coor, double &rgma, int ut,
	int doy, double &cf) {
	/*C----------------------------------------------------------------------
	C      Fortran code to obtain the foF2 storm-time correction factor at
	C      a given location and time, using the current and the 12 previous
	C      ap values as input.
	C
	C      ap ---> (13 elements integer array). Array with the preceeding
	C              13 value of the 3-hourly ap index. The 13th value
	C              in the array will contain the ap at the UT of interest,
	C              the 12th value will contain the 1st three hourly interval
	C              preceeding the time of interest, and so on to the 1st
	C              ap value at the earliest time.
	C     coor --> (integer). If coor = 2, rga should contain the
	C                         geomagnetic latitude.
	C                         If coor = 1, rga should contain the
	C                         geographic latitude.
	C     rga ---> (real, -90 to 90) geographic or geomagnetic latitude.
	C     rgo ---> (real, 0 to 360, positive east from Greenwich.)
	C                           geographic longitude, only used if coor=1.
	C     ut  ---> (integer, hours 00 to 23) Universal Time of interest.
	C     doy ---> (integer, 1 to 366)Day of the year.
	C     cf  ---> (real) The output; the storm-time correction factor used
	C              to scale foF2, foF2 * cf.
	C     rgma --> corrected magnetic latitude calculated from rga and rgo
	C
	C     This model and computer code was developed by E. Araujo-Pradere,
	C     T. Fuller-Rowell and M. Condrescu, SEC, NOAA, Boulder, USA
	C     Ref:
	C     T. Fuller-Rowell, E. Araujo-Pradere, and M. Condrescu, An
	C       Empirical Ionospheric Storm-Time Ionospheric Correction Model,
	C       Adv. Space Res. 8, 8, 15-24, 2000.
	C----------------------------------------------------------------------
	C     DIMENSIONS AND COEFFICIENTS VALUES*/
	double c4[20] = { 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
		0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
		0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00 },

		c3[20] = { 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, -9.44e-12,
		0.00e+00, 3.04e-12, 0.00e+00, 9.32e-12, -1.07e-11, 0.00e+00, 0.00e+00,
		0.00e+00, 1.09e-11, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, -1.01e-11 },

		c2[20] = { 1.16e-08, 0.00e+00, 0.00e+00, -1.46e-08, 0.00e+00, 9.86e-08,
		2.25e-08, -1.67e-08, -1.62e-08, -9.42e-08, 1.17e-07, 4.32e-08, 3.97e-08,
		3.13e-08, -8.04e-08, 3.91e-08, 2.58e-08, 3.45e-08, 4.76e-08, 1.13e-07 },

		c1[20] = { -9.17e-05, -1.37e-05, 0.00e+00, 7.14e-05, 0.00e+00, -3.21e-04,
		-1.66e-04, -4.10e-05, 1.36e-04, 2.29e-04, -3.89e-04, -3.08e-04,
		-2.81e-04, -1.90e-04, 4.76e-05, -2.80e-04, -2.07e-04, -2.91e-04,
		-3.30e-04, -4.04e-04 },

		c0[20] = { 1.0136e+00, 1.0478e+00, 1.00e+00, 1.0258e+00, 1.00e+00,
		1.077e+00, 1.0543e+00, 1.0103e+00, 9.9927e-01, 9.6876e-01, 1.0971e+00,
		1.0971e+00, 1.0777e+00, 1.1134e+00, 1.0237e+00, 1.0703e+00, 1.0248e+00,
		1.0945e+00, 1.1622e+00, 1.1393e+00 },

		fap[36] = { 0., 0., 0.037037037, 0.074074074, 0.111111111, 0.148148148,
		0.185185185, 0.222222222, 0.259259259, 0.296296296, 0.333333333,
		0.37037037, 0.407407407, 0.444444444, 0.481481481, 0.518518519,
		0.555555556, 0.592592593, 0.62962963, 0.666666667, 0.703703704,
		0.740740741, 0.777777778, 0.814814815, 0.851851852, 0.888888889,
		0.925925926, 0.962962963, 1., 0.66666667, 0.33333334, 0., 0.333333,
		0.666666, 1.,0.7 };

	int code_tmp[] = { 3, 4, 5, 4, 3, 2, 1, 2, 3, 2, 1, 2, 3, 4, 5, 4, 8, 7, 6, 7, 8, 9,
		10, 9, 13, 12, 11, 12, 13, 14, 15, 14, 18, 17, 16, 17, 18, 19, 20, 19, 18, 17, 16,
		17, 18, 19, 20, 19 };
	int code[8][6];
	for (size_t i = 0; i < 48; i++)
		code[i % 8][i / 8] = code_tmp[i];
	//  CALLING THE PROGRAM TO CONVERT TO GEOMAGNETIC COORDINATES
	if (coor == 1)
		rgma = conver(rga, rgo);
	else if (coor == 2)
		rgma = rga;
	else 
		std::cout << "Wrong Coordinates Selection -------- >> " << coor << std::endl;
	int ape[39], i;

	// FROM 3-HOURLY TO HOURLY ap (New, interpolates between the three hourly ap values)
	ape[0] = ap[0];
	ape[1] = ap[0];
	ape[37] = ap[12];
	ape[38] = ap[12];

	for (size_t k = 1; k <= 13; k++) {
		i = (k * 3) - 1;
		ape[i - 1] = ap[k - 1];
	}
	for (size_t k = 1; k <= 12; k++) {
		i = k * 3;
		ape[i - 1] = (ap[k - 1] * 2 + ap[k]) / 3;
	}
	for (size_t k = 2; k <= 13; k++) {
		i = k * 3 - 2;
		ape[i - 1] = (ap[k - 2] + ap[k - 1] * 2) / 3;
	}
	/*C     TO OBTAIN THE INTEGRAL OF ap.
	C     INTEGRAL OF ap*/
	int k;
	if (ut == 24) ut = 0;
	if (ut % 3 == 0) k = 1;
	else if (ut % 3 == 1) k = 2;
	else if (ut % 3 == 2) k = 3;
	else {
		std::cout << "Wrong Universal Time value -------- >> " << ut << std::endl;
		return;
	}
	int rap = 0;
	for (size_t j = 1; j <= 36; j++)
		rap += fap[j - 1] * ape[k + j - 1];
	if (rap <= 200) {
		cf = 1.;
		return;
	}
	if ((doy > 366) || (doy < 1)) {
		std::cout << "Wrong Day of Year value --- >> " << doy << std::endl;
		return;
	}
	if ((rgma > 90) || (rgma < -90)) {
		std::cout << "Wrong GEOMAGNETIC LATITUDE value --- >> " << rgma << std::endl;
		return;
	}
	int dayno = doy;
	if (rgma < 0) {
		dayno = doy + 172;
		if (dayno > 365) dayno -= 365;
	}
	double rs;
	if (dayno >= 82) rs = (dayno - 82.) / 45.6 + 1;
	else rs = (dayno + 283.) / 45.6 + 1;
	int s1 = int(rs),
		s2 = s1 + 1;
	double facs = rs - s1;
	if (s2 == 9) s2 = 1;

	rgma = abs(rgma);

	double rl = (rgma + 10.) / 20. + 1;
	if (rl == 6.0) rl = 5.9;
	int l1 = int(rl),
		l2 = l1 + 1;
	double facl = rl - l1;

	// FACTORS CALCULATIONS
	double rapf, cf1, cf2, cf3, cf4, cf300;
	int n1, n2, n3, n4;
	if (rap < 300.) {
		rapf = 300.;
		n1 = code[s1 - 1][l1 - 1];
		cf1 = c4[n1 - 1] * pow(rapf, 4) + c3[n1 - 1] * pow(rapf, 3) +
			c2[n1 - 1] * pow(rapf, 2) + c1[n1 - 1] * rapf + c0[n1 - 1];
		n2 = code[s1 - 1][l2 - 1];
		cf2 = c4[n2 - 1] * pow(rapf, 4) + c3[n2 - 1] * pow(rapf, 3) +
			c2[n2 - 1] * pow(rapf, 2) + c1[n2 - 1] * rapf + c0[n2 - 1];
		n3 = code[s2 - 1][l1 - 1];
		cf3 = c4[n3 - 1] * pow(rapf, 4) + c3[n3 - 1] * pow(rapf, 3) +
			c2[n3 - 1] * pow(rapf, 2) + c1[n3 - 1] * rapf + c0[n3 - 1];
		n4 = code[s2 - 1][l2 - 1];
		cf4 = c4[n4 - 1] * pow(rapf, 4) + c3[n4 - 1] * pow(rapf, 3) +
			c2[n4 - 1] * pow(rapf, 2) + c1[n4 - 1] * rapf + c0[n4 - 1];

		//  INTERPOLATION
		cf300 = cf1 * (1 - facs) * (1 - facl) + cf2 * (1 - facs) * facl +
			cf3 * facs * (1 - facl) + cf4 * facs * facl;
		cf = (cf300 - 1.0) * rap / 100. - 2. * cf300 + 3.;
		return;
	}
	n1 = code[s1 - 1][l1 - 1];
	cf1 = c4[n1 - 1] * pow(rap, 4) + c3[n1 - 1] * pow(rap, 3) +
		c2[n1 - 1] * pow(rap, 2) + c1[n1 - 1] * rap + c0[n1 - 1];
	n2 = code[s1 - 1][l2 - 1];
	cf2 = c4[n2 - 1] * pow(rap, 4) + c3[n2 - 1] * pow(rap, 3) +
		c2[n2 - 1] * pow(rap, 2) + c1[n2 - 1] * rap + c0[n2 - 1];
	n3 = code[s2 - 1][l1 - 1];
	cf3 = c4[n3 - 1] * pow(rap, 4) + c3[n3 - 1] * pow(rap, 3) +
		c2[n3 - 1] * pow(rap, 2) + c1[n3 - 1] * rap + c0[n3 - 1];
	n4 = code[s2 - 1][l2 - 1];
	cf4 = c4[n4 - 1] * pow(rap, 4) + c3[n4 - 1] * pow(rap, 3) +
		c2[n4 - 1] * pow(rap, 2) + c1[n4 - 1] * rap + c0[n4 - 1];

	//  INTERPOLATION
	cf = cf1 * (1 - facs) * (1 - facl) + cf2 * (1 - facs) * facl +
		cf3 * facs * (1 - facl) + cf4 * facs * facl;
}

double conver(double rga, double rgo) {
/*C     This subroutine converts a geographic latitude and longitude
C     location to a corrected geomagnetic latitude.
C
C     INPUT: 
C       geographic latitude   -90. to +90.
C       geographic longitude  0. to 360. positive east from Greenwich.
C
C     OUTPUT:
C       corrected geomagnetic latitude	-90. to +90.*/
	double cormag[20][91];
	std::string filename = "data/conv.dat";
	std::fstream myfile(filename, std::ios_base::in);

	if (myfile.is_open())
	{
		for (size_t j = 0; j < 91; j++)
			for (size_t i = 0; i < 20; i++)
				myfile >> cormag[i][j];

		myfile.close();
	}
	else
		std::cout << "ERROR. FILE " << filename << "WASN'T OPENED!" << std::endl;

	double rlan = rga,
		rlo = rgo;
	/*C     From "normal" geographic latitude 
	C     to angle from South Pole. */
	double rla = rlan + 90;
	if (rlo == 360) rlo = 0;

	// PROXIMITY
	// coefficients of the latitudinal points
	int la1 = (int(rla / 2) + 1),
		la2 = la1 + 1;
	if (la2 > 91) la2 = 91;

	// coefficients of the longitudinal points
	int lo1 = (int(rlo / 18) + 1),
		lo2 = lo1 % 20 + 1;

	// Four points of Geomagnetic Coordinates
	double gm1 = cormag[lo1 - 1][la1 - 1],
		gm2 = cormag[lo1 - 1][la2 - 1],
		gm3 = cormag[lo2 - 1][la1 - 1],
		gm4 = cormag[lo2 - 1][la2 - 1];

	// latitudinal points		
	double x = rla / 2.0 - int(rla / 2.0);

	// longitudinal points
	double y = rlo / 18.0 - int(rlo / 18.0);

	// INTERPOLATION
	double gmla = gm1 * (1 - x) * (1 - y) + gm2 * (1 - y) * x + 
		gm3 * y * (1 - x) + gm4 * x * y;
	return 90 - gmla;
}

double fmodip(double xlat) {
	double dec, dip, dipl, ymodip;
	igrf_dip(xlat, findrlat.flon, findrlat.ryear, 300., dec, dip, dipl, ymodip);
	return ymodip;
}

void regfa1(double x11, double x22, double fx11, double fx22,
	double eps, double fw, double(*f)(double), bool schalt, double &x) {
/*C REGULA-FALSI-PROCEDURE TO FIND X WITH F(X)-FW=0. X1,X2 ARE THE                
C STARTING VALUES. THE COMUTATION ENDS WHEN THE X-INTERVAL                      
C HAS BECOME LESS THEN EPS . IF SIGN(F(X1)-FW)= SIGN(F(X2)-FW)                  
C THEN SCHALT=.TRUE.*/
	schalt = false;
	double ep = eps,
		x1 = x11,
		x2 = x22,
		f1 = fx11 - fw,
		f2 = fx22 - fw;
	int ng = 2,
		lfd = 0;
	if (f1 * f2 > 0) {
		x = 0.;
		schalt = true;
	}
	double fx, dx;
	bool l1, goto300 = false,
		links, k = false;
	while (true) {
		if (goto300) {
			goto300 = false;
			l1 = links;
			dx = (x2 - x1) / ng;
			if (!links) dx *= ng - 1;
			x = x1 + dx;
		}
		else
			x = (x1*f2 - x2 * f1) / (f2 - f1);
		fx = f(x) - fw;
		lfd++;
		if (lfd > 20) {
			ep *= 10.;
			lfd = 0;
		}
		links = (f1 * fx > 0.0);
		k = !k;
		if (links) {
			x1 = x;
			f1 = fx;
		}
		else {
			x2 = x;
			f2 = fx;
		}
		if (abs(x2 - x1) <= ep) break;
		if (k) {
			goto300 = true;
			continue;
		}
		if ((links && !l1) || (!links && l1))  ng *= 2;
	}
}

void shamdhmf2(double rlat, double flon, double t, double rz, double &hmf2) {
/*C-------------------------------------------------------------------
C	COMPUTES THE HOURLY VALUES OF hmF2 FROM A SET OF SH COEFFICIENTS
C	IN A POINT OF A GIVEN GEOCENTRIC LATITUDE AND LONGITUDE OF THE 
C	EARTH'S SURFACE FOR A GIVEN MONTH AND A GIVEN SUNSPOT NUMBER.
C   PARAMETERS AND COEFFICIENTS ARE GIVEN IN DATA STATEMENTS.
C
C INPUT:	RLAT    The geographic latitude on the FLON meridian  
C                   where the modified dip latitude is equal to the 
C                   modip value for the considered point in space.
C           FLON    =LONGITUDE+15.*UT(hours)
C           T       Month as a REAL number (1.0 to 12.0)
C           RZ      12-month running mean
C OUTPUT:   HMF2    F2 peak altitude in km
C
C Altadill, D., S. Magdaleno, J.M. Torta, and E. Blanch 
C     Adv. Space Res. 52, 1756-1769, 2013.
C-------------------------------------------------------------------*/
	double theta = 180.;
	amtb.re = 6371.2;
	amtb.tzero = 1.0;
	amtb.ifit = -1;
	int icen = 0, iref = 0;
	amtb.ib = 2;
	amtb.kint = 8;
	amtb.lint = 4;
	amtb.kext = 0;
	amtb.lext = -1;
	double const_[9][9],
		const_tmp[] = { 1., 1., 1., 1., 1., 1.73205, 0.866025, 1., 2.44949,
		1.93649, 0.7905691, 1., 3.16228, 3.35410, 2.09165, 0.739510, 1.,
		3.87298, 5.12348, 4.18330, 2.21853, 0.701561, 1., 4.58258, 7.24569,
		7.24569, 4.96078, 2.32681, 0.671693, 1., 5.29150, 9.72111, 11.4564,
		9.49918, 5.69951, 2.42182, 0.647260, 1., 6., 12.5499, 16.9926, 16.4531,
		11.8645, 6.40755, 2.50683, 0.626707 };
	size_t i = 0;
	for (size_t n = 0; n < 9; n++)
		for (size_t m = 0; m <= n; m++)	{
			const_[n][m] = const_tmp[i];
			i++;
		}
	std::string filename = "data/shamdhmf2.dat";
	std::fstream myfile(filename, std::ios_base::in);
	double ganm[9][9][5], gbnm[9][9][5], hanm[9][9][5],
		hbnm[9][9][5], gnm[9][9][5], hnm[9][9][5];
	if (myfile.is_open())
	{
		for (size_t n = 0; n < 9; n++)
			for (size_t m = 0; m <= n; m++)
				for (size_t j = 0; j < 5; j++)	{
					myfile >> ganm[n][m][j];
					myfile >> gbnm[n][m][j];
					myfile >> hanm[n][m][j];
					myfile >> hbnm[n][m][j];
				}
		myfile.close();
	}
	else
		std::cout << "ERROR. FILE " << filename << "WASN'T OPENED!" << std::endl;
	amtb.kmax = amtb.kint >= amtb.kext ? amtb.kint : amtb.kext;
	int kt, lj;
	if (amtb.kmax <= 8) {
		kt = amtb.lint >= amtb.lext ? amtb.lint : amtb.lext;
		if (kt <= 4) {
			for (size_t n = 0; n <= amtb.kmax; n++) {
				for (size_t m = 0; m <= n; m++) {
					for (size_t j = 0; j <= kt; j++) {
						gnm[n][m][j] = ganm[n][m][j] + gbnm[n][m][j] * rz;
						hnm[n][m][j] = hanm[n][m][j] + hbnm[n][m][j] * rz;
					}
					if (n <= amtb.kint) {
						lj = amtb.lint;
						amtb.fn[n][m] = n;
						if (m == 0) {
							for (size_t j = 0; j <= kt; j++)
								amtb.bint[n][m][j] = gnm[n][m][j];
						}
						else {
							for (size_t j = 0; j <= lj; j++) {
								amtb.bint[n][m][j] = gnm[n][m][j];
								amtb.bint[m - 1][n][j] = hnm[n][m][j];
							}
						}
					}	
				}
			}
			double dum;
			schnevpdh(rz, rlat, flon, dum, t, -1, dum, dum, hmf2);
		}
	}
}

void schnevpdh(double rz, double flat, double flon, double r, double t, int l, 
	double & bn, double & be, double & bv) {
/*C------------------------------------------------------------------------
C WHEN L IS POSITIVE:
C COMPUTES SPHERICAL CAP HARMONIC (GEOCENTRIC) FIELD COMPONENTS
C HORIZONTAL NORTH BN,HORIZONTAL EAST BE,AND VERTICAL DOWNWARD BV.
C WHEN L IS NEGATIVE:
C COMPUTES GENERAL FUNCTION BV, ITS HORIZONTAL NORTH DERIVATIVE BN,
C AND ITS HORIZONTAL EAST DERIVATIVE BE, ON SPHERICAL CAP SURFACE.
C NOTE THAT THESE ARE METRICAL DERIVATIVES, AND BE IS THE
C LONGITUDINAL DERIVATIVE DIVIDED BY SIN(COLATITUDE).
C	  INPUT:
C FLAT,FLON,R ARE GEOCENTRIC SPHERICAL CAP LATITUDE,LONGITUDE,RADIAL
C DISTANCE; T IS TIME.
C L =  0  ON FIRST CALL:  RETURNS SPHERICAL CAP POLE POSITION FLATO,FLONO
C         AND HALF-ANGLE THETA AS BN,BE, AND BV AFTER INITIALIZATION.
C         ON SUBSEQUENT CALLS:  ACTS AS L=1.
C      1  COMPUTES POTENTIAL FIELD COMPONENTS FROM INTERNAL COEFFICIENTS.
C      2  COMPUTES POTENTIAL FIELD COMPONENTS FROM EXTERNAL COEFFICIENTS.
C      3  COMPUTES FIELD FROM BOTH INTERNAL AND EXTERNAL COEFFICIENTS.
C     -1  COMPUTES GENERAL FUNCTION BV AND DERIVATIVES BN WITH RESPECT TO
C         LATITUDE AND BE WITH RESPECT TO LONGITUDE DIVIDED BY COS(LAT)
C         (R IS DUMMY VARIABLE IN THIS CASE).
C NOTE:   SUBROUTINE IS INITIALIZED DURING FIRST CALL REGARDLESS OF L.
C
C SUBROUTINE USED:  LEGFUN
C
C PARAMS & COEFFS TRANSFERRED FROM MAIN PROGRAM IN COMMON/AMTB/
C
C ADAPTED FROM SUBROUTINE SCHNEV OF G.V. HAINES (COMPUTERS & GEOSCIENCES, 
C 14, 413-447, 1988)
C------------------------------------------------------------------------*/
	double constp[9][9],
		const_tmp[] = { 1., 1., 1., 1., 1., 1.73205, 0.866025, 1., 2.44949,
		1.93649, 0.7905691, 1., 3.16228, 3.35410, 2.09165, 0.739510, 1.,
		3.87298, 5.12348, 4.18330, 2.21853, 0.701561, 1., 4.58258, 7.24569,
		7.24569, 4.96078, 2.32681, 0.671693, 1., 5.29150, 9.72111, 11.4564,
		9.49918, 5.69951, 2.42182, 0.647260, 1., 6., 12.5499, 16.9926, 16.4531,
		11.8645, 6.40755, 2.50683, 0.626707 };
	size_t i = 0;
	for (size_t n = 0; n < 9; n++)
		for (size_t m = 0; m <= n; m++) {
			constp[n][m] = const_tmp[i];
			i++;
		}
	/*C     IBF   =  0   TO USE ORDINARY POLYNOMIALS AS BASIS FUNCTIONS
	C              1          LEGENDRE POLYNOMIALS
	C              2          FOURIER SERIES
	C              3          COSINE SERIES
	C              4          SINE SERIES
	C     NOTE:    TZERO AND THINT MAY DEPEND ON IBF.*/
	int ibf = 2, kt, incept;
	double t1 = 1., t2 = 12., thint, aor, ar, aor3, delt[5], del, st, recip, ct;
	tbfit(t1, t2, ibf, thint, amtb.tzero);

	if (l >= 0) {
		if (amtb.ifit < 0) return;
		aor = amtb.re / r;
		ar = aor * aor;
		if (l > 1) {
			if (amtb.kext > 0) aor3 = aor * ar;
			if (l > 2) kt = amtb.lint >= amtb.lext ? amtb.lint : amtb.lext;
			else kt = amtb.lext;
		}
		else
			kt = amtb.lint;
	}
	else {
		if (amtb.ifit >= 0) return;
		ar = -1;
		kt = amtb.lint;
	}
	delt[0] = 1.;
	if (kt > 0) {
		del = (t - amtb.tzero) / thint;
		for (size_t i = 1; i <= kt; i++) {
			if (i == 1) {
				if (ibf <= 1)
					delt[i] = del;
				else if (ibf == 2) {
					st = sin(del);
					delt[i] = st;
				}
				else if (ibf == 3)
					delt[i] = cos(del);
				else
					delt[i] = sin(del);
				continue;
			}
			if (ibf == 0)
				delt[i] = delt[i - 1] * del;
			else if (ibf == 1) {
				recip = 1. / i;
				delt[i] = (2 - recip) * delt[i - 1] * del - (1 - recip)*delt[i - 2];
			}
			else if (ibf == 2)
				if (i % 2 == 0)
					if (i == 2) {
						ct = cos(del);
						delt[i] = ct;
					}
					else
						delt[i] = delt[i - 2] * ct - delt[i - 3] * st;
				else
					delt[i] = delt[i - 2] * ct - delt[i - 1] * st;
			else if (ibf == 3)
				delt[i] = cos(i * del);
			else if (ibf == 4)
				delt[i] = sin(i * del);
			else
				return;
		} //cycle 
		incept = 0;
			if (((ibf == 2) || (ibf == 3)) && incept == 1)
				for (size_t i = 2; i <= amtb.lint; i += 4 - ibf)
					delt[i] = 1 - delt[i];
	} // 103
	double x = 0, y = 0, ffn, p, dp, pms, gte, hti, hte, ra,
		z = 0, gti, colat, cl, sl, cml[8], sml[8], bt, bt1, bt2, bt3;
	int nmm;
	if ((l != 2) && (amtb.kint >= 0)) {
		gti = 0;
		for (size_t i = 0; i <= amtb.lint; i++)
			gti += amtb.bint[0][0][i] * delt[i];
		z = -ar * gti;
	}
	colat = 90 - flat;
	for (size_t n = 1; n <= amtb.kmax; n++)	{
		if (n == 1) {
			cl = cos(flon * constants::umr);
			sl = sin(flon * constants::umr);
			cml[0] = cl;
			sml[0] = sl;
		}
		else {
			sml[n - 1] = sl * cml[n - 2] + cl * sml[n - 2];
			cml[n - 1] = cl * cml[n - 2] - sl * sml[n - 2];
		}
		for (size_t m = 0; m <= n; m++)	{
			if (amtb.ib != 2) {
				nmm = n - m;
				if (nmm % 2) continue;
			}
			ffn = amtb.fn[n][m];
			legfun(m, ffn, constp[n][m], colat, p, dp, pms);
			if (l >= 0)
				ar = pow(aor, ffn + 2);
			else {
				ar = 1;
				ffn = -2;
				dp = -dp;
				pms = -pms;
			}
			if (!m) {
				bt1 = 0;
				bt3 = 0;
				bt = 0;
				if ((l != 2) && (n <= amtb.kint)) {
					gti = 0;
					for (size_t i = 0; i <= amtb.lint; i++)
						gti += amtb.bint[n][m][i];
					bt1 = ar * gti;
					bt3 = bt1;
				}
				if ((l < 1) && (n <= amtb.kext)) {
					gte = 0;
					for (size_t i = 0; i <= amtb.lext; i++)
						gte = gte + amtb.bext[n][m][i] * delt[i];
					bt = aor3 / ar * gte;
					bt1 = bt1 + bt;
				}
				x += bt1 * dp;
				z -= (ffn * (bt3 - bt) + bt3) * p;
				break;
			} // 130
			bt1 = 0.;
			bt2 = 0.;
			bt3 = 0.;
			bt = 0.;
			if ((l != 2) && (n <= amtb.kint)) {
				gti = 0;
				hti = 0;
				for (size_t i = 0; i <= amtb.lint; i++)	{
					gti += amtb.bint[n][m][i] * delt[i];
					hti += amtb.bint[m - 1][n][i] * delt[i];
				}
				bt1 = ar * (gti * cml[m - 1] + hti * sml[m - 1]);
				bt2 = ar * (gti * sml[m] - hti * cml[m]);
				bt3 = bt1;
			} // 133
			if ((l > 1) && (n <= amtb.kext)) {
				gte = 0;
				hte = 0;
				for (size_t i = 0; i <= amtb.lext; i++)	{
					gte += amtb.bext[n][m][i] * delt[i];
					hte += amtb.bext[m - 1][n][i] * delt[i];
				}
				ra = aor3 / ar;
				bt = ra * (gte * cml[m - 1] + hte * sml[m - 1]);
				bt1 = bt1 + bt;
				bt2 = bt2 + ra * (gte * sml[m - 1] - hte * cml[m - 1]);
			} // 135
			x += bt1 * dp;
			y += bt2 * pms;
			z -= (ffn * (bt3 - bt) + bt3) * p;
		}
	}
	bn = x;
	be = y;
	bv = z;
}

void tbfit(double t1, double t2, int ibf,
	double &thint, double &tzero) {
/*C-------------------------------------------------------------------
C	COURTESY OF G.V. HAINES
C
C     T2    =  BEGINNING OF TIME INTERVAL.
C     T1    =  END OF TIME INTERVAL.
C     IBF   =  0   TO USE ORDINARY POLYNOMIALS AS TEMPORAL BASIS FUNCTIONS
C              1          LEGENDRE POLYNOMIALS
C              2          FOURIER SERIES
C              3          COSINE SERIES
C              4          SINE SERIES
C              5          COSINE + SINE SERIES
C     TZERO =  TIME-TRANSLATION PARAMETER:
C              FOR IBF.LE.1, CHOOSE TZERO = CENTER OF TIME INTERVAL
C              FOR IBF.GE.2, CHOOSE TZERO = BEGINNING OF TIME INTERVAL.
C     THINT =  TIME-SCALING PARAMETER. THINT = HALF OF TIME INTERVAL T2-T1
C              FOR IBF.LE.1; PI*HALF OF TIME INTERVAL FOR IBF.EQ.2;
C              AND PI*TIME INTERVAL FOR IBF.GE.3.
C     NOTE:    CHOOSING TZERO AND THINT IN THIS WAY SCALES TIME
C              TO (-1,1) FOR IBF.LE.1;  TO (0,2PI) FOR IBF.EQ.2;
C              AND TO (0,PI) FOR IBF.GE.3.
C-------------------------------------------------------------------*/
	if (ibf <= 1)
		tzero = (t2 + t1) / 2.0;
	else
		tzero = t1;
	thint = t2 - t1;
	if (ibf <= 2)  thint = thint / 2.0;
}

void legfun(int m, double fn, double constp, double colat,
	double &p, double &dp, double &pms) {
/*C-------------------------------------------------------------------
C     SERIES FORM FOR ASSOCIATED LEGENDRE FUNCTION P, ITS DERIVATIVE DP,
C     AND THE FUNCTION PMS=P*M/SIN(COLAT), IN POWERS OF (1-COS(COLAT))/2.
C     INTEGRAL ORDER M, REAL DEGREE FN, NORMALIZING CONSTANT CONST.
C     COLATITUDE COLAT IN DEGREES.
C     IPRT = 0     NO PRINT-OUT
C            1     PRINT PARAMETERS AND P SERIES
C            2     PRINT PARAMETERS AND DP SERIES
C            3     PRINT PARAMETERS AND BOTH P AND DP SERIES
C           -1     PRINT PARAMETERS ONLY
C     INPUT M,FN,CONST,COLAT,IPRT.   OUTPUT P,DP,PMS
C	ADAPTED FROM G.V. HAINES (COMPUTERS & GEOSCIENCES, 14, 413-447, 1988).
C-------------------------------------------------------------------*/
	double fnn = fn * (fn + 1.), x, c, al;
	if (colat < 60) {
		x = sin(colat / 2 * constants::umr);
		c = 1 - 2 * x;
	}
	else {
		c = cos(colat * constants::umr);
		x = (1. - c) / 2.;
	}
	double s = sin(colat * constants::umr);
	if (m > 1)
		al = pow(constp * s, m - 1);
	else {
		if (m < 0) return;
		al = constp;
	}
	double pnm = al,
	dpnm = 0, b, a, absa, absb;
	int j = 0, jpm;
	do {
		j++;
		jpm = j + m;
		b = al * (jpm - 1 - fnn / jpm);
		dpnm = dpnm + b;
		a = (b * x) / j;
		pnm += a;
		al = a;
		absa = abs(a);
		absb = abs(b);
		if ((absa <= 1e-7) && (absb <= 1e-7))
			break;
		else if ((absb >= 1e+13) || (absa >= 1e+13))
			return; // error
	} while (j < 60);
	// SERIES TRUNCATED SUCCESSFULLY.
	double ps = pnm,
		dps = dpnm;
	if (m) {
		pms = ps * m;
		p = ps * s;
		dp = dps * s * s / 2 + c * pms;
	}
	else {
		pms = 0;
		p = ps;
		dp = dps * s / 2;
	}
}

double model_hmf2(int day, int month, double ut, double xmodip, double xlong,
	double f107_81) {
/*c--------------------------------------------------------------------------
c  Input:
c    day, month (integer)
c    UT         (real) - universal time
c    xmodip     (real) - modified dip latitude in degrees
c    long       (real) - geodatic longitude    in degrees
c    F107_81    (real) - F10.7 average over 3 solar rotations
c                        (81 days, centered on the current day)
c
c Output:
c    hmF2       (real) - F2-layer peak height in km
c--------------------------------------------------------------------------*/
	double hmf2_0 = sdmf2(ut, month, f107_81, xmodip, xlong),
		hmf2med, hmf2_m, hmf2_p;
	int monthr, montha;
	if (day <= 15)	{
		if (day == 15)
			hmf2med = hmf2_0;
		else {
			monthr = month - 1;
			if (monthr == 0) monthr = 12;
			hmf2_m = sdmf2(ut, monthr, f107_81, xmodip, xlong);
			hmf2med = (hmf2_0 - (day - 15) * (hmf2_m - hmf2_0) / 30.);
		}
	}
	else {
		montha = month % 12 + 1;
		hmf2_p = sdmf2(ut, montha, f107_81, xmodip, xlong);
		hmf2med = (hmf2_0 + (day - 15) * (hmf2_p - hmf2_0) / 30.);
	}

	return hmf2med;
}

double sdmf2(double ut, int monthut, double f107a, double xmodip, double xlong) {
/*c--------------------------------------------------------------------------
c    Global median model of the F2-layer peak height
c
c  Requires the following subroutines and functions:
c     SDMF2, hmF2_med_SD, read_data_SD, fun_hmF2_SD,
c     fun_Gk, Legendre, fun_hmF2UT, Koeff_UT, fun_Akp_UT,
c     fun_Fk_UT, fun_Gk_UT  
c
c Author of the code:
c         Valentin Shubin
c         Pushkov Institute of Terrestrial Magnetism,
c         Ionosphere and Radio wave propagation (IZMIRAN)
c         Moscow, Troitsk, 142190, Russia
c         e-mail: shubin@izmiran.ru
c         
c     [Ref. V.N. Shubin. Global median model of the F2-layer
c     peak height based on ionospheric radio-occultation and
c     ground-based Digisonde observations. Advances in Space Research (2015)
c     http://dx.doi.org/10.1016/j.asr.2015.05.029]
c
c  Input:
c    UT      - universal time (real)
c    monthut - month (integer)
c    xmodip  - modified dip latitude in degrees (real)
c    long    - geodatic longitude    in degrees (real)
c    F107A   - F10.7 index averaged over the 3 Sun rotations 
c              in units of 10^-22 W/(m^2 Hz) (real)
c
c  Output:
c    hmF2  - F2-layer peak height in km (real)
c--------------------------------------------------------------------------*/
	for (size_t i = 0; i < 24; i++)
		hmf2ut.hmf2_ut[i] = hmf2_med_sd(i, monthut, f107a, xmodip, xlong);

	return fun_hmf2ut(ut);
}

double hmf2_med_sd(int iut, int monthut, double f107a, double xmodip, double xlong) {
/*c---------------------------------------------------------------------
c    Input: 
c      iUT     - universal time (real)
c      monthut - month (integer)
c      F107A   - F10.7 index averaged over the 3 Sun rotations 
c                in units of 10^-22 W/(m^2 Hz) (real)
c      xmodip  - modified dip latitude in degrees (real)
c      long    - geodatic longitude    in degrees (real)
c
c    function to interpolate hmF2 between the two levels of solar activity
c    used the following auxiliary subroutines and functions:
c    read_data_SD, fun_hmF2_SD
c---------------------------------------------------------------------*/
/*c    Arrays ft1 (12) and ft2 (12) are the median values of F10.7A,
c    which were used as the margins for
c    the low and high solar activity levels */
	//                Jan   Feb   Mar   Apr   May   Jun
	double ft1[] = { 73.6, 72.3, 71.8, 70.9, 73.6, 73.0,
					 71.1, 69.4, 69.1, 70.9, 72.3, 74.6, };
	double ft2[] = {144.2,142.9,167.2,125.3,124.4,127.9,
					142.0,165.9,132.6,142.0,145.6,143.0, };
	//				  Jul   Aug   Sep   Oct   Nov   Dec

	double teta = 90.0 - xmodip,
		hmf2_1, hmf2_2;

	double **coeff_month = new double*[149];
	for (size_t i = 0; i < 149; i++)
		*coeff_month = new double[48];
	read_data_sd(monthut, coeff_month);

	double *kf = new double[149];
	for (size_t i = 0; i < 149; i++)
		kf[i] = coeff_month[i][iut];
	hmf2_1 = fun_hmf2_sd(teta, xlong, kf);
	for (size_t i = 0; i < 149; i++)
		kf[i] = coeff_month[i][iut + 24];
	hmf2_2 = fun_hmf2_sd(teta, xlong, kf);

	double cov = f107a,
		cov1 = ft1[monthut - 1],
		cov2 = ft2[monthut - 1],
		a = (hmf2_2 - hmf2_1) / log(cov2 / cov1),
		b = hmf2_2 - a * log(cov2);

	return a * log(cov) + b;
}

void read_data_sd(int month, double **coeff_month) {
	/*c------------------------------------------------------------------
	c    subroutine to read arrays mcsat11.datÖ mcsat22.dat
	c    with coefficients of hmF2 spatial decomposition
	c    for for 12 month, 24 UT hour and two solar activity levels
	c------------------------------------------------------------------*/
	std::string filename = "data/mcsat" + std::to_string(month + 10) + ".dat";
	std::fstream infile(filename, std::ios_base::in);
	std::string line, number_s;
	int k = 0,
		ii, jj;
	if (infile.is_open())
	{
		// Read each line into a string
		while (std::getline(infile, line))
			for (size_t i = 0; i < line.size(); i += 12) {
				number_s = line.substr(i, 12);
				number_s.replace(8, 1, "E");
				jj = k % 149;
				ii = (k - jj) / 149;
				coeff_month[jj][ii] = stod(number_s);
				k++;
			}
		infile.close();
	}
	else
		std::cout << "ERROR. FILE WASN'T OPENED!" << std::endl;
}

double fun_hmf2_sd(double teta, double xlong, double *kf) {
/*c---------------------------------------------------------------------
c    Input: 
c      teta - (90-modip) in degrees
c      long - geodatic longitude in degrees
c      Kf   - coefficients of hmF2 spatial decomposition 
c
c    function to calculate spherical harmonics decomposition
c    for the spatial dependencies of hmF2
c    used the following auxiliary subroutines and functions:
c    fun_Gk, Legendre
c---------------------------------------------------------------------*/
	double *gk = new double[149],
		hmf2 = 0.;
	fun_gk(teta, xlong, gk);
	for (size_t k = 0; k < 149; k++)
		hmf2 += kf[k] * gk[k];

	return hmf2;
}

void fun_gk(double teta, double xlong, double *gk) {
	double **pl_mn = new double*[9];
	for (size_t i = 0; i < 9; i++)
		pl_mn[i] = new double[13];
	int mm = 8,
		nn = 12,
		k = 0;
	legendre(mm, nn, pl_mn, teta);
	for (size_t m = 0; m <= mm; m++)
		if (m == 0)
			for (size_t n = 0; n <= nn; n++) {
				gk[k] = pl_mn[m][n];
				k++;
			}
		else
			for (size_t n = m; n <= nn; n++) {
				gk[k] = pl_mn[m][n] * cos(m * xlong * constants::umr);
				gk[k + 1] = pl_mn[m][n] * sin(m * xlong * constants::umr);
				k += 2;
			}
}

void legendre(int mm, int nn, double **p, double teta) {
/*c---------------------------------------------------------------------
c     Input:
c       mm       - harmonics for longitude
c       nn       - harmonics for the modified dip latitude (modip) 
c       teta     - (90-modip) in degrees
c    Output:
c       P(mm,nn) - associated Legendre function    
c
c    subroutine to calculate associated Legendre function P(mm,nn)
c    with Schmidt normalization
c---------------------------------------------------------------------*/
	double z = cos(constants::umr * teta);
	p[0][0] = 1;
	p[0][1] = z;
	if (mm) p[1][1] = sin(constants::umr * teta);

	for (size_t j = 2; j <= mm; j++)
		p[j][j] = (2 * j - 1) * p[j - 1][j - 1] * p[1][1];

	for (size_t m = 0; m <= mm; m++)
		for (size_t n = 1; n <= nn; n++) {
			if (m > n)	{
				p[m][n] = 0.0;
				continue;
			}
			if (n + 1 > nn) break;
			if (n + 1 == m) continue;
			if (m > (n - 1))
				p[m][n + 1] = (2 * n + 1) * z * p[m][n] / (n + 1 - m);
			else
				p[m][n + 1] = ((2 * n + 1) * z * p[m][n] - (n + m) * p[m][n - 1]) / (n + 1 - m);
		}
	int s;
	for (size_t n = 1; n <= nn; n++) {
		for (size_t m = 1; m <= mm; m++) {
			if (m > n) {
				p[m][n] = 0.0;
				break;
			}
			s = 1;
			for (size_t l = n - m + 1; l <= n + m; l++)
				s *= l;
			p[m][n] *= sqrt(2.0 / s);
		}
	}
}

double fun_hmf2ut(double t) {
/*c---------------------------------------------------------------------
c    Input: T - universal time        
c 
c    function to calculate Fourier  decomposition
c    for the temporal variations of hmF2
c    used the following auxiliary subroutines:
c    Koeff_UT, fun_Akp_UT, fun_Fk_UT, fun_Gk_UT
c---------------------------------------------------------------------*/
	int mm = 3,
		mk = mm * 2;
	double *gk_ut = new double[7],
		*kf_ut = new double[7];
	koeff_ut(mm, mk, kf_ut);
	fun_gk_ut(mm, mk, t, gk_ut);
	double hmf2 = 0.;
	for (size_t k = 0; k <= mk; k++)
		hmf2 += kf_ut[k] * gk_ut[k];
	return hmf2;
}

void koeff_ut(int mm, int mk, double *kf_ut) {
	double **akp_ut = new double*[mk + 1],
		*dk_ut = new double[mk + 1];
	for (size_t i = 0; i <= mk; i++)
		akp_ut[i] = new double[mk + 1];
	fun_akp_ut(mm, mk, akp_ut, dk_ut);
	double sum_d;
	for (size_t k = mk; k >= 0; k--) {
		sum_d = 0;
		for (size_t m = k+1; m <= mk; m++)
			sum_d += akp_ut[m][k] * kf_ut[m];
		kf_ut[k] = sum_d + dk_ut[k];
	}
}

void fun_akp_ut(int mm, int mk, double **akp_ut, double *dk_ut) {
	double *gk_ut = new double[mk + 1],
		*fk_ut = new double[mk + 1];
	gk_ut[0] = 1;
	fk_ut[0] = 1;
	double sum_dn, sum_dd, sum_an, sum_ad;
	for (size_t p = 0; p <= mk; p++) {
		sum_dn = 0;
		sum_dd = 0;
		for (size_t k = p + 1; k <= mk; k++) {
			sum_an = 0;
			sum_ad = 0;
			for (size_t i = 0; i <= 23; i++) {
				fun_gk_ut(mm, mk, i, gk_ut);
				fun_fk_ut(mk, gk_ut, akp_ut, fk_ut);
				sum_an += gk_ut[k] * fk_ut[p];
				sum_ad += fk_ut[p] * fk_ut[p];
				if (p == (k - 1)) {
					sum_dn += hmf2ut.hmf2_ut[i] * fk_ut[p];
					sum_dd += fk_ut[p] * fk_ut[p];
				}
			}
			akp_ut[k][p] = -sum_an / sum_ad;
		}
		if (p < mk)
			dk_ut[p] = sum_dn / sum_dd;
	}
	int p = mk;
	sum_dn = 0;
	sum_dd = 0;
	for (size_t i = 0; i <= 23; i++) {
		fun_gk_ut(mm, mk, i, gk_ut);
		fun_fk_ut(mk, gk_ut, akp_ut, fk_ut);
		sum_dn += hmf2ut.hmf2_ut[i] * fk_ut[p];
		sum_ad += fk_ut[p] * fk_ut[p];
	}
	dk_ut[p] = sum_dn / sum_dd;
}

void fun_fk_ut(int mk, double *gk_ut, double **akp_ut, double *fk_ut) {
	double sum_g;
	for (size_t k = 0; k <= mk; k++) {
		sum_g = 0;
		for (size_t p = 0; p <= k; p++)	{
			if (k == p) continue;
			sum_g += akp_ut[k][p] * fk_ut[p];
		}
		fk_ut[k] = sum_g + gk_ut[k];
 	}
}

void fun_gk_ut(int mm, int mk, double t, double *gk_ut) {
	size_t k = 0;
	for (size_t m = 0; m <= mm; m++) {
		if (m == 0) {
			gk_ut[k] = 1;
			k++;
		}
		else {
			gk_ut[k] = cos(m * t * constants::humr);
			gk_ut[k + 1] = sin(m * t * constants::humr);
			k += 2;
		}
	}
}

double b0_98(double hour, double sax, double sux, int nseasn, double r,
	double zlo, double zmodip) {
/*C-----------------------------------------------------------------
C Interpolation procedure for bottomside thickness parameter B0.
C Array B0F(ILT,ISEASON,IR,ILATI) distinguishes between day and
C night (ILT=1,2), four seasons (ISEASON is northern season with
C ISEASON=1 northern spring), low and high solar activity Rz12=10,
C 100 (IR=1,2), and modified dip latitudes of 0, 18 and 45
C degress (ILATI=1,2,3). In the DATA statement the first value
C corresponds to B0F(1,1,1,1), the second to B0F(2,1,1,1), the
C third to B0F(1,2,1,1) and so on.
C
C input:
C       hour    LT in decimal hours
C       SAX     time of sunrise in decimal hours
C       SUX     time of sunset in decimal hours
C       nseasn  season in northern hemisphere (1=spring)
C       R       12-month running mean of sunspot number
C       ZLO     longitude
C       ZMODIP  modified dip latitude
C
C JUNE 1989 --------------------------------------- Dieter Bilitza
C
C Updates (B0_new -> B0_98):
C
C 01/98 corrected to include a smooth transition at the modip equator
C       and no discontinuity at the equatorial change in season.
C 09/98 new B0 values incl values at the magnetic equator
C 10/98 longitude as input to determine if magnetic equator in northern 
C         or southern hemisphere
C-------------------------------------------------------------------*/
	// jseasn is southern hemisphere season
	int jseasn = nseasn + 2,
		num_lat = 3;
	if (jseasn > 4) jseasn -= 4;
	double zz = zmodip + 90,
		zz0 = 0,
		bfr[2][2][3],
		bfd[2][3];

	int b0f_[] = {	
		201, 68, 210, 61, 192, 68, 199, 67, 240, 80, 245, 83, 
		233, 71, 230, 65, 108, 65, 142, 81, 110, 68, 77, 75, 
		124, 98, 164, 100, 120, 94, 96, 112, 78, 81, 94, 84, 
		81, 81, 65, 70, 102, 87, 127, 91, 109, 88, 81, 78 
	};

	int m = 0;
	double b0f[2][4][2][3];
	for (size_t i = 0; i < 3; i++)
		for (size_t j = 0; j < 2; j++)
			for (size_t k = 0; k < 4; k++)
				for (size_t l = 0; l < 2; l++) {
					b0f[l][k][j][i] = b0f_[m];
					m++;
				}
	double dayval, nitval;
	// Interpolation in Rz12: linear from 10 to 100
	for (size_t isl = 1; isl <= num_lat; isl++)	{
		for (size_t isd = 1; isd <= 2; isd++) {
			bfr[isd - 1][0][isl - 1] = b0f[isd - 1][nseasn - 1][0][isl - 1] +
				(b0f[isd - 1][nseasn - 1][1][isl - 1] - b0f[isd - 1][nseasn - 1][0][isl - 1]) /
				90 * (r - 10);
			bfr[isd - 1][1][isl - 1] = b0f[isd - 1][jseasn - 1][0][isl - 1] +
				(b0f[isd - 1][jseasn - 1][1][isl - 1] - b0f[isd - 1][jseasn - 1][0][isl - 1]) /
				90 * (r - 10);
		}

		/*Interpolation day/night with transitions at SAX (sunrise)
		and SUX (sunset) for northern/southern hemisphere iss=1/2*/
		for (size_t iss = 1; iss <= 2; iss++) {
			dayval = bfr[0][iss - 1][isl - 1];
			nitval = bfr[1][iss - 1][isl - 1];
			bfd[iss - 1][isl - 1] = hpol(hour, dayval, nitval, sax, sux, 1., 1.);
		}
	}
	/*C Interpolation with epstein-transitions in modified dip latitude.
	C Transitions at +/-18 and +/-45 degrees; constant above +/-45.
	C
	C g(1:5) are the latitudinal slopes of B0;
	C       g(1) is for the region from -90 to -45 degrees
	C       g(2) is for the region from -45 to -18 degrees
	C       g(3) is for the region from -18 to   0 degrees
	C       g(4) is for the region from   0 to  18 degrees
	C       g(5) is for the region from  18 to  45 degrees
	C       g(6) is for the region from  45 to  90 degrees
	C
	C B0 =  bfd(2,3) at modip = -45,
	C       bfd(2,2) at modip = -18,
	C       bfd(2,1) or bfd(1,1) at modip = 0,
	C       bfd(1,2) at modip = 20,
	C       bfd(1,3) at modip = 45.
	C If the Longitude is between 200 and 320 degrees then the modip 
	C equator is in the southern hemisphere and bfd(2,1) is used at the 
	C equator, otherwise bfd(1,1) is used.*/
	double zx1 = bfd[1][2],
		zx2 = bfd[1][1],
		zx3 = bfd[0][0];
	if ((zlo > 200) && (zlo < 320)) zx3 = bfd[1][0];
	double zx4 = bfd[0][1],
		zx5 = bfd[0][2],
		g[] = {
		0,
		(zx2 - zx1) / 27,
		(zx3 - zx2) / 18,
		(zx4 - zx3) / 18,
		(zx5 - zx4) / 27,
		0
	};
	double sum = zx1, aa, bb, dsum,
		zx[] = { 45.,72.,90.,108.,135. },
		dd[] = { 3.0, 3.0, 3.0, 3.0, 3.0 };
	for (size_t i = 1; i <= 5; i++) {
		aa = eptr(zz, dd[i-1], zx[i-1]);
		bb = eptr(zz0, dd[i - 1], zx[i - 1]);
		dsum = (g[i] - g[i - 1]) * (aa - bb) * dd[i - 1];
		sum += dsum;
	}
	return sum;
}

double eptr(double x, double sc, double hx) {
	double d1 = (x - hx) / sc,
		res;
	if (abs(d1) < constants::argmax)
		res = log(1 + exp(d1));
	else if (d1 > 0)
		res = d1;
	else
		res = 0;

	return res;
}

void shamdb0d(double rlat, double flon, double t, double rz, double &b) {
/*C-------------------------------------------------------------------
C	COMPUTES THE HOURLY VALUES OF B0 FROM A SET OF SH COEFFICIENTS
C	IN A POINT OF A GIVEN GEOCENTRIC LATITUDE AND LONGITUDE
C	OF THE EARTH'S SURFACE FOR A GIVEN MONTH AND A GIVEN SUSPOT NUMER
C
C INPUT:	RLAT    The geogrphic latitude on the meridian given by 
C					the local time (FLON), where the modified dip
C                   latitude is the same as of the orginal site.
C			FLON	=LONGITUDE+15.*UT(hours)
C			T		Month as a REAL number (1.0 to 12.0)
C			RZ		12-month running mean
C OUTOUT	B		=B0
C
C  Blanch E., D. Arrazola, D. Altadill, D. Buresova, M. Mosert, 
C     Adv. Space Res. 39, 701-710, 2007.
C  Altadill, D., D. Arrazola, E. Blanch, D. Buresova, 
C     Adv. Space Res. 42, 610-616, 2008.
C  Altadill, D., J.M. Torta, and E. Blanch, 
C     Adv. Space Res. 43,1825-1834, 2009.
C-------------------------------------------------------------------*/
	const size_t kdim = 6, ldim = 4;
	double theta = 180.;
	amtb.re = 6371.2;
	amtb.tzero = 1.0;
	amtb.ifit = -1;
	int icen = 0, iref = 0;
	amtb.ib = 2;
	amtb.kint = 6;
	amtb.lint = 4;
	amtb.kext = 0;
	amtb.lext = -1;
	double const_[kdim + 1][kdim + 1],
		const_tmp[] = { 1., 1., 1., 1., 1., 1.73205, 0.866025, 1., 2.44949,
		1.93649, 0.7905691, 1., 3.16228, 3.35410, 2.09165, 0.739510, 1.,
		3.87298, 5.12348, 4.18330, 2.21853, 0.701561, 1., 4.58258, 7.24569,
		7.24569, 4.96078, 2.32681, 0.671693 };
	size_t i = 0;
	for (size_t n = 0; n <= kdim; n++)
		for (size_t m = 0; m <= n; m++) {
			const_[n][m] = const_tmp[i];
			i++;
		}
	std::string filename = "data/shamdb0d.dat";
	std::fstream myfile(filename, std::ios_base::in);
	double ganm[kdim + 1][kdim + 1][ldim + 1], gbnm[kdim + 1][kdim + 1][ldim + 1], 
		hanm[kdim + 1][kdim + 1][ldim + 1], hbnm[kdim + 1][kdim + 1][ldim + 1],
		gnm[kdim + 1][kdim + 1][ldim + 1], hnm[kdim + 1][kdim + 1][ldim + 1];
	if (myfile.is_open())
	{
		for (size_t n = 0; n <= kdim; n++)
			for (size_t m = 0; m <= n; m++)
				for (size_t j = 0; j <= ldim; j++) {
					myfile >> ganm[n][m][j];
					myfile >> gbnm[n][m][j];
					myfile >> hanm[n][m][j];
					myfile >> hbnm[n][m][j];
				}
		myfile.close();
	}
	else
		std::cout << "ERROR. FILE " << filename << "WASN'T OPENED!" << std::endl;
	amtb.kmax = amtb.kint >= amtb.kext ? amtb.kint : amtb.kext;
	int kt = amtb.lint >= amtb.lext ? amtb.lint : amtb.lext,
		lj;

	for (size_t n = 0; n <= amtb.kmax; n++)
		for (size_t m = 0; m <= n; m++) {
			for (size_t j = 0; j <= kt; j++) {
				gnm[n][m][j] = ganm[n][m][j] + gbnm[n][m][j] * rz;
				hnm[n][m][j] = hanm[n][m][j] + hbnm[n][m][j] * rz;
			}
			if (n <= amtb.kint) {
				lj = amtb.lint;
				amtb.fn[n][m] = n;
				if (m == 0)
					for (size_t j = 0; j <= kt; j++)
						amtb.bint[n][m][j] = gnm[n][m][j];
				else
					for (size_t j = 0; j <= lj; j++) {
						amtb.bint[n][m][j] = gnm[n][m][j];
						amtb.bint[m - 1][n][j] = hnm[n][m][j];
					}
			}
		}
	double dum;
	schnevpd(rz, rlat, flon, dum, t, -1, dum, dum, b);
}

void schnevpd(double rz, double flat, double flon, double r, double t,
	int l, double &bn, double &be, double &bv) {
/*C-------------------------------------------------------------------
C     WHEN L IS POSITIVE:
C     COMPUTES SPHERICAL CAP HARMONIC (GEOCENTRIC) FIELD COMPONENTS
C     HORIZONTAL NORTH BN,HORIZONTAL EAST BE,AND VERTICAL DOWNWARD BV.
C     WHEN L IS NEGATIVE:
C     COMPUTES GENERAL FUNCTION BV, ITS HORIZONTAL NORTH DERIVATIVE BN,
C     AND ITS HORIZONTAL EAST DERIVATIVE BE, ON SPHERICAL CAP SURFACE.
C     NOTE THAT THESE ARE METRICAL DERIVATIVES, AND BE IS THE
C     LONGITUDINAL DERIVATIVE DIVIDED BY SIN(COLATITUDE).

C     FLAT,FLON,R ARE GEOCENTRIC SPHERICAL CAP LATITUDE,LONGITUDE,RADIAL
C     DISTANCE; T IS TIME.

C     L =  0  ON FIRST CALL:  RETURNS SPHERICAL CAP POLE POSITION FLATO,FLONO
C             AND HALF-ANGLE THETA AS BN,BE, AND BV AFTER INITIALIZATION.
C             ON SUBSEQUENT CALLS:  ACTS AS L=1.
C          1  COMPUTES POTENTIAL FIELD COMPONENTS FROM INTERNAL COEFFICIENTS.
C          2  COMPUTES POTENTIAL FIELD COMPONENTS FROM EXTERNAL COEFFICIENTS.
C          3  COMPUTES FIELD FROM BOTH INTERNAL AND EXTERNAL COEFFICIENTS.
C         -1  COMPUTES GENERAL FUNCTION BV AND DERIVATIVES BN WITH RESPECT TO
C             LATITUDE AND BE WITH RESPECT TO LONGITUDE DIVIDED BY COS(LAT)
C             (R IS DUMMY VARIABLE IN THIS CASE).
C     NOTE:   SUBROUTINE IS INITIALIZED DURING FIRST CALL REGARDLESS OF L.

C     SUBPROGRAM USED:  LEGFUN

C	***** PARAMS & COEFFS TRANSFERRED FROM MAIN PROGRAM *****

C	ADAPTED FROM SUBROUTINE SCHNEV OF G.V. HAINES (COMPUTERS & GEOSCIENCES, 
C      14, 413-447, 1988)
C-------------------------------------------------------------------*/
	const size_t kdim = 6, ldim = 4;
	double constp[kdim + 1][kdim + 1],
		const_tmp[] = { 1., 1., 1., 1., 1., 1.73205, 0.866025, 1., 2.44949,
		1.93649, 0.7905691, 1., 3.16228, 3.35410, 2.09165, 0.739510, 1.,
		3.87298, 5.12348, 4.18330, 2.21853, 0.701561, 1., 4.58258, 7.24569,
		7.24569, 4.96078, 2.32681, 0.671693 };
	size_t i = 0;
	for (size_t n = 0; n <= kdim; n++)
		for (size_t m = 0; m <= n; m++) {
			constp[n][m] = const_tmp[i];
			i++;
		}
	/*C     IBF   =  0   TO USE ORDINARY POLYNOMIALS AS BASIS FUNCTIONS
	C              1          LEGENDRE POLYNOMIALS
	C              2          FOURIER SERIES
	C              3          COSINE SERIES
	C              4          SINE SERIES
	C     NOTE:    TZERO AND THINT MAY DEPEND ON IBF.*/
	int ibf = 2, kt, incept;
	double t1 = 1., t2 = 12., thint, aor, ar, aor3, delt[5], del, st, recip, ct;
	tbfit(t1, t2, ibf, thint, amtb.tzero);

	if (l >= 0) {
		if (amtb.ifit < 0) return;
		aor = amtb.re / r;
		ar = aor * aor;
		if (l > 1) {
			if (amtb.kext > 0) aor3 = aor * ar;
			if (l > 2) kt = amtb.lint >= amtb.lext ? amtb.lint : amtb.lext;
			else kt = amtb.lext;
		}
		else
			kt = amtb.lint;
	}
	else {
		if (amtb.ifit >= 0) return;
		ar = -1;
		kt = amtb.lint;
	}
	delt[0] = 1.;
	if (kt > 0) {
		del = (t - amtb.tzero) / thint;
		for (size_t i = 1; i <= kt; i++) {
			if (i == 1) {
				if (ibf <= 1)
					delt[i] = del;
				else if (ibf == 2) {
					st = sin(del);
					delt[i] = st;
				}
				else if (ibf == 3)
					delt[i] = cos(del);
				else
					delt[i] = sin(del);
				continue;
			}
			if (ibf == 0)
				delt[i] = delt[i - 1] * del;
			else if (ibf == 1) {
				recip = 1. / i;
				delt[i] = (2 - recip) * delt[i - 1] * del - (1 - recip)*delt[i - 2];
			}
			else if (ibf == 2)
				if (i % 2 == 0)
					if (i == 2) {
						ct = cos(del);
						delt[i] = ct;
					}
					else
						delt[i] = delt[i - 2] * ct - delt[i - 3] * st;
				else
					delt[i] = delt[i - 2] * ct - delt[i - 1] * st;
			else if (ibf == 3)
				delt[i] = cos(i * del);
			else if (ibf == 4)
				delt[i] = sin(i * del);
			else
				return;
		} //cycle 
		incept = 0;
		if (((ibf == 2) || (ibf == 3)) && incept == 1)
			for (size_t i = 2; i <= amtb.lint; i += 4 - ibf)
				delt[i] = 1 - delt[i];
	} // 103
	double x = 0, y = 0, ffn, p, dp, pms, gte, hti, hte, ra,
		z = 0, gti, colat, cl, sl, cml[8], sml[8], bt, bt1, bt2, bt3;
	int nmm;
	if ((l != 2) && (amtb.kint >= 0)) {
		gti = 0;
		for (size_t i = 0; i <= amtb.lint; i++)
			gti += amtb.bint[0][0][i] * delt[i];
		z = -ar * gti;
	}
	colat = 90 - flat;
	for (size_t n = 1; n <= amtb.kmax; n++) {
		if (n == 1) {
			cl = cos(flon * constants::umr);
			sl = sin(flon * constants::umr);
			cml[0] = cl;
			sml[0] = sl;
		}
		else {
			sml[n - 1] = sl * cml[n - 2] + cl * sml[n - 2];
			cml[n - 1] = cl * cml[n - 2] - sl * sml[n - 2];
		}
		for (size_t m = 0; m <= n; m++) {
			if (amtb.ib != 2) {
				nmm = n - m;
				if (nmm % 2) continue;
			}
			ffn = amtb.fn[n][m];
			legfun(m, ffn, constp[n][m], colat, p, dp, pms);
			if (l >= 0)
				ar = pow(aor, ffn + 2);
			else {
				ar = 1;
				ffn = -2;
				dp = -dp;
				pms = -pms;
			}
			if (!m) {
				bt1 = 0;
				bt3 = 0;
				bt = 0;
				if ((l != 2) && (n <= amtb.kint)) {
					gti = 0;
					for (size_t i = 0; i <= amtb.lint; i++)
						gti += amtb.bint[n][m][i];
					bt1 = ar * gti;
					bt3 = bt1;
				}
				if ((l < 1) && (n <= amtb.kext)) {
					gte = 0;
					for (size_t i = 0; i <= amtb.lext; i++)
						gte = gte + amtb.bext[n][m][i] * delt[i];
					bt = aor3 / ar * gte;
					bt1 = bt1 + bt;
				}
				x += bt1 * dp;
				z -= (ffn * (bt3 - bt) + bt3) * p;
				break;
			} // 130
			bt1 = 0.;
			bt2 = 0.;
			bt3 = 0.;
			bt = 0.;
			if ((l != 2) && (n <= amtb.kint)) {
				gti = 0;
				hti = 0;
				for (size_t i = 0; i <= amtb.lint; i++) {
					gti += amtb.bint[n][m][i] * delt[i];
					hti += amtb.bint[m - 1][n][i] * delt[i];
				}
				bt1 = ar * (gti * cml[m - 1] + hti * sml[m - 1]);
				bt2 = ar * (gti * sml[m] - hti * cml[m]);
				bt3 = bt1;
			} // 133
			if ((l > 1) && (n <= amtb.kext)) {
				gte = 0;
				hte = 0;
				for (size_t i = 0; i <= amtb.lext; i++) {
					gte += amtb.bext[n][m][i] * delt[i];
					hte += amtb.bext[m - 1][n][i] * delt[i];
				}
				ra = aor3 / ar;
				bt = ra * (gte * cml[m - 1] + hte * sml[m - 1]);
				bt1 = bt1 + bt;
				bt2 = bt2 + ra * (gte * sml[m - 1] - hte * cml[m - 1]);
			} // 135
			x += bt1 * dp;
			y += bt2 * pms;
			z -= (ffn * (bt3 - bt) + bt3) * p;
		}
	}
	bn = x;
	be = y;
	bv = z;
}

void shab1d(double flat, double flon, double t, double rz, double &b) {
/*C-------------------------------------------------------------------
C	COMPUTES THE HOURLY VALUES OF B1 FROM A SET OF SH COEFFICIENTS
C	IN A POINT OF A GIVEN GEOCENTRIC LATITUDE AND LONGITUDE
C	OF THE EARTH'S SURFACE FOR A GIVEN MONTH AND A GIVEN SUSPOT NUMER
C
C   PARAMETERS ARE THE SAME AS IN SHAMDB0D, EXCEPT:
C		FLAT	Geographic latitude
C		B		=B1
C
C	***** PARAMS & COEFFS IN DATA SENTENCES *****
C-------------------------------------------------------------------*/
	const size_t kdim = 6, ldim = 4;
	double theta = 180.;
	amtb.re = 6371.2;
	amtb.tzero = 1.0;
	amtb.ifit = -1;
	int icen = 0, iref = 0;
	amtb.ib = 2;
	amtb.kint = 6;
	amtb.lint = 4;
	amtb.kext = 0;
	amtb.lext = -1;
	double const_[kdim + 1][kdim + 1],
		const_tmp[] = { 1., 1., 1., 1., 1., 1.73205, 0.866025, 1., 2.44949,
		1.93649, 0.7905691, 1., 3.16228, 3.35410, 2.09165, 0.739510, 1.,
		3.87298, 5.12348, 4.18330, 2.21853, 0.701561, 1., 4.58258, 7.24569,
		7.24569, 4.96078, 2.32681, 0.671693 };
	size_t i = 0;
	for (size_t n = 0; n <= kdim; n++)
		for (size_t m = 0; m <= n; m++) {
			const_[n][m] = const_tmp[i];
			i++;
		}
	std::string filename = "data/shab1d.dat";
	std::fstream myfile(filename, std::ios_base::in);
	double ganm[kdim + 1][kdim + 1][ldim + 1], gbnm[kdim + 1][kdim + 1][ldim + 1],
		hanm[kdim + 1][kdim + 1][ldim + 1], hbnm[kdim + 1][kdim + 1][ldim + 1],
		gnm[kdim + 1][kdim + 1][ldim + 1], hnm[kdim + 1][kdim + 1][ldim + 1];
	if (myfile.is_open())
	{
		for (size_t n = 0; n <= kdim; n++)
			for (size_t m = 0; m <= n; m++)
				for (size_t j = 0; j <= ldim; j++) {
					myfile >> ganm[n][m][j];
					myfile >> gbnm[n][m][j];
					myfile >> hanm[n][m][j];
					myfile >> hbnm[n][m][j];
				}
		myfile.close();
	}
	else
		std::cout << "ERROR. FILE " << filename << "WASN'T OPENED!" << std::endl;
	amtb.kmax = amtb.kint >= amtb.kext ? amtb.kint : amtb.kext;
	int kt = amtb.lint >= amtb.lext ? amtb.lint : amtb.lext,
		lj;

	for (size_t n = 0; n <= amtb.kmax; n++)
		for (size_t m = 0; m <= n; m++) {
			for (size_t j = 0; j <= kt; j++) {
				gnm[n][m][j] = ganm[n][m][j] + gbnm[n][m][j] * rz;
				hnm[n][m][j] = hanm[n][m][j] + hbnm[n][m][j] * rz;
			}
			if (n <= amtb.kint) {
				lj = amtb.lint;
				amtb.fn[n][m] = n;
				if (m == 0)
					for (size_t j = 0; j <= kt; j++)
						amtb.bint[n][m][j] = gnm[n][m][j];
				else
					for (size_t j = 0; j <= lj; j++) {
						amtb.bint[n][m][j] = gnm[n][m][j];
						amtb.bint[m - 1][n][j] = hnm[n][m][j];
					}
			}
		}
	double dum;
	schnevpdb1(rz, flat, flon, dum, t, -1, dum, dum, b);
}

void schnevpdb1(double rz, double flat, double flon, double r, double t,
	int l, double &bn, double &be, double &bv) {
/*C-------------------------------------------------------------------
C     WHEN L IS POSITIVE:
C     COMPUTES SPHERICAL CAP HARMONIC (GEOCENTRIC) FIELD COMPONENTS
C     HORIZONTAL NORTH BN,HORIZONTAL EAST BE,AND VERTICAL DOWNWARD BV.
C     WHEN L IS NEGATIVE:
C     COMPUTES GENERAL FUNCTION BV, ITS HORIZONTAL NORTH DERIVATIVE BN,
C     AND ITS HORIZONTAL EAST DERIVATIVE BE, ON SPHERICAL CAP SURFACE.
C     NOTE THAT THESE ARE METRICAL DERIVATIVES, AND BE IS THE
C     LONGITUDINAL DERIVATIVE DIVIDED BY SIN(COLATITUDE).

C     FLAT,FLON,R ARE GEOCENTRIC SPHERICAL CAP LATITUDE,LONGITUDE,RADIAL
C     DISTANCE; T IS TIME.

C     L =  0  ON FIRST CALL:  RETURNS SPHERICAL CAP POLE POSITION FLATO,FLONO
C             AND HALF-ANGLE THETA AS BN,BE, AND BV AFTER INITIALIZATION.
C             ON SUBSEQUENT CALLS:  ACTS AS L=1.
C          1  COMPUTES POTENTIAL FIELD COMPONENTS FROM INTERNAL COEFFICIENTS.
C          2  COMPUTES POTENTIAL FIELD COMPONENTS FROM EXTERNAL COEFFICIENTS.
C          3  COMPUTES FIELD FROM BOTH INTERNAL AND EXTERNAL COEFFICIENTS.
C         -1  COMPUTES GENERAL FUNCTION BV AND DERIVATIVES BN WITH RESPECT TO
C             LATITUDE AND BE WITH RESPECT TO LONGITUDE DIVIDED BY COS(LAT)
C             (R IS DUMMY VARIABLE IN THIS CASE).
C     NOTE:   SUBROUTINE IS INITIALIZED DURING FIRST CALL REGARDLESS OF L.

C     SUBPROGRAM USED:  LEGFUN

C	***** PARAMS & COEFFS TRANSFERRED FROM MAIN PROGRAM *****

C	ADAPTED FROM SUBROUTINE SCHNEV OF G.V. HAINES (COMPUTERS & GEOSCIENCES, 
C      14, 413-447, 1988)
C-------------------------------------------------------------------*/
	const size_t kdim = 6, ldim = 4;
	double constp[kdim + 1][kdim + 1],
		const_tmp[] = { 1., 1., 1., 1., 1., 1.73205, 0.866025, 1., 2.44949,
		1.93649, 0.7905691, 1., 3.16228, 3.35410, 2.09165, 0.739510, 1.,
		3.87298, 5.12348, 4.18330, 2.21853, 0.701561, 1., 4.58258, 7.24569,
		7.24569, 4.96078, 2.32681, 0.671693 };
	size_t i = 0;
	for (size_t n = 0; n <= kdim; n++)
		for (size_t m = 0; m <= n; m++) {
			constp[n][m] = const_tmp[i];
			i++;
		}
	/*C     IBF   =  0   TO USE ORDINARY POLYNOMIALS AS BASIS FUNCTIONS
	C              1          LEGENDRE POLYNOMIALS
	C              2          FOURIER SERIES
	C              3          COSINE SERIES
	C              4          SINE SERIES
	C     NOTE:    TZERO AND THINT MAY DEPEND ON IBF.*/
	int ibf = 2, kt, incept;
	double t1 = 1., t2 = 12., thint, aor, ar, aor3, delt[5], del, st, recip, ct;
	tbfit(t1, t2, ibf, thint, amtb.tzero);

	if (l >= 0) {
		if (amtb.ifit < 0) return;
		aor = amtb.re / r;
		ar = aor * aor;
		if (l > 1) {
			if (amtb.kext > 0) aor3 = aor * ar;
			if (l > 2) kt = amtb.lint >= amtb.lext ? amtb.lint : amtb.lext;
			else kt = amtb.lext;
		}
		else
			kt = amtb.lint;
	}
	else {
		if (amtb.ifit >= 0) return;
		ar = -1;
		kt = amtb.lint;
	}
	delt[0] = 1.;
	if (kt > 0) {
		del = (t - amtb.tzero) / thint;
		for (size_t i = 1; i <= kt; i++) {
			if (i == 1) {
				if (ibf <= 1)
					delt[i] = del;
				else if (ibf == 2) {
					st = sin(del);
					delt[i] = st;
				}
				else if (ibf == 3)
					delt[i] = cos(del);
				else
					delt[i] = sin(del);
				continue;
			}
			if (ibf == 0)
				delt[i] = delt[i - 1] * del;
			else if (ibf == 1) {
				recip = 1. / i;
				delt[i] = (2 - recip) * delt[i - 1] * del - (1 - recip)*delt[i - 2];
			}
			else if (ibf == 2)
				if (i % 2 == 0)
					if (i == 2) {
						ct = cos(del);
						delt[i] = ct;
					}
					else
						delt[i] = delt[i - 2] * ct - delt[i - 3] * st;
				else
					delt[i] = delt[i - 2] * ct - delt[i - 1] * st;
			else if (ibf == 3)
				delt[i] = cos(i * del);
			else if (ibf == 4)
				delt[i] = sin(i * del);
			else
				return;
		} //cycle 
		incept = 0;
		if (((ibf == 2) || (ibf == 3)) && incept == 1)
			for (size_t i = 2; i <= amtb.lint; i += 4 - ibf)
				delt[i] = 1 - delt[i];
	} // 103
	double x = 0, y = 0, ffn, p, dp, pms, gte, hti, hte, ra,
		z = 0, gti, colat, cl, sl, cml[8], sml[8], bt, bt1, bt2, bt3;
	int nmm;
	if ((l != 2) && (amtb.kint >= 0)) {
		gti = 0;
		for (size_t i = 0; i <= amtb.lint; i++)
			gti += amtb.bint[0][0][i] * delt[i];
		z = -ar * gti;
	}
	colat = 90 - flat;
	for (size_t n = 1; n <= amtb.kmax; n++) {
		if (n == 1) {
			cl = cos(flon * constants::umr);
			sl = sin(flon * constants::umr);
			cml[0] = cl;
			sml[0] = sl;
		}
		else {
			sml[n - 1] = sl * cml[n - 2] + cl * sml[n - 2];
			cml[n - 1] = cl * cml[n - 2] - sl * sml[n - 2];
		}
		for (size_t m = 0; m <= n; m++) {
			if (amtb.ib != 2) {
				nmm = n - m;
				if (nmm % 2) continue;
			}
			ffn = amtb.fn[n][m];
			legfun(m, ffn, constp[n][m], colat, p, dp, pms);
			if (l >= 0)
				ar = pow(aor, ffn + 2);
			else {
				ar = 1;
				ffn = -2;
				dp = -dp;
				pms = -pms;
			}
			if (!m) {
				bt1 = 0;
				bt3 = 0;
				bt = 0;
				if ((l != 2) && (n <= amtb.kint)) {
					gti = 0;
					for (size_t i = 0; i <= amtb.lint; i++)
						gti += amtb.bint[n][m][i];
					bt1 = ar * gti;
					bt3 = bt1;
				}
				if ((l < 1) && (n <= amtb.kext)) {
					gte = 0;
					for (size_t i = 0; i <= amtb.lext; i++)
						gte = gte + amtb.bext[n][m][i] * delt[i];
					bt = aor3 / ar * gte;
					bt1 = bt1 + bt;
				}
				x += bt1 * dp;
				z -= (ffn * (bt3 - bt) + bt3) * p;
				break;
			} // 130
			bt1 = 0.;
			bt2 = 0.;
			bt3 = 0.;
			bt = 0.;
			if ((l != 2) && (n <= amtb.kint)) {
				gti = 0;
				hti = 0;
				for (size_t i = 0; i <= amtb.lint; i++) {
					gti += amtb.bint[n][m][i] * delt[i];
					hti += amtb.bint[m - 1][n][i] * delt[i];
				}
				bt1 = ar * (gti * cml[m - 1] + hti * sml[m - 1]);
				bt2 = ar * (gti * sml[m] - hti * cml[m]);
				bt3 = bt1;
			} // 133
			if ((l > 1) && (n <= amtb.kext)) {
				gte = 0;
				hte = 0;
				for (size_t i = 0; i <= amtb.lext; i++) {
					gte += amtb.bext[n][m][i] * delt[i];
					hte += amtb.bext[m - 1][n][i] * delt[i];
				}
				ra = aor3 / ar;
				bt = ra * (gte * cml[m - 1] + hte * sml[m - 1]);
				bt1 = bt1 + bt;
				bt2 = bt2 + ra * (gte * sml[m - 1] - hte * cml[m - 1]);
			} // 135
			x += bt1 * dp;
			y += bt2 * pms;
			z -= (ffn * (bt3 - bt) + bt3) * p;
		}
	}
	bn = x;
	be = y;
	bv = z;
}

double fof1ed(double ylati, double r, double chi) {
/*c--------------------------------------------------------------
C CALCULATES THE F1 PEAK PLASMA FREQUENCY (FOF1/MHZ)
C INPUT:   
C		YLATI	ABSOLUT VALUE OF DIP-LATITUDE IN DEGREE
c       R		12-MONTH RUNNING MEAN OF SUNSPOT NUMBER 
c       CHI		SOLAR ZENITH ANGLE IN DEGREE
C REFERENCE: 
c       E.D.DUCHARME ET AL., RADIO SCIENCE 6, 369-378, 1971
C                                      AND 8, 837-839, 1973
c       HOWEVER WITH MAGNETIC DIP LATITUDE INSTEAD OF GEOMAGNETIC
c       DIPOLE LATITUDE, EYFRIG, 1979                    
C--------------------------------------------- D. BILITZA, 1988.*/
	if (chi > 90) return 0;
	double dla = ylati,
		f0 = 4.35 + dla * (0.0058 - 1.2e-4 * dla),
		f100 = 5.348 + dla * (0.011 - 2.3e-4 * dla),
		fs = f0 + (f100 - f0) * r / 100.0,
		xmue = 0.093 + dla * (0.0046 - 5.4e-5 * dla) + 3.0e-4 * r,
		fof1 = pow(fs * cos(chi * constants::umr), xmue),
		chi0 = 49.84733 + 0.349504 * dla,
		chi100 = 38.96113 + 0.509932 * dla,
		chim = (chi0 + (chi100 - chi0) * r / 100.);
	if (chi > chim) fof1 = -fof1;
	return fof1;
}

double f1_c1(double xmodip, double hour, double suxnon, double saxnon) {
/*c F1 layer shape parameter C1 after Reinisch and Huang, Advances in
c Space Research, Volume 25, Number 1, 81-88, 2000.*/
	double absmdp = abs(xmodip),
		dela = 4.32;
	if (absmdp >= 18) dela = 1 + exp(-(absmdp - 30.0) / 10.0);
	double c1old = 0.09 + 0.11 / dela,
		c1;
	if (suxnon == saxnon)
		c1 = 2.5 * c1old;
	else
		c1 = 2.5 * c1old * cos((hour - 12.) / (suxnon - saxnon) * constants::pi);
	return c1 < 0 ? 0 : c1;
}

void f1_prob(double sza, double glat, double rz12,
	double &f1prob, double &f1probl) {
/*c--------------------------------------------------------------------------
c Occurrence probability of F1 layer after Scotto et al., Advances in
c Space Research, Volume 20, Number 9, 1773-1775, 1997.
c
c Input: 	sza		solar zenith angle in degrees 
c 			glat	geomagnetic latitude in degrees
C			rz12	12-month running mean of sunspot number
c Output: 	f1prob	F1 occurrence probability without L-condition cases 
c 			f1probl	F1 occurrence probability with L-condition cases
c--------------------------------------------------------------------------*/
	double xarg = 0.5 + 0.5 * cos(sza*constants::umr),
		a = 2.98 + 0.0854 * rz12,
		b = 0.0107 - 0.0022 * rz12,
		c = -0.000256 + 0.0000147 * rz12,
		gamma = a + (b + c * glat) * glat,
		f1pr = pow(xarg, gamma);
	if (f1pr < 1e-3) f1pr = 0;
	f1prob = f1pr;
	double f1prl = pow(xarg, 2.36);
	if (f1prl < 1e-3) f1prl = 0.0;
	f1probl = f1prl;
}

double hpol(double hour, double tw, double xnw, double sa, double su,
	double dsa, double dsu) {
/*C-------------------------------------------------------
C PROCEDURE FOR SMOOTH TIME-INTERPOLATION USING EPSTEIN  
C STEP FUNCTION AT SUNRISE (SA) AND SUNSET (SU). THE 
C STEP-WIDTH FOR SUNRISE IS DSA AND FOR SUNSET DSU.
C TW,NW ARE THE DAY AND NIGHT VALUE OF THE PARAMETER TO 
C BE INTERPOLATED. SA AND SU ARE TIME OF SUNRIES AND 
C SUNSET IN DECIMAL HOURS.
C BILITZA----------------------------------------- 1979.*/
	double res;
	if (abs(su) > 25) {
		if (su > 0)
			res = tw;
		else
			res = xnw;
	}
	else
		res = xnw + (tw - xnw) * epst(hour, dsa, sa) + (xnw - tw) * epst(hour, dsu, su);
	return res;
}

double epst(double x, double sc, double hx) {
	double d1 = (x - hx) / sc;
	if (abs(d1) < constants::argmax)
		return 1. / (1. + exp(-d1));
	else
		return d1 > 0 ? 1 : 0;
}

void tal(double shabr, double sdelta, double shbr, double sdtdh0, bool aus6, double * spt) {
/*C-----------------------------------------------------------
C CALCULATES THE COEFFICIENTS SPT FOR THE POLYNOMIAL
C Y(X)=1+SPT(1)*X**2+SPT(2)*X**3+SPT(3)*X**4+SPT(4)*X**5               
C TO FIT THE VALLEY IN Y, REPRESENTED BY:                
C Y(X=0)=1, THE X VALUE OF THE DEEPEST VALLEY POINT (SHABR),                    
C THE PRECENTAGE DEPTH (SDELTA), THE WIDTH (SHBR) AND THE                       
C DERIVATIVE DY/DX AT THE UPPER VALLEY BOUNDRY (SDTDH0).                        
C IF THERE IS AN UNWANTED ADDITIONAL EXTREMUM IN THE VALLEY                     
C REGION, THEN AUS6=.TRUE., ELSE AUS6=.FALSE..     
C FOR -SDELTA THE COEFF. ARE CALCULATED FOR THE FUNCTION                        
C Y(X)=EXP(SPT(1)*X**2+...+SPT(4)*X**5).           
C-----------------------------------------------------------*/
	aus6 = false;
	if (shbr <= 0) {
		aus6 = true;
		return;
	}
	double z1 = -sdelta / (100.0 * shabr * shabr);
	if (sdelta <= 0) {
		sdelta = -sdelta;
		z1 = log(1. - sdelta / 100.) / (shabr * shabr);
	}
	double z3 = sdtdh0 / (2. * shbr),
		z4 = shabr - shbr;
	spt[3] = 2.0 * (z1 * (shbr - 2.0 * shabr) * shbr + z3 * z4 * shabr) /
		(shabr * shbr * z4 * z4 * z4);
	spt[2] = z1 * (2.0 * shbr - 3.0 * shabr) / (shabr * z4 * z4) -
		(2.*shabr + shbr) * spt[3];
	spt[1] = -2.0 * z1 / shabr - 2.0 * shabr * spt[2] - 3.0 * shabr * shabr * spt[3];
	spt[0] = z1 - shabr * (spt[1] + shabr * (spt[2] + shabr * spt[3]));
	double b = 4. * spt[2] / (5. * spt[3]) + shabr,
		c = -2. * spt[0] / (5 * spt[3] * shabr),
		z2 = b * b / 4. - c;
	if (z2 < 0) return;
	z3 = sqrt(z2);
	z1 = b / 2.;
	z2 = -z1 + z3;
	if ((z2 > 0) && (z2 < shbr)) aus6 = true;
	if (abs(z3) > 1e-15) {
		z2 = -z1 - z3;
		if ((z2 > 0) && (z2 < shbr)) aus6 = true;
		return;
	}
	z2 = c / z2;
	if ((z2 > 0) && (z2 < shbr)) aus6 = true;
}

double xmded(double xhi, double r, double yw) {
/*C D. BILITZA, 1978, CALCULATES ELECTRON DENSITY OF D MAXIMUM.                   
C XHI/DEG. IS SOLAR ZENITH ANGLE, R SMOOTHED ZURICH SUNSPOT NUMBER              
C AND YW/M-3 THE ASSUMED CONSTANT NIGHT VALUE.     
C [REF.: D.BILITZA, WORLD DATA CENTER A REPORT UAG-82,7,BOULDER,1981]
C corrected 4/25/97 - D. Bilitza*/
	if (xhi >= 90) return yw;
	double y = 6.05e8 + 0.088e8 * r,
		yy = cos(xhi * constants::umr),
		yyy = -0.1 / pow(yy, 2.7), ymd;
	if (yyy < -40)
		ymd = 0;
	else
		ymd = y * exp(yyy);
	if (ymd < yw) ymd = yw;
	return ymd;
}

void apfmsis(int isdate, double hour, double *iapo) {
/*c-----------------------------------------------------------------------
c Finds 3-hourly Ap indices for NRLMSIS00 model 
c    INPUTS: 	ISDATE		Array-index from APF_ONLY
c				HOUR		UT in decimal hours
c    OUTPUT:   	IAPO(1:7)	3-hourly Ap index
C  
C IAPO(1) DAILY AP
C IAPO(2) 3-HR AP INDEX FOR CURRENT TIME			   	
C IAPO(3) 3-HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME	
C IAPO(4) 3-HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME	
C IAPO(5) 3-HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME  
C IAPO(6) AVERAGE OF EIGHT 3-HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C         TO CURRENT TIME
C IAPO(7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
C         TO CURRENT TIME
c
c The 3-hour UT intervals during the day are: (0-3),)3-6),)6-9),)9-12),
c )12-15),)15-18),)18-21),)21-24(.
c 
c If date is outside the range of the Ap indices file then IAPO(2)=-5  
c-----------------------------------------------------------------------*/
	int is = isdate,
		ihour = int(hour / 3.) + 1;
	if (ihour > 8) ihour = 8;

	iapo[0] = apfa.aap[is - 1][8];
	// There must be at least 20 indices available
	if ((is - 1) * 8 + ihour < 20) {
		iapo[1] = -5.;
		return;
	}
	// assemble Ap values as needed by MSIS
	int j1 = ihour + 1, iap[20];
	for (size_t i = 1; i <= ihour; i++)
		iap[j1 - i - 1] = apfa.aap[is - 1][i - 1];

	j1 = ihour + 9;
	for (size_t i = 1; i <= 8; i++)
		iap[j1 - i - 1] = apfa.aap[is - 2][i - 1];

	j1 = ihour + 17;
	int j2 = 8 - (20 - ihour - 8) + 1;
	if (j2 < 1) j2 = 1;
	for (size_t i = j2; i <= 8; i++)
		iap[j1 - i - 1] = apfa.aap[is - 3][i - 1];

	if (ihour < 4) {
		j1 = ihour + 25;
		j2 = 8 - (20 - ihour - 16) + 1;
		for (size_t i = 1; i <= 4; i++)
			iap[j1 - i - 1] = apfa.aap[is - 4][i - 1];
	}

	for (size_t i = 1; i <= 4; i++)
		iapo[i] = iap[i - 1];
	
	double sum1 = 0,
		sum2 = 0;
	for (size_t i = 1; i <= 8; i++)	{
		sum1 += iap[3 + i];
		sum2 += iap[11 + i];
	}
	iapo[5] = sum1 / 8;
	iapo[6] = sum2 / 8;
}

void tselec(double *sv) {
/*C-----------------------------------------------------------------------
C        SET SWITCHES
C        Output in  COMMON/CSW/SW(25),ISW,SWC(25)
C        SW FOR MAIN TERMS, SWC FOR CROSS TERMS
C  
C        TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SV),
C        WHERE SV IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
C        FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C
C        To get current values of SW: CALL TRETRV(SW)
C-----------------------------------------------------------------------
*/
	double *sav = new double[25];

	for (size_t i = 1; i <= 25; i++) {
		sav[i - 1] = sv[i - 1];
		csw.sw[i - 1] = remainder(sv[i - 1], 2);
		if ((abs(sv[i - 1]) == 1) || (abs(sv[i - 1]) == 2))
			csw.swc[i - 1] = 1;
		else
			csw.swc[i - 1] = 0;
	}
	csw.isw = 64999;
}

void gtd7(int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
	double f107, double *ap, int mass, double *d, double *t) {
/*C-----------------------------------------------------------------------
C
C     NRLMSISE-00
C     -----------
C        Neutral Atmosphere Empirical Model from the surface to lower
C        exosphere
C        J.M. Picone, A.E. Hedin, D.P. Drob, and A.C. Aikin, NRLMSISE-00 
C             empirical model of the atmosphere: Statistical comparisons 
C             and scientific issues, J. Geophys. Res., 107(A12), 1468, 
C             doi:10.1029/2002JA009430, 2002.
C
C        NEW FEATURES:
C          *Extensive satellite drag database used in model generation
C          *Revised O2 (and O) in lower thermosphere
C          *Additional nonlinear solar activity term
C          *"ANOMALOUS OXYGEN" NUMBER DENSITY, OUTPUT D(9)
C           At high altitudes (> 500 km), hot atomic oxygen or ionized
C           oxygen can become appreciable for some ranges of subroutine
C           inputs, thereby affecting drag on satellites and debris. We
C           group these species under the term "anomalous oxygen," since
C           their individual variations are not presently separable with
C           the drag data used to define this model component.
C
C        SUBROUTINES FOR SPECIAL OUTPUTS:
C        
C        HIGH ALTITUDE DRAG: EFFECTIVE TOTAL MASS DENSITY 
C        (SUBROUTINE GTD7D, OUTPUT D(6))
C           For atmospheric drag calculations at altitudes above 500 km,
C           call SUBROUTINE GTD7D to compute the "effective total mass
C           density" by including contributions from "anomalous oxygen."
C           See "NOTES ON OUTPUT VARIABLES" below on D(6).
C
C        PRESSURE GRID (SUBROUTINE GHP7)
C          See subroutine GHP7 to specify outputs at a pressure level
C          rather than at an altitude.
C
C        OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
C 
C     INPUT VARIABLES:
C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
C              (Year ignored in current model)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C                 MASS 17 IS Anomalous O ONLY.)
C
C     NOTES ON INPUT VARIABLES: 
C        UT, Local Time, and Longitude are used independently in the
C        model and are not of equal importance for every situation.  
C        For the most physically realistic calculation these three
C        variables should be consistent (STL=SEC/3600+GLONG/15).
C        The Equation of Time departures from the above formula
C        for apparent local time can be included if available but
C        are of minor importance.
c
C        F107 and F107A values used to generate the model correspond
C        to the 10.7 cm radio flux at the actual distance of the Earth
C        from the Sun rather than the radio flux at 1 AU. The following
C        site provides both classes of values:
C        ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
C
C        F107, F107A, and AP effects are neither large nor well
C        established below 80 km and these parameters should be set to
C        150., 150., and 4. respectively.
C
C     OUTPUT VARIABLES:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)                       
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
C     NOTES ON OUTPUT VARIABLES:
C        TO GET OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.) 
C
C        O, H, and N are set to zero below 72.5 km
C
C        T(1), Exospheric temperature, is set to global average for
C        altitudes below 120 km. The 120 km gradient is left at global
C        average value for altitudes below 72 km.
C
C        D(6), TOTAL MASS DENSITY, is NOT the same for subroutines GTD7 
C        and GTD7D
C
C          SUBROUTINE GTD7 -- D(6) is the sum of the mass densities of the
C          species labeled by indices 1-5 and 7-8 in output variable D.
C          This includes He, O, N2, O2, Ar, H, and N but does NOT include
C          anomalous oxygen (species index 9).
C
C          SUBROUTINE GTD7D -- D(6) is the "effective total mass density
C          for drag" and is the sum of the mass densities of all species
C          in this model, INCLUDING anomalous oxygen.
C        
C     SWITCHES: The following is for test and special purposes:
C          
C        TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW),
C        WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
C        FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C        FOR THE FOLLOWING VARIATIONS
C               1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
C               3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
C               5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
C               7 - DIURNAL               8 - SEMIDIURNAL
C               9 - DAILY AP             10 - ALL UT/LONG EFFECTS
C              11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
C              13 - MIXED AP/UT/LONG     14 - TERDIURNAL
C              15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
C              16 - ALL TINF VAR         17 - ALL TLB VAR
C              18 - ALL TN1 VAR           19 - ALL S VAR
C              20 - ALL TN2 VAR           21 - ALL NLB VAR
C              22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR
C        NOTES: The default should be sw(1:25)=1.0; Uses the full Ap history 
C        if sw(9) = -1.0; sw(9)=-1 uses the full ap history and is supposed 
C		 to give the best results for geomagentic storms; sw(9) = 0 turns off 
C        all ap variations.
C
C        To get current values of SW: CALL TRETRV(SW)
C
C Simplified version for F77 compilers and IRI (Sep 27,2010, dbilitza):
C		GTS7, GLOBE7: AP(1), P(1) -> AP(7), P(150)
C   		PD(150,9) -> PDA1(150),..,PDA9(150) 
C-----------------------------------------------------------------------*/
	double zn3[] = { 32.5, 20., 15., 10., 0. },
		zn2[] = { 72.5, 55., 45., 32.5 },
		zmix = 62.5,
		alast = 99999,
		mssl = -999;
	int mn3 = 5, mn2 = 4;
	// Put identification data into common/datime/
	for (size_t i = 1; i <= 3; i++)
		datime.isdate[i - 1] = datim7.isd[i - 1];
	for (size_t i = 1; i <= 2; i++) {
		datime.istime[i - 1] = datim7.ist[i - 1];
		datime.name[i - 1] = datim7.nam[i - 1];
	}
	// Test for changed input
	double v1 = vtst7(iyd, sec, glat, glong, stl, f107a, f107, ap, 1);
	// Latitude variation of gravity (none for SW(2)=0)
	double xlat = glat;
	if (csw.sw[1] == 0) xlat = 45;
	double gsurf, re, xmm = lower7.pdm[4][2];
	glatf(xlat, gsurf, re);
	// THERMOSPHERE/MESOSPHERE (above ZN2(1))
	double altt = alt ? alt >= zn2[0] : zn2[0],
		mss = mass;
	// Only calculate N2 in thermosphere if alt in mixed region
	if ((alt < zmix) && (mass > 0)) mss = 28;
	// Only calculate thermosphere if input parameters changed
	// or altitude above ZN2(1) in mesosphere
	double *ds = new double[9],
		*ts = new double[2], dm28m;
	if ((v1 == 1) || (alt > zn2[0]) || (alast > zn2[0]) || (mss != mssl)) {
		ds[0] = d[0]; // switches on pgr mod
		gts7(iyd, sec, altt, glat, glong, stl, f107a, f107, ap, mss, ds, ts);
		dm28m = dmix.dm28;
		// metric adjustment
		if (metsel.imr == 1) dm28m = dmix.dm28 * 1.e6;
		mssl = mss;
	}
}

double vtst7(int iyd, double sec, double glat, double glong, double stl, double f107a,
	double f107, double *ap, int ic) {
	double iydl[] = { -999, -999 },
		secl[] = { -999, -999 },
		glatl[] = { -999, -999 },
		gll[] = { -999, -999 },
		stll[] = { -999, -999 },
		fal[] = { -999, -999 },
		fl[] = { -999, -999 },
		apl[7][2], swl[25][2], swcl[25][2],
		vtst7 = 0;
	for (size_t i = 0; i < 7; i++)
		for (size_t j = 0; j < 2; j++)
			apl[i][j] = -999;
	for (size_t i = 0; i < 25; i++)
		for (size_t j = 0; j < 2; j++) {
			swl[i][j] = -999;
			swcl[i][j] = -999;
		}
	bool flag_t = false;
	if ((iyd != iydl[ic - 1]) || (sec != secl[ic - 1]) ||
		(glat != glatl[ic - 1]) || (glong != gll[ic - 1]) ||
		(stl != stll[ic - 1]) || (f107a != fal[ic - 1]) ||
		(f107 != fl[ic - 1]))
		flag_t = true;

	for (size_t i = 1; i <= 7; i++)
		if (ap[i - 1] != apl[i - 1][ic - 1]) {
			flag_t = true;
			break;
		}
	for (size_t i = 1; i <= 25; i++) {
		if ((csw.sw[i - 1] != swl[i - 1][ic - 1]) ||
			(csw.swc[i - 1] != swcl[i - 1][ic - 1])) {
			flag_t = true;
			break;
		}
	}
	if (flag_t)	{
		vtst7 = 1;
		iydl[ic - 1] = iyd;
		secl[ic - 1] = sec;
		glatl[ic - 1] = glat;
		gll[ic - 1] = glong;
		stll[ic - 1] = stl;
		fal[ic - 1] = f107a;
		fl[ic - 1] = f107;
		for (size_t i = 1; i <= 7; i++)
		apl[i - 1][ic - 1] = ap[i];
		for (size_t i = 1; i <= 25; i++) {
			swl[i - 1][ic - 1] = csw.sw[i - 1];
			swcl[i - 1][ic - 1] = csw.swc[i - 1];
		}
	}
}

void glatf(double lat, double &gv, double &reff) {
	/*C-----------------------------------------------------------------------
	C      CALCULATE LATITUDE VARIABLE GRAVITY (GV) AND EFFECTIVE
	C      RADIUS (REFF)
	C-----------------------------------------------------------------------*/
	double dgtr = 1.74533e-2,
		c2 = cos(2. * dgtr * lat);
	gv = 980.616 * (1. - .0026373 * c2);
	reff = 2. * gv / (3.085462e-6 + 2.27e-9 * c2) * 1e-5;
}

void gts7(int iyd, double sec, double alt, double glat, double glong, double stl,
	double f107a, double f107, double *ap, int mass, double *d, double *t) {
	/*C-----------------------------------------------------------------------
	C
	C     Thermospheric portion of NRLMSISE-00
	C     See GTD7 for more extensive comments
	C
	C        OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
	C
	C     INPUT VARIABLES:
	C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
	C              (Year ignored in current model)
	C        SEC - UT(SEC)
	C        ALT - ALTITUDE(KM) (>72.5 km)
	C        GLAT - GEODETIC LATITUDE(DEG)
	C        GLONG - GEODETIC LONGITUDE(DEG)
	C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
	C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
	C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
	C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
	C           - ARRAY CONTAINING:
	C             (1) DAILY AP
	C             (2) 3 HR AP INDEX FOR CURRENT TIME
	C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
	C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
	C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
	C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
	C                    TO CURRENT TIME
	C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
	C                    TO CURRENT TIME
	C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
	C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
	C                 MASS 17 IS Anomalous O ONLY.)
	C
	C     NOTES ON INPUT VARIABLES:
	C        UT, Local Time, and Longitude are used independently in the
	C        model and are not of equal importance for every situation.
	C        For the most physically realistic calculation these three
	C        variables should be consistent (STL=SEC/3600+GLONG/15).
	C        The Equation of Time departures from the above formula
	C        for apparent local time can be included if available but
	C        are of minor importance.
	c
	C        F107 and F107A values used to generate the model correspond
	C        to the 10.7 cm radio flux at the actual distance of the Earth
	C        from the Sun rather than the radio flux at 1 AU. The following
	C        site provides both classes of values:
	C        ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
	C
	C        F107, F107A, and AP effects are neither large nor well
	C        established below 80 km and these parameters should be set to
	C        150., 150., and 4. respectively.
	C
	C     OUTPUT VARIABLES:
	C        D(1) - HE NUMBER DENSITY(CM-3)
	C        D(2) - O NUMBER DENSITY(CM-3)
	C        D(3) - N2 NUMBER DENSITY(CM-3)
	C        D(4) - O2 NUMBER DENSITY(CM-3)
	C        D(5) - AR NUMBER DENSITY(CM-3)
	C        D(6) - TOTAL MASS DENSITY(GM/CM3) [Anomalous O NOT included]
	C        D(7) - H NUMBER DENSITY(CM-3)
	C        D(8) - N NUMBER DENSITY(CM-3)
	C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
	C        T(1) - EXOSPHERIC TEMPERATURE
	C        T(2) - TEMPERATURE AT ALT
	C-----------------------------------------------------------------------*/
	int tnmod = 0; // for switching on mod MSIS
	if (d[0] < 0) tnmod = -d[0];
	double zn1[] = { 120., 110., 100.,90., 72.5 };
	// Test for changed input
	double v2 = vtst7(iyd, sec, glat, glong, stl, f107a, f107, ap, 2),
		yrd = iyd;
	gts3c.za = parm7.pdl[15][1];
	zn1[0] = gts3c.za;
	for (size_t i = 0; i < 9; i++) // ???
		d[i] = 0;
	double alast = -999;
	// VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1)
	if (alt > zn1[0]) {
		if ((v2 == 1) || (alast <= zn1[0]))
			ttest.tinf = lower7.ptm[0] * parm7.pt[0] *
			(1. + csw.sw[15] * globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, parm7.pt));
	}
	else
		ttest.tinf = lower7.ptm[0] * parm7.pt[0];
	if (tnmod > 600) ttest.tinf = tnmod;
	t[0] = ttest.tinf;
	// GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5)
	if (alt > zn1[4]) {
		if ((v2 == 1) || (alast <= zn1[4]))
			gts3c.g0 = lower7.ptm[3] * parm7.ps[0] *
			(1 + csw.sw[18] * globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, parm7.ps));
	}
	else
		gts3c.g0 = lower7.ptm[3] * parm7.ps[0];
	// Calculate these temperatures only if input changed
	if ((v2 == 1) || (alt < 300))
		gts3c.tlb = lower7.ptm[1] * (1 + csw.sw[16] *
			globe7(yrd, sec, glat, glong, stl, f107a, f107, ap, parm7.pda4)) * parm7.pda4[0];
	gts3c.s = gts3c.g0 / (ttest.tinf - gts3c.tlb);
	// lower thermosphere temp variations not significant for
	// density above 300 km
	if(alt < 300) {
		if ((v2 == 1) || (alast >= 300)) {
			meso7.tn1[1] = lower7.ptm[6] * parm7.ptl[0][0] / (1 - csw.sw[17] * glob7s(&parm7.ptl[0][0]));
			meso7.tn1[2] = ptm(3)*ptl(1, 2) / (1. - sw(18)*glob7s(ptl(1, 2)))
				tn1(4) = ptm(8)*ptl(1, 3) / (1. - sw(18)*glob7s(ptl(1, 3)))
				tn1(5) = ptm(5)*ptl(1, 4) / (1. - sw(18)*sw(20)*glob7s(ptl(1, 4)))
				tgn1(2) = ptm(9)*pma(1, 9)*(1. + sw(18)*sw(20)*glob7s(pma(1, 9)))
				$   *tn1(5)*tn1(5) / (ptm(5)*ptl(1, 4))**2
		}
		endif
		else
			tn1(2) = ptm(7)*ptl(1, 1)
			tn1(3) = ptm(3)*ptl(1, 2)
			tn1(4) = ptm(8)*ptl(1, 3)
			tn1(5) = ptm(5)*ptl(1, 4)
			tgn1(2) = ptm(9)*pma(1, 9)
			$  *tn1(5)*tn1(5) / (ptm(5)*ptl(1, 4))**2
			endif
			c
			z0 = zn1(4)
			t0 = tn1(4)
			tr12 = 1.;
}

double globe7(double yrd, double sec, double xlat, double xlong,
	double tloc, double f107a, double f107, double *ap, double p[]) {
/*C-----------------------------------------------------------------------
C       CALCULATE G(L) FUNCTION 
C       Upper Thermosphere Parameters
C-----------------------------------------------------------------------*/
	// 3hr Magnetic activity functions
	// Eq.A24d
	std::function<double(double)> g0 = [p](double a) {
		return (a - 4 + (p[25] - 1) *
			(a - 4 + (exp(-abs(p[24]) * (a - 4)) - 1.) /
				abs(p[24]))); 
	};
	// Eq. A24c
	std::function<double(double)> sumex = [](double ex) {
		return 1 + (1 - pow(ex, 19)) / (1 - ex) * pow(ex, 5);
	};
	// Eq. A24a
	std::function<double(double)> sg0 = [ap, g0, sumex](double ex) {
		return (g0(ap[1]) + (g0(ap[2]) * ex + g0(ap[3]) * pow(ex, 2) + 
			g0(ap[4]) * pow(ex, 3) + (g0(ap[5]) * pow(ex, 4) +
			g0(ap[6]) * pow(ex, 12)) * (1 - pow(ex, 8)) / (1. - ex))) / sumex(ex);
	};
	for (size_t i = 0; i < 15; i++)
		ttest.t[i] = 0;
	double sw8;
	if (csw.sw[8] > 0) sw8 = 1;
	if (csw.sw[8] < 0) sw8 = -1;
	double iyr = yrd / 1000,
		day = yrd - iyr * 1000, xl=1000,
		dgtr = 1.74533e-2, tll=1000,
		hr=0.2618;
	// Eq. A22 (remainder of code)
	if (xl != xlat) {
		double c = sin(xlat * dgtr),
			s = cos(xlat * dgtr),
			c2 = c * c,
			c4 = c2 * c2,
			s2 = s * s;
		lpoly.plg[1][0] = c;
		lpoly.plg[2][0] = 0.5 * (3. * c2 - 1.);
		lpoly.plg[3][0] = 0.5 * (5. * c * c2 - 3. * c);
		lpoly.plg[4][0] = (35. * c4 - 30. * c2 + 3.) / 8.;
		lpoly.plg[5][0] = (63. * c2 * c2 * c - 70. * c2 * c + 15 * c) / 8.;
		lpoly.plg[6][0] = (11. * c * lpoly.plg[5][0] - 5. * lpoly.plg[4][0]) / 6.;
		lpoly.plg[1][1] = s;
		lpoly.plg[2][1] = 3 * c * s;
		lpoly.plg[3][1] = 1.5 * (5 * c2 - 1) * s;
		lpoly.plg[4][1] = 2.5 * (7. * c2 * c - 3. * c) * s;
		lpoly.plg[5][1] = 1.875 * (21. * c4 - 14. * c2 + 1.) * s;
		lpoly.plg[6][1] = (11. * c * lpoly.plg[5][1] - 6. * lpoly.plg[4][1]) / 5.;
		lpoly.plg[2][2] = 3 * s2;
		lpoly.plg[3][2] = 15 * s2 * c;
		lpoly.plg[4][2] = 7.5 * (7 * c2 - 1.) * s2;
		lpoly.plg[5][2] = 3 * c * lpoly.plg[4][2] - 2 * lpoly.plg[3][2];
		lpoly.plg[6][2] = (11 * c * lpoly.plg[5][2] - 7 * lpoly.plg[4][2]) / 4;
		lpoly.plg[7][2] = (13 * c * lpoly.plg[6][2] - 8 * lpoly.plg[5][2]) / 5;
		lpoly.plg[3][3] = 15 * s2 * s;
		lpoly.plg[4][3] = 105 * s2 * s * c;
		lpoly.plg[5][3] = (9 * c * lpoly.plg[4][3] - 7 * lpoly.plg[3][3]) / 2;
		lpoly.plg[6][3] = (11 * c * lpoly.plg[5][3] - 8 * lpoly.plg[4][3]) / 3;
		xl = xlat;
	}
	if ((tll != tloc) && ((csw.sw[6] != 0) || (csw.sw[7] != 0) || (csw.sw[13] != 0))) {
		lpoly.stloc = sin(hr * tloc);
		lpoly.ctloc = cos(hr * tloc);
		lpoly.s2tloc = sin(2 * hr * tloc);
		lpoly.c2tloc = cos(2 * hr * tloc);
		lpoly.s3tloc = sin(3 * hr * tloc);
		lpoly.c3tloc = cos(3 * hr * tloc);
		tll = tloc;
	}
	double dayl = -1, dr = 1.72142e-2, p13 = -1000, cd13,
		p17 = -1000, cd17, p31 = -1000, cd31, p38=-1000, cd38,
		f1, f2;
	if ((day != dayl) || (p[13] != p13)) cd13 = cos(dr * (day - p[13]));
	if ((day != dayl) || (p[17] != p17)) cd17 = cos(2 * dr * (day - p[17]));
	if ((day != dayl) || (p[31] != p31)) cd31 = cos(dr * (day - p[31]));
	if ((day != dayl) || (p[38] != p38)) cd38 = cos(2 * dr * (day - p[38]));
	dayl = day;
	p13 = p[13];
	p17 = p[17];
	p31 = p[31];
	p38 = p[38];
	// F10.7 EFFECT
	lpoly.df = f107 - f107a;
	lpoly.dfa = f107a - 150;
	ttest.t[0] = p[19] * lpoly.df * (1 + p[59] * lpoly.dfa) + 
		p[20] * lpoly.df * lpoly.df + p[21] * lpoly.dfa +
		p[29] * lpoly.dfa * lpoly.dfa;
	f1 = 1 + (p[47] * lpoly.dfa + p[19] * lpoly.df + p[20] * lpoly.df * lpoly.df) * csw.swc[0];
	f2 = 1. + (p[49] * lpoly.dfa + p[19] * lpoly.df + p[20] * lpoly.df * lpoly.df) * csw.swc[0];
	// TIME INDEPENDENT
	ttest.t[1] = (p[1] * lpoly.plg[2][0] + p[2] * lpoly.plg[4][0] + p[22] * lpoly.plg[6][0]) +
		(p[14] * lpoly.plg[2][0]) * lpoly.dfa * csw.sw[0];
	// SYMMETRICAL ANNUAL
	ttest.t[2] = p[18] * cd31;
	// SYMMETRICAL SEMIANNUAL
	ttest.t[3] = (p[15] + p[16] * lpoly.plg[2][0]) * cd17;
	// ASYMMETRICAL ANNUAL
	ttest.t[4] = f1 * (p[9] * lpoly.plg[1][0] + p[10] * lpoly.plg[3][0]) * cd13;
	// ASYMMETRICAL SEMIANNUAL
	ttest.t[5] = p[37] * lpoly.plg[1][0] * cd38;
	// DIURNAL
	double t71, t72;
	if (csw.sw[6] != 0) {
		t71 = (p[11] * lpoly.plg[2][1]) * cd13 * csw.swc[4];
		t72 = (p[12] * lpoly.plg[2][1]) * cd13 * csw.swc[4];
		ttest.t[6] = f2 * ((p[3] * lpoly.plg[1][1] + p[4] * lpoly.plg[3][1] +
			p[27] * lpoly.plg[5][1] + t71) * lpoly.ctloc +
			(p[6] * lpoly.plg[1][1] + p[7] * lpoly.plg[3][1] + 
			p[28] * lpoly.plg[5][1] + t72)	* lpoly.stloc);
	}
	// SEMIDIURNAL
	double t81, t82;
	if (csw.sw[7] != 0) {
		t81 = (p[23] * lpoly.plg[3][2] + p[35] * lpoly.plg[5][2]) * cd13 * csw.swc[4];
		t82 = (p[33] * lpoly.plg[3][2] + p[36] * lpoly.plg[5][2]) * cd13 * csw.swc[4];
		ttest.t[7] = f2 * ((p[5] * lpoly.plg[2][2] + p[41] * lpoly.plg[4][2] + t81) 
			* lpoly.c2tloc + 3 + (p[8] * lpoly.plg[2][2] + p[42] * lpoly.plg[4][2] + t82) *
			lpoly.s2tloc);
	}
	// TERDIURNAL
	if (csw.sw[13] != 0)
		ttest.t[13] = f2 *
			((p[39] * lpoly.plg[3][3] + (p[93] * lpoly.plg[4][3] + p[46] * lpoly.plg[6][3]) *
			cd13 * csw.swc[4]) * lpoly.s3tloc +
			(p[40] * lpoly.plg[3][3] + (p[94] * lpoly.plg[4][3] + p[48] * lpoly.plg[6][3]) *
			cd13 * csw.swc[4]) * lpoly.c3tloc);
	// MAGNETIC ACTIVITY BASED ON DAILY AP
	double p43, p44, exp1;
	if (sw8 != -1) {
		lpoly.apd = ap[0] - 4;
		p43 = p[43];
		p44 = p[44];
		if (p43 < 0) p43 = 1e-5;
		lpoly.apdf = lpoly.apd + (p44 - 1.) * (lpoly.apd + (exp(-p44 * lpoly.apd) - 1.) / p44);
		if (csw.sw[8] != 0)
			ttest.t[8] = lpoly.apdf * (p[32] + p[45] * lpoly.plg[2][0] + p[34] * lpoly.plg[4][0] +
				(p[100] * lpoly.plg[1][0] + p[101] * lpoly.plg[3][0] + p[102] * lpoly.plg[5][0]) *
				cd13 * csw.swc[4] + (p[121] * lpoly.plg[1][1] + p[122] * lpoly.plg[3][1] +
					p[123] * lpoly.plg[5][1]) * csw.swc[6] * cos(hr*(tloc - p[124])));
	}
	else {
		if (p[51] != 0) {
			exp1 = exp(-10800 * abs(p[51]) / (1 + p[138] * (45. - abs(xlat))));
			if (exp1 > 0.99999) exp1 = 0.99999;
			if (p[24] < 1e-4) p[24] = 1e-4;
			lpoly.apt[0] = sg0(exp1);
			if (csw.sw[8] != 0)
				ttest.t[8] = lpoly.apt[0] * (p[50] + p[96] * lpoly.plg[2][0] + 
					p[54] * lpoly.plg[4][0] + (p[125] * lpoly.plg[2][0] + p[126] *
						lpoly.plg[3][0] + p[127] * lpoly.plg[5][0]) * cd13 * csw.swc[4] +
					(p[128] * lpoly.plg[1][1] + p[129] * lpoly.plg[3][1] + p[130] * lpoly.plg[5][1]) *
					csw.swc[6] * cos(hr * (tloc - p[131])));
		}
	}
	double sr = 7.2722e-5;
	if ((csw.sw[9] != 0) && (xlong > -1000)) {
		// LONGITUDINAL
		if (csw.sw[10] != 0) {
			ttest.t[10] = (1 + p[80] * lpoly.dfa * csw.swc[0]) *
				((p[64] * lpoly.plg[2][1] +
					p[65] * lpoly.plg[4][1] +
					p[66] * lpoly.plg[6][1] +
					p[103] * lpoly.plg[1][1] +
					p[104] * lpoly.plg[3][1] +
					p[105] * lpoly.plg[5][1] +
					csw.swc[4] * (p[109] * lpoly.plg[1][1] +
						p[110] * lpoly.plg[3][1] +
						p[111] * lpoly.plg[5][1]) * cd13) * cos(dgtr * xlong) +
						(p[90] * lpoly.plg[2][1] +
							p[91] * lpoly.plg[4][1] +
							p[92] * lpoly.plg[6][1] +
							p[106] * lpoly.plg[1][1] +
							p[107] * lpoly.plg[3][1] +
							p[108] * lpoly.plg[5][1] +
							csw.swc[4] * (p[112] * lpoly.plg[1][1] +
								p[113] * lpoly.plg[3][1] +
								p[114] * lpoly.plg[5][1]) * cd13) * sin(dgtr * xlong));
		} // 230
		// UT AND MIXED UT, LONGITUDE
		if (csw.swc[1] != 0) {
			ttest.t[11] = (1 + p[95] * lpoly.plg[1][0]) *
				(1 + p[81] * lpoly.dfa * csw.swc[0]) *
				(1 + p[119] * lpoly.plg[1][0] * csw.swc[4] * cd13) *
				((p[68] * lpoly.plg[1][0] +
					p[69] * lpoly.plg[3][0] +
					p[70] * lpoly.plg[5][0]) * cos(sr * (sec - p[71])));
			ttest.t[11] += csw.swc[10] *
				(p[76] * lpoly.plg[3][2] +
					p[77] * lpoly.plg[5][2] +
					p[78] * lpoly.plg[7][2]) * cos(sr * (sec - p[79]) +
						2 * dgtr * xlong) * (1 + p[137] * lpoly.dfa * csw.swc[0]);
		} // 240
		// UT,LONGITUDE MAGNETIC ACTIVITY
		if (csw.swc[12] != 0) {
			if (sw8 != -1) {
				ttest.t[12] = lpoly.apdf * csw.swc[10] * (1 + p[120] * lpoly.plg[1][0]) *
					((p[60] * lpoly.plg[2][1] +
						p[61] * lpoly.plg[5][1] +
						p[62] * lpoly.plg[6][1]) * cos(dgtr * (xlong - p[63]))) +
					lpoly.apdf * csw.swc[10] * csw.swc[4] *
					(p[115] * lpoly.plg[1][1] +
						p[116] * lpoly.plg[3][1] +
						p[117] * lpoly.plg[5][1]) * cd13 * cos(dgtr * (xlong - p[118])) +
					lpoly.apdf * csw.swc[11] *
					(p[83] * lpoly.plg[1][0] +
						p[84] * lpoly.plg[3][0] +
						p[85] * lpoly.plg[5][0]) * cos(sr * (sec - p[75]));
			} 
			else { // 45
				if (p[51] != 0) {
					ttest.t[12] = lpoly.apt[0] * csw.swc[10] * (1 +
						p[132] * lpoly.plg[1][0]) *
						((p[52] * lpoly.plg[2][1] +
							p[98] * lpoly.plg[4][1] +
							p[67] * lpoly.plg[6][1]) * cos(dgtr*(xlong - p[97]))) +
						lpoly.apt[0] * csw.swc[10] * csw.swc[4] *
						(p[133] * lpoly.plg[1][1] +
							p[134] * lpoly.plg[3][1] +
							p[135] * lpoly.plg[5][1]) * cd13 * cos(dgtr * (xlong - p[136])) +
						lpoly.apt[0] * csw.swc[1] *
						(p[55] * lpoly.plg[1][0] +
							p[56] * lpoly.plg[3][0] +
							p[57] * lpoly.plg[5][0]) * cos(sr * (sec - p[58]));
				}
			}
		} // 48
		//  PARMS NOT USED: 83, 90,100,140-150
	} // 49
	ttest.tinf = p[30];
	for (size_t i = 1; i <= 14; i++) {
		ttest.tinf += abs(csw.sw[i - 1]) * ttest.t[i - 1];
	}
	return ttest.tinf;
}

double glob7s(double *p) {
/*C-----------------------------------------------------------------------
C      VERSION OF GLOBE FOR LOWER ATMOSPHERE 10/26/99
C-----------------------------------------------------------------------*/
	double pset = 2, dayl = -1,
		p31 = -1000, p17 = -1000, p13 = -1000, p38 = -1000,
		cd31, cd17, cd13, cd38, dr = 1.72142e-2;
	if (p[99] == 0) p[99] = pset;
	if (p[99] != pset)
		return;
	double *t = new double[14];
	for (size_t i = 0; i < 14; i++)
		t[i] = 0;
	if ((lpoly.day != dayl) || (p31 != p[31])) cd31 = cos(dr * (lpoly.day - p[31]));
	if ((lpoly.day != dayl) || (p17 != p[17])) cd17 = cos(2 * dr * (lpoly.day - p[17]));
	if ((lpoly.day != dayl) || (p13 != p[13])) cd13 = cos(dr * (lpoly.day - p[13]));
	if ((lpoly.day != dayl) || (p38 != p[38])) cd17 = cos(2 * dr * (lpoly.day - p[39]));
	dayl = lpoly.day;
	p31 = p[31];
	p17 = p[17];
	p13 = p[13];
	p38 = p[38];
	// F10.7
	t[0] = p[21] * lpoly.dfa;
	// time independent
	t[1] = p[1] * lpoly.plg[2][0] + p[2] * lpoly.plg[4][0] +
		p[22] * lpoly.plg[6][0] + p[26] * lpoly.plg[1][0] +
		p[14] * lpoly.plg[3][0] + p[59] * lpoly.plg[5][0];
	// symmetrical annual
	t[2] = (p[18] + p[47] * lpoly.plg[2][0] + p[29] * lpoly.plg[4][0]) * cd31;
	// symmetrical semiannual
	t[3] = (p[15] + p[16] * lpoly.plg[2][0] + p[30] * lpoly.plg[4][0]) * cd17;
	// asymmetrical annual
	t[4] = (p[9] * lpoly.plg[1][0] + p[10] * lpoly.plg[3][0] +
		p[20] * lpoly.plg[5][0]) * cd13;
	// asymmetrical semiannual
	t[5] = (p[37] * lpoly.plg[1][0]) * cd38;
	// diurnal
	double t71, t72;
	if (csw.sw[6] != 0) {
		t71 = p[11] * lpoly.plg[2][1] * cd13 * csw.swc[4];
		t72 = p[12] * lpoly.plg[2][1] * cd13 * csw.swc[4];
			t[6] =((p[3] * lpoly.plg[1][1] + p[4] * lpoly.plg[3][1] + t71) * lpoly.ctloc +
				(p[6] * lpoly.plg[1][1] + p[7] * lpoly.plg[3][1] + t72) * lpoly.stloc);
	}
	// semidiurnal
	double t81, t82;
	if (csw.sw[7] != 0) {
		t81 = (p[23] * lpoly.plg[3][2] + p[35] * lpoly.plg[5][2]) * cd13 * csw.swc[4];
		t82 = (p[33] * lpoly.plg[3][2] + p[36] * lpoly.plg[5][2]) * cd13 * csw.swc[4];
		t[7] =((p[5] * lpoly.plg[2][2] + p[41] * lpoly.plg[4][2] + t81) * lpoly.c2tloc +
			(p[8] * lpoly.plg[2][2] + p[42] * lpoly.plg[4][2] + t82) * lpoly.s2tloc);
	}
	// terdiurnal
	if (csw.sw[13] != 0)
		t[13] = p[39] * lpoly.plg[3][3] * lpoly.s3tloc +
			p[40] * lpoly.plg[3][3] * lpoly.c3tloc;
	// magnetic activity
	if (csw.sw[8] != 0) {
		if (csw.sw[8] == 1)
			t[8] = lpoly.apdf * (p[32] + p[45] * lpoly.plg[2][0] * csw.swc[1]);
		else if (csw.sw[8] == -1)
			t[8] = (p[50] * lpoly.apt[0] + 
				p[96] * lpoly.plg[2][0] * lpoly.apt[0] * csw.swc[1]);
	}
	double dgtr = 1.74533e-2;
	// longitudinal
	if ((csw.sw[9] != 0) && (csw.sw[10] != 0) && (lpoly.xlong > -1000.)) {
		t[10] = (1 + lpoly.plg[1][0] * (p[80] * csw.swc[4] * cos(dr * (lpoly.day - p[81])) +
			p[85] * csw.swc[5] * cos(2 * dr * (lpoly.day - p[86]))) +
			p[83] * csw.swc[2] * cos(dr * (lpoly.day - p[84])) +
			p[87] * csw.swc[3] * cos(2 * dr * (lpoly.day - p[88]))) *
			((p[64] * lpoly.plg[2][1] + p[65] * lpoly.plg[4][1] + p[66] * lpoly.plg[6][1] +
				p[74] * lpoly.plg[1][1] + p[75] * lpoly.plg[3][1] + p[76] * lpoly.plg[5][1]) *
				cos(dgtr * lpoly.xlong) + (p[90] * lpoly.plg[2][1] +
					p[91] * lpoly.plg[4][1] + p[92] * lpoly.plg[6][1] +
					p[77] * lpoly.plg[1][1] + p[78] * lpoly.plg[3][1] +
					p[79] * lpoly.plg[5][1]) * sin(dgtr * lpoly.xlong));
	}
	double tt = 0;
	for (size_t i = 0; i < 14; i++)
		tt += abs(csw.sw[i]) * t[i];
	return tt;
}