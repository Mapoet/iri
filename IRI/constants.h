#pragma once

#ifndef CONSTANTS_H
#define CONSTANTS_H
/*C Initialize parameters for COMMON/IGRF1/
C   ERA		EARTH RADIUS (WGS-84: 6371.137 KM) 
C   EREQU   MAJOR HALF AXIS FOR EARTH ELLIPSOID (6378.160 KM)
C   ERPOL   MINOR HALF AXIS FOR EARTH ELLIPSOID (6356.775 KM)
C   AQUAD   SQUARE OF MAJOR HALF AXIS FOR EARTH ELLIPSOID
C   BQUAD   SQUARE OF MINOR HALF AXIS FOR EARTH ELLIPSOID
C   EEXC	Eccentricity of Earth's orbit
C   DIMO    Earth's dipole moment in Gauss 
C ERA, EREQU and ERPOL as recommended by the INTERNATIONAL 
C ASTRONOMICAL UNION .*/
namespace constants {
	const double pi = 3.141592653589793,
		umr = pi / 180,
		humr = pi / 12.,
		dumr = pi / 182.5,
		era = 6371.2,
		erequ = 6378.16,
		erpol = 6356.775,
		aquad = erequ * erequ,
		bquad = erpol * erpol,
		eexc = 0.01675,
		dimo = 0.311653,
		argmax = 88.0;
}
#endif