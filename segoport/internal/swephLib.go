package internal

import (
	"bufio"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

/* DeltaT = Ephemeris Time - Universal Time, in days.
 * Before 1955 we use the data developed by Stephenson, Morrison, and Hohenkerk (2016),
 * 1955 - today + a couple of years:
 * ---------------------------------
 * The tabulated values of deltaT from the Astronomical Almaanc (AA 1997 etc. pp. K8-K9) are used. Some more recent
 * values have been taken from IERS (http://maia.usno.navy.mil/ser7/deltat.data).
 * Bessel's interpolation formula is implemented to obtain fourth order interpolated values at intermediate times.
 * The values are adjusted depending on the ephemeris used and its inherent value of secular tidal acceleration ndot.
 *
 * future:
 * ---------------------------------
 * For the time after the last tabulated value, we use the formula of Stephenson (1997; p. 507), with a modification
 * that avoids a jump at the end of the tabulated period. A linear term is added that makes a slow transition from the
 * table to the formula over a period of 100 years. (Need not be updated, when table will be enlarged.)
 *
 * References:
 *
 * Stephenson, F. R., and L. V. Morrison, "Long-term changes in the rotation of the Earth: 700 B.C. to A.D. 1980,"
 * Philosophical Transactions of the Royal Society of London Series A 313, 47-70 (1984)
 *
 * Borkowski, K. M., "ELP2000-85 and the Dynamical Time - Universal Time relation," Astronomy and Astrophysics
 * 205, L8-L10 (1988)
 * Borkowski's formula is derived from partly doubtful eclipses going back to 2137 BC and uses lunar position based on
 * tidal coefficient of -23.9 arcsec/cy^2.
 *
 * Chapront-Touze, Michelle, and Jean Chapront, _Lunar Tables and Programs from 4000 B.C. to A.D. 8000_, Willmann-Bell
 * 1991. Their table agrees with the one here, but the entries are rounded to the nearest whole second.
 *
 * Stephenson, F. R., and M. A. Houlden, _Atlas of Historical Eclipse Maps_, Cambridge U. Press (1986)
 *
 * Stephenson, F.R. & Morrison, L.V., "Long-Term Fluctuations in the Earth's Rotation: 700 BC to AD 1990", Philosophical
 * Transactions of the Royal Society of London, Ser. A, 351 (1995), 165-202.
 *
 * Stephenson, F. Richard, _Historical Eclipses and Earth's Rotation_, Cambridge U. Press (1997)
 *
 * Morrison, L. V., and F.R. Stephenson, "Historical Values of the Earth's Clock Error DT and the Calculation of Eclipses",
 * JHA xxxv (2004), pp.327-336
 *
 * Stephenson, F.R., Morrison, L.V., and Hohenkerk, C.Y., "Measurement of the Earth's Rotation: 720 BC to AD 2015",
 * Royal Society Proceedings A  7 Dec 2016, http://rspa.royalsocietypublishing.org/lookup/doi/10.1098/rspa.2016.0404
 *
 * Table from AA for 1620 through today. Note, Stephenson and Morrison's table starts at the year 1630.
 * The Chapronts' table does not agree with the Almanac prior to 1630. The actual accuracy decreases rapidly prior to 1780.
 *
 * Jean Meeus, Astronomical Algorithms, 2nd edition, 1998.
 *
 * For a comprehensive collection of publications and formulae, see:
 * http://www.phys.uu.nl/~vgent/deltat/deltat_modern.htm
 * http://www.phys.uu.nl/~vgent/deltat/deltat_old.htm
 *
 * For future values of delta t, the following data from the
 * Earth Orientation Department of the US Naval Observatory can be used:
 * (TAI-UTC) from: ftp://maia.usno.navy.mil/ser7/tai-utc.dat
 * (UT1-UTC) from: ftp://maia.usno.navy.mil/ser7/finals.all (cols. 59-68)
 *             or: ftp://ftp.iers.org/products/eop/rapid/standard/finals.data
 * file description in: ftp://maia.usno.navy.mil/ser7/readme.finals
 * Delta T = TAI-UT1 + 32.184 sec = (TAI-UTC) - (UT1-UTC) + 32.184 sec
 *
 * Also, there is the following file:
 * http://maia.usno.navy.mil/ser7/deltat.data, but it is about 3 months behind (on 3 feb 2009); and predictions:
 * http://maia.usno.navy.mil/ser7/deltat.preds
 *
 * Last update of table dt[]: Dieter Koch, 18 dec 2013.
 * ATTENTION: Whenever updating this table, do not forget to adjust
 * the macros TABEND and TABSIZ !
 */

// Constants
const (
	TABSTART              = 1620
	TABEND                = 2028
	TABSIZ                = TABEND - TABSTART + 1
	TABSIZ_SPACE          = TABSIZ + 100
	TAB2_SIZ              = 27
	TAB2_START            = -1000
	TAB2_END              = 1600
	TAB2_STEP             = 100
	TAB97_SIZ             = 43
	TAB97_START           = -500
	TAB97_END             = 1600
	TAB97_STEP            = 50
	LTERM_EQUATION_YSTART = 1820
	LTERM_EQUATION_COEFF  = 32
	//J2000                 = 2451545.0
)

// Global variables
var dt = []float64{
	// 1620.0 - 1659.0
	124.00, 119.00, 115.00, 110.00, 106.00, 102.00, 98.00, 95.00, 91.00, 88.00,
	85.00, 82.00, 79.00, 77.00, 74.00, 72.00, 70.00, 67.00, 65.00, 63.00,
	62.00, 60.00, 58.00, 57.00, 55.00, 54.00, 53.00, 51.00, 50.00, 49.00,
	48.00, 47.00, 46.00, 45.00, 44.00, 43.00, 42.00, 41.00, 40.00, 38.00,
	/* 1660.0 - 1699.0 */
	37.00, 36.00, 35.00, 34.00, 33.00, 32.00, 31.00, 30.00, 28.00, 27.00,
	26.00, 25.00, 24.00, 23.00, 22.00, 21.00, 20.00, 19.00, 18.00, 17.00,
	16.00, 15.00, 14.00, 14.00, 13.00, 12.00, 12.00, 11.00, 11.00, 10.00,
	10.00, 10.00, 9.00, 9.00, 9.00, 9.00, 9.00, 9.00, 9.00, 9.00,
	/* 1700.0 - 1739.0 */
	9.00, 9.00, 9.00, 9.00, 9.00, 9.00, 9.00, 9.00, 10.00, 10.00,
	10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 11.00, 11.00, 11.00,
	11.00, 11.00, 11.00, 11.00, 11.00, 11.00, 11.00, 11.00, 11.00, 11.00,
	11.00, 11.00, 11.00, 11.00, 12.00, 12.00, 12.00, 12.00, 12.00, 12.00,
	/* 1740.0 - 1779.0 */
	12.00, 12.00, 12.00, 12.00, 13.00, 13.00, 13.00, 13.00, 13.00, 13.00,
	13.00, 14.00, 14.00, 14.00, 14.00, 14.00, 14.00, 14.00, 15.00, 15.00,
	15.00, 15.00, 15.00, 15.00, 15.00, 16.00, 16.00, 16.00, 16.00, 16.00,
	16.00, 16.00, 16.00, 16.00, 16.00, 17.00, 17.00, 17.00, 17.00, 17.00,
	/* 1780.0 - 1799.0 */
	17.00, 17.00, 17.00, 17.00, 17.00, 17.00, 17.00, 17.00, 17.00, 17.00,
	17.00, 17.00, 16.00, 16.00, 16.00, 16.00, 15.00, 15.00, 14.00, 14.00,
	/* 1800.0 - 1819.0 */
	13.70, 13.40, 13.10, 12.90, 12.70, 12.60, 12.50, 12.50, 12.50, 12.50,
	12.50, 12.50, 12.50, 12.50, 12.50, 12.50, 12.50, 12.40, 12.30, 12.20,
	/* 1820.0 - 1859.0 */
	12.00, 11.70, 11.40, 11.10, 10.60, 10.20, 9.60, 9.10, 8.60, 8.00,
	7.50, 7.00, 6.60, 6.30, 6.00, 5.80, 5.70, 5.60, 5.60, 5.60,
	5.70, 5.80, 5.90, 6.10, 6.20, 6.30, 6.50, 6.60, 6.80, 6.90,
	7.10, 7.20, 7.30, 7.40, 7.50, 7.60, 7.70, 7.70, 7.80, 7.80,
	/* 1860.0 - 1899.0 */
	7.88, 7.82, 7.54, 6.97, 6.40, 6.02, 5.41, 4.10, 2.92, 1.82,
	1.61, .10, -1.02, -1.28, -2.69, -3.24, -3.64, -4.54, -4.71, -5.11,
	-5.40, -5.42, -5.20, -5.46, -5.46, -5.79, -5.63, -5.64, -5.80, -5.66,
	-5.87, -6.01, -6.19, -6.64, -6.44, -6.47, -6.09, -5.76, -4.66, -3.74,
	/* 1900.0 - 1939.0 */
	-2.72, -1.54, -.02, 1.24, 2.64, 3.86, 5.37, 6.14, 7.75, 9.13,
	10.46, 11.53, 13.36, 14.65, 16.01, 17.20, 18.24, 19.06, 20.25, 20.95,
	21.16, 22.25, 22.41, 23.03, 23.49, 23.62, 23.86, 24.49, 24.34, 24.08,
	24.02, 24.00, 23.87, 23.95, 23.86, 23.93, 23.73, 23.92, 23.96, 24.02,
	/* 1940.0 - 1949.0 */
	24.33, 24.83, 25.30, 25.70, 26.24, 26.77, 27.28, 27.78, 28.25, 28.71,
	/* 1950.0 - 1959.0 */
	29.15, 29.57, 29.97, 30.36, 30.72, 31.07, 31.35, 31.68, 32.18, 32.68,
	/* 1960.0 - 1969.0 */
	33.15, 33.59, 34.00, 34.47, 35.03, 35.73, 36.54, 37.43, 38.29, 39.20,
	/* 1970.0 - 1979.0 */
	/* from 1974 on values (with 4-digit precision) were calculated from IERS data */
	40.18, 41.17, 42.23, 43.37, 44.4841, 45.4761, 46.4567, 47.5214, 48.5344, 49.5862,
	/* 1980.0 - 1989.0 */
	50.5387, 51.3808, 52.1668, 52.9565, 53.7882, 54.3427, 54.8713, 55.3222, 55.8197, 56.3000,
	/* 1990.0 - 1999.0 */
	56.8553, 57.5653, 58.3092, 59.1218, 59.9845, 60.7854, 61.6287, 62.2951, 62.9659, 63.4673,
	/* 2000.0 - 2009.0 */
	63.8285, 64.0908, 64.2998, 64.4734, 64.5736, 64.6876, 64.8452, 65.1464, 65.4574, 65.7768,
	/* 2010.0 - 2018.0 */
	66.0699, 66.3246, 66.6030, 66.9069, 67.2810, 67.6439, 68.1024, 68.5927, 68.9676, 69.2202,
	/* 2020.0 - 2023.0        */
	69.3612, 69.3593, 69.2945, 69.1833,
	/* Extrapolated values:
	 * 2024 - 2028 */
	69.10, 69.00, 68.90, 68.80, 68.80,
}

var dt2 = []int32{
	// -1000 to 1600
	25400, 23700, 22000, 21000, 19040, 17190, 15530, 14080, 12790, 11640,
	10580, 9600, 8640, 7680, 6700, 5710, 4740, 3810, 2960, 2200,
	1570, 1090, 740, 490, 320, 200, 120,
}

var dt97 = []int32{
	// -500 to 1600
	16800, 16000, 15300, 14600, 14000, 13400, 12800, 12200, 11600, 11100,
	10600, 10100, 9600, 9100, 8600, 8200, 7700, 7200, 6700, 6200,
	5700, 5200, 4700, 4300, 3800, 3400, 3000, 2600, 2200, 1900,
	1600, 1350, 1100, 900, 750, 600, 470, 380, 300, 230,
	180, 140, 110,
}

const NDTCF16 = 54

// dtcf16 contains spline coefficients from Stephenson, Morrison, & Hohenkerk (2016)
// representing delta T measurements from 720 BC to AD 2015.
// Source: http://rspa.royalsocietypublishing.org/lookup/doi/10.1098/rspa.2016.0404
// Year numbers have been replaced by Julian day numbers
var dtcf16 = [NDTCF16][6]float64{
	//    JD begin      JD end       Coefficient1  Coefficient2  Coefficient3  Coefficient4    Year range
	/*00*/ {1458085.5, 1867156.5, 20550.593, -21268.478, 11863.418, -4541.129}, // -720 to  400
	/*01*/ {1867156.5, 2086302.5, 6604.404, -5981.266, -505.093, 1349.609}, //  400 to 1000
	/*02*/ {2086302.5, 2268923.5, 1467.654, -2452.187, 2460.927, -1183.759}, // 1000 to 1500
	/*03*/ {2268923.5, 2305447.5, 292.635, -216.322, -43.614, 56.681}, // 1500 to 1600
	/*04*/ {2305447.5, 2323710.5, 89.380, -66.754, 31.607, -10.497}, // 1600 to 1650
	/*05*/ {2323710.5, 2349276.5, 43.736, -49.043, 0.227, 15.811}, // 1650 to 1720
	/*06*/ {2349276.5, 2378496.5, 10.730, -1.321, 62.250, -52.946}, // 1720 to 1800
	/*07*/ {2378496.5, 2382148.5, 18.714, -4.457, -1.509, 2.507}, // 1800 to 1810
	/*08*/ {2382148.5, 2385800.5, 15.255, 0.046, 6.012, -4.634}, // 1810 to 1820
	/*09*/ {2385800.5, 2389453.5, 16.679, -1.831, -7.889, 3.799}, // 1820 to 1830
	/*10*/ {2389453.5, 2393105.5, 10.758, -6.211, 3.509, -0.388}, // 1830 to 1840
	/*11*/ {2393105.5, 2396758.5, 7.668, -0.357, 2.345, -0.338}, // 1840 to 1850
	/*12*/ {2396758.5, 2398584.5, 9.317, 1.659, 0.332, -0.932}, // 1850 to 1855
	/*13*/ {2398584.5, 2400410.5, 10.376, -0.472, -2.463, 1.596}, // 1855 to 1860
	/*14*/ {2400410.5, 2402237.5, 9.038, -0.610, 2.325, -2.497}, // 1860 to 1865
	/*15*/ {2402237.5, 2404063.5, 8.256, -3.450, -5.166, 2.729}, // 1865 to 1870
	/*16*/ {2404063.5, 2405889.5, 2.369, -5.596, 3.020, -0.919}, // 1870 to 1875
	/*17*/ {2405889.5, 2407715.5, -1.126, -2.312, 0.264, -0.037}, // 1875 to 1880
	/*18*/ {2407715.5, 2409542.5, -3.211, -1.894, 0.154, 0.562}, // 1880 to 1885
	/*19*/ {2409542.5, 2411368.5, -4.388, 0.101, 1.841, -1.438}, // 1885 to 1890
	/*20*/ {2411368.5, 2413194.5, -3.884, -0.531, -2.473, 1.870}, // 1890 to 1895
	/*21*/ {2413194.5, 2415020.5, -5.017, 0.134, 3.138, -0.232}, // 1895 to 1900
	/*22*/ {2415020.5, 2416846.5, -1.977, 5.715, 2.443, -1.257}, // 1900 to 1905
	/*23*/ {2416846.5, 2418672.5, 4.923, 6.828, -1.329, 0.720}, // 1905 to 1910
	/*24*/ {2418672.5, 2420498.5, 11.142, 6.330, 0.831, -0.825}, // 1910 to 1915
	/*25*/ {2420498.5, 2422324.5, 17.479, 5.518, -1.643, 0.262}, // 1915 to 1920
	/*26*/ {2422324.5, 2424151.5, 21.617, 3.020, -0.856, 0.008}, // 1920 to 1925
	/*27*/ {2424151.5, 2425977.5, 23.789, 1.333, -0.831, 0.127}, // 1925 to 1930
	/*28*/ {2425977.5, 2427803.5, 24.418, 0.052, -0.449, 0.142}, // 1930 to 1935
	/*29*/ {2427803.5, 2429629.5, 24.164, -0.419, -0.022, 0.702}, // 1935 to 1940
	/*30*/ {2429629.5, 2431456.5, 24.426, 1.645, 2.086, -1.106}, // 1940 to 1945
	/*31*/ {2431456.5, 2433282.5, 27.050, 2.499, -1.232, 0.614}, // 1945 to 1950
	/*32*/ {2433282.5, 2434378.5, 28.932, 1.127, 0.220, -0.277}, // 1950 to 1953
	/*33*/ {2434378.5, 2435473.5, 30.002, 0.737, -0.610, 0.631}, // 1953 to 1956
	/*34*/ {2435473.5, 2436569.5, 30.760, 1.409, 1.282, -0.799}, // 1956 to 1959
	/*35*/ {2436569.5, 2437665.5, 32.652, 1.577, -1.115, 0.507}, // 1959 to 1962
	/*36*/ {2437665.5, 2438761.5, 33.621, 0.868, 0.406, 0.199}, // 1962 to 1965
	/*37*/ {2438761.5, 2439856.5, 35.093, 2.275, 1.002, -0.414}, // 1965 to 1968
	/*38*/ {2439856.5, 2440952.5, 37.956, 3.035, -0.242, 0.202}, // 1968 to 1971
	/*39*/ {2440952.5, 2442048.5, 40.951, 3.157, 0.364, -0.229}, // 1971 to 1974
	/*40*/ {2442048.5, 2443144.5, 44.244, 3.198, -0.323, 0.172}, // 1974 to 1977
	/*41*/ {2443144.5, 2444239.5, 47.291, 3.069, 0.193, -0.192}, // 1977 to 1980
	/*42*/ {2444239.5, 2445335.5, 50.361, 2.878, -0.384, 0.081}, // 1980 to 1983
	/*43*/ {2445335.5, 2446431.5, 52.936, 2.354, -0.140, -0.166}, // 1983 to 1986
	/*44*/ {2446431.5, 2447527.5, 54.984, 1.577, -0.637, 0.448}, // 1986 to 1989
	/*45*/ {2447527.5, 2448622.5, 56.373, 1.649, 0.709, -0.277}, // 1989 to 1992
	/*46*/ {2448622.5, 2449718.5, 58.453, 2.235, -0.122, 0.111}, // 1992 to 1995
	/*47*/ {2449718.5, 2450814.5, 60.677, 2.324, 0.212, -0.315}, // 1995 to 1998
	/*48*/ {2450814.5, 2451910.5, 62.899, 1.804, -0.732, 0.112}, // 1998 to 2001
	/*49*/ {2451910.5, 2453005.5, 64.082, 0.675, -0.396, 0.193}, // 2001 to 2004
	/*50*/ {2453005.5, 2454101.5, 64.555, 0.463, 0.184, -0.008}, // 2004 to 2007
	/*51*/ {2454101.5, 2455197.5, 65.194, 0.809, 0.161, -0.101}, // 2007 to 2010
	/*52*/ {2455197.5, 2456293.5, 66.063, 0.828, -0.142, 0.168}, // 2010 to 2013
	/*53*/ {2456293.5, 2457388.5, 66.917, 1.046, 0.360, -0.282}, // 2013 to 2016
}

// Helper functions
//func adjustForTidAcc(ans float64, Y float64, tidAcc float64, tidalCoeff int, flag bool) float64 {
//	// TODO Placeholder for the actual implementation
//	return ans
//}

// const DEMO = 0

/* Transpiled from calc_deltat, original comments:
 * returns DeltaT (ET - UT) in days
 * double tjd 	= 	julian day in UT
 * delta t is adjusted to the tidal acceleration that is compatible with the ephemeris flag contained in iflag and with
 * the ephemeris files made accessible through swe_set_ephe_path() or swe_set_jplfile().
 * If iflag = -1, then the default tidal acceleration is used (i.e. that of DE431).
 */

// CalcDeltat calculates the value for deltaT/
// Input tjd : Julian day, iflag: flag for ephemeris file, if -1 uses default tidal acceleration of DE431
func CalcDeltat(tjd float64, iflag int32, deltat *float64, serr *string) int32 {
	var ans float64 = 0
	var B, Y, Ygreg, dd float64
	var iy int
	var retc int32

	deltaModel := sweData.AstroModels[SE_MODEL_DELTAT]
	var tidAcc float64
	var denum, denumret int32
	epheflag := iflag & SEFLG_EPHMASK
	otherflag := iflag & ^SEFLG_EPHMASK

	if deltaModel == 0 {
		deltaModel = SEMOD_DELTAT_DEFAULT
	}

	// with iflag == -1, we use default tid_acc
	if iflag == -1 {
		// for default tid_acc
		retc = SwiGetTidAcc(tjd, 0, 9999, &denumret, &tidAcc, serr)
	} else {
		// otherwise we use tid_acc consistent with epheflag
		denum = swed.jpldenum
		if (epheflag & SEFLG_SWIEPH) != 0 {
			denum = sweData.Fidat[SEI_FILE_MOON].SwephDenum
		}
		if SwiInitSwedIfStart() && (epheflag&SEFLG_MOSEPH) == 0 {
			if serr != nil {
				*serr = "Please call swe_set_ephe_path() or swe_set_jplfile() before calling swe_deltat_ex()"
			}
			retc = SwiSetTidAcc(tjd, epheflag, denum, nil) // _set_ saves tid_acc in swed
		} else {
			retc = SwiSetTidAcc(tjd, epheflag, denum, serr) // _set_ saves tid_acc in swed
		}
		tidAcc = swed.tidAcc
	}

	iflag = otherflag | retc
	Y = 2000.0 + (tjd-J2000)/365.25
	Ygreg = 2000.0 + (tjd-J2000)/365.2425

	// Stephenson/Morrison/Hohenkerk 2016
	if deltaModel == SEMOD_DELTAT_STEPHENSON_ETC_2016 && tjd < 2435108.5 {
		*deltat = deltaStephensonEtc2016(tjd, tidAcc)
		if tjd >= 2434108.5 {
			*deltat += (1.0 - (2435108.5-tjd)/1000.0) * 0.6610218 / 86400.0
		}
		return iflag
	}

	// Espenak & Meeus 2006
	if deltaModel == SEMOD_DELTAT_ESPENAK_MEEUS_2006 && tjd < 2317746.13090277789 {
		*deltat = deltaEspanakMeeus1620(tjd, tidAcc)
		return iflag
	}

	// Stephenson & Morrison 2004
	if deltaModel == SEMOD_DELTAT_STEPHENSON_MORRISON_2004 && Y < TABSTART {
		if Y < TAB2_END {
			*deltat = deltatStephensonMorrison2004_1600(tjd, tidAcc)
			return iflag
		} else if Y >= TAB2_END {
			B = TABSTART - TAB2_END
			iy = (TAB2_END - TAB2_START) / TAB2_STEP
			dd = (Y - TAB2_END) / B
			ans = float64(dt2[iy]) + dd*(dt[0]-float64(dt2[iy]))
			ans = adjustForTidacc(ans, Ygreg, tidAcc, SE_TIDAL_26, false)
			*deltat = ans / 86400.0
			return iflag
		}
	}

	// Stephenson 1997
	if deltaModel == SEMOD_DELTAT_STEPHENSON_1997 && Y < TABSTART {
		if Y < TAB97_END {
			*deltat = deltatStephensonMorrison1997_1600(tjd, tidAcc)
			return iflag
		} else if Y >= TAB97_END {
			B = TABSTART - TAB97_END
			iy = int((TAB97_END - TAB97_START) / TAB97_STEP)
			dd = (Y - TAB97_END) / B
			ans = float64(dt97[iy]) + dd*(dt[0]-float64(dt97[iy]))
			ans = adjustForTidacc(ans, Ygreg, tidAcc, SE_TIDAL_26, false)
			*deltat = ans / 86400.0
			return iflag
		}
	}

	// Stephenson/Morrison 1984 with Borkowski 1988
	if deltaModel == SEMOD_DELTAT_STEPHENSON_MORRISON_1984 && Y < TABSTART {
		if Y >= 948.0 {
			B = 0.01 * (Y - 2000.0)
			ans = (23.58*B+100.3)*B + 101.6
		} else {
			B = 0.01*(Y-2000.0) + 3.75
			ans = 35.0*B*B + 40.0
		}
		*deltat = ans / 86400.0
		return iflag
	}

	// 1620 - today + a few years (tabend)
	if Y >= TABSTART {
		*deltat = deltaAA(tjd, tidAcc)
		return iflag
	}

	*deltat = ans / 86400.0
	return iflag
}

// SweDeltatEx calculates delta t with extended parameters
func SweDeltatEx(tjd float64, iflag int32, serr *string) float64 {
	// Assuming sweData is a global struct with these fields defined somewhere
	if sweData.DeltaTUserdefIsSet {
		return sweData.DeltaTUserdef
	}

	if serr != nil {
		*serr = ""
	}

	// Assuming calcDeltat is defined elsewhere
	var deltat float64
	CalcDeltat(tjd, iflag, &deltat, serr)
	return deltat
}

// SweDeltat calculates delta t with default parameters
func SweDeltat(tjd float64) float64 {
	// Assuming swiGuessEpheFlag is defined elsewhere
	iflag := SwiGuessEpheFlag()
	return SweDeltatEx(tjd, iflag, nil)
}

/* Transpiled from deltat_aa. Original comments:
 * The tabulated values of deltaT, in hundredths of a second, were taken from The Astronomical Almanac 1997etc.,
 * pp. K8-K9. Some more recent values are taken from IERS http://maia.usno.navy.mil/ser7/deltat.data .
 * Bessel's interpolation formula is implemented to obtain fourth order interpolated values at intermediate times.
 * The values are adjusted depending on the ephemeris used and its inherent value of secular tidal acceleration ndot.
 * Note by Dieter Jan. 2017:
 * Bessel interpolation assumes equidistant sampling points. However the sampling points are not equidistant, because
 * they are for first January of every year and years can have either 365 or 366 days. The interpolation uses
 * a step width of 365.25 days. As a consequence, in three out of four years the interpolation does not reproduce the
 * exact values of the sampling points on the days they refer to.  */

func deltaAA(tjd, tidAcc float64) float64 {
	var ans, ans2, ans3 float64
	var p, B, B2, Y, dd float64
	d := make([]float64, 6)

	// read additional values from swedelta.txt
	tabsiz := initDt()
	tabend := TABSTART + tabsiz - 1
	deltaModel := sweData.AstroModels[SE_MODEL_DELTAT]
	if deltaModel == 0 {
		deltaModel = SEMOD_DELTAT_DEFAULT
	}

	Y = 2000.0 + (tjd-2451544.5)/365.25
	if Y <= float64(tabend) {
		// Index into the table
		p = math.Floor(Y)
		iy := int(p - TABSTART)

		// Zeroth order estimate is value at start of year
		ans = dt[iy]
		k := iy + 1
		if k >= tabsiz {
			goto done // No data, can't go on
		}

		// The fraction of tabulation interval
		p = Y - p
		// First order interpolated value
		ans += p * (dt[k] - dt[iy])
		if (iy-1 < 0) || (iy+2 >= tabsiz) {
			goto done // can't do second differences
		}

		// Make table of first differences
		k = iy - 2
		for i := 0; i < 5; i++ {
			if (k < 0) || (k+1 >= tabsiz) {
				d[i] = 0
			} else {
				d[i] = dt[k+1] - dt[k]
			}
			k++
		}

		// Compute second differences
		for i := 0; i < 4; i++ {
			d[i] = d[i+1] - d[i]
		}

		B = 0.25 * p * (p - 1.0)
		ans += B * (d[1] + d[2])

		if iy+2 >= tabsiz {
			goto done
		}

		// Compute third differences
		for i := 0; i < 3; i++ {
			d[i] = d[i+1] - d[i]
		}
		B = 2.0 * B / 3.0
		ans += (p - 0.5) * B * d[1]

		if (iy-2 < 0) || (iy+3 > tabsiz) {
			goto done
		}

		// Compute fourth differences
		for i := 0; i < 2; i++ {
			d[i] = d[i+1] - d[i]
		}
		B = 0.125 * B * (p + 1.0) * (p - 2.0)
		ans += B * (d[0] + d[1])

	done:
		ans = adjustForTidacc(ans, Y, tidAcc, SE_TIDAL_26, false)
		return ans / 86400.0
	}

	// today - future:
	// 3rd degree polynomial based on data given by
	// Stephenson/Morrison/Hohenkerk 2016
	if deltaModel == SEMOD_DELTAT_STEPHENSON_ETC_2016 {
		B = (Y - 2000)
		if Y < 2500 {
			ans = B*B*B*121.0/30000000.0 + B*B/1250.0 + B*521.0/3000.0 + 64.0
			// for slow transition from tablulated data
			B2 = float64(tabend - 2000)
			ans2 = B2*B2*B2*121.0/30000000.0 + B2*B2/1250.0 + B2*521.0/3000.0 + 64.0
		} else {
			// we use a parable after 2500
			B = 0.01 * (Y - 2000)
			ans = B*B*32.5 + 42.5
		}
	} else {
		// Formula Stephenson (1997; p. 507),
		// with modification to avoid jump at end of AA table,
		// similar to what Meeus 1998 had suggested.
		// Slow transition within 100 years.
		B = 0.01 * (Y - 1820)
		ans = -20 + 31*B*B
		// for slow transition from tablulated data
		B2 = 0.01 * (float64(tabend) - 1820)
		ans2 = -20 + 31*B2*B2
	}

	// slow transition from tabulated values to Stephenson formula:
	if Y <= float64(tabend+100) {
		ans3 = dt[tabsiz-1]
		dd = (ans2 - ans3)
		ans += dd * (Y - float64(tabend+100)) * 0.01
	}

	return ans / 86400.0
}

func deltatLongtermMorrisonStephenson(tjd float64) float64 {
	Ygreg := 2000.0 + (tjd-J2000)/365.2425
	u := (Ygreg - 1820) / 100.0
	return (-20 + 32*u*u)
}

func deltatStephensonMorrison1997_1600(tjd, tidAcc float64) float64 {
	var ans, ans2, ans3 float64
	var p, B, Y, dd float64

	Y = 2000.0 + (tjd-J2000)/365.25

	// before -500:
	// formula by Stephenson (1997; p. 508) but adjusted to fit the starting
	// point of table dt97 (Stephenson 1997).
	if Y < TAB97_START {
		B = (Y - 1735) * 0.01
		ans = -20 + 35*B*B
		ans = adjustForTidacc(ans, Y, tidAcc, SE_TIDAL_26, false)

		// transition from formula to table over 100 years
		if Y >= TAB97_START-100 {
			// starting value of table dt97:
			ans2 = adjustForTidacc(float64(dt97[0]), TAB97_START, tidAcc, SE_TIDAL_26, false)
			// value of formula at epoch TAB97_START
			B = (TAB97_START - 1735) * 0.01
			ans3 = -20 + 35*B*B
			ans3 = adjustForTidacc(ans3, Y, tidAcc, SE_TIDAL_26, false)
			dd = ans3 - ans2
			B = (Y - (TAB97_START - 100)) * 0.01
			// fit to starting point of table dt97
			ans = ans - dd*B
		}
	}

	// between -500 and 1600:
	// linear interpolation between values of table dt97 (Stephenson 1997)
	if Y >= TAB97_START && Y < TAB2_END {
		p = math.Floor(Y)
		iy := int((p - TAB97_START) / 50.0)
		dd = (Y - (TAB97_START + 50*float64(iy))) / 50.0
		ans = float64(dt97[iy]) + float64(dt97[iy+1]-dt97[iy])*dd
		// correction for tidal acceleration used by our ephemeris
		ans = adjustForTidacc(ans, Y, tidAcc, SE_TIDAL_26, false)
	}

	ans /= 86400.0
	return ans
}

func deltatStephensonMorrison2004_1600(tjd, tidAcc float64) float64 {
	var ans, ans2, ans3 float64
	var p, B, dd float64
	var tjd0 float64

	Y := 2000.0 + (tjd-J2000)/365.2425

	// before -1000:
	// formula by Stephenson & Morrison (2004; p. 335) but adjusted to fit the
	// starting point of table dt2.
	if Y < TAB2_START {
		ans = deltatLongtermMorrisonStephenson(tjd)
		ans = adjustForTidacc(ans, Y, tidAcc, SE_TIDAL_26, false)

		// transition from formula to table over 100 years
		if Y >= TAB2_START-100 {
			// starting value of table dt2:
			ans2 = adjustForTidacc(float64(dt2[0]), TAB2_START, tidAcc, SE_TIDAL_26, false)

			// value of formula at epoch TAB2_START
			tjd0 = (TAB2_START-2000)*365.2425 + J2000
			ans3 = deltatLongtermMorrisonStephenson(tjd0)
			ans3 = adjustForTidacc(ans3, Y, tidAcc, SE_TIDAL_26, false)
			dd = ans3 - ans2
			B = (Y - (TAB2_START - 100)) * 0.01
			// fit to starting point of table dt2
			ans = ans - dd*B
		}
	}

	// between -1000 and 1600:
	// linear interpolation between values of table dt2 (Stephenson & Morrison 2004)
	if Y >= TAB2_START && Y < TAB2_END {
		Yjul := 2000 + (tjd-2451557.5)/365.25
		p = math.Floor(Yjul)
		iy := int((p - TAB2_START) / TAB2_STEP)
		dd = (Yjul - (TAB2_START + TAB2_STEP*float64(iy))) / TAB2_STEP
		ans = float64(dt2[iy]) + float64(dt2[iy+1]-dt2[iy])*dd
		// correction for tidal acceleration used by our ephemeris
		ans = adjustForTidacc(ans, Y, tidAcc, SE_TIDAL_26, false)
	}

	ans /= 86400.0
	return ans
}

func deltaStephensonEtc2016(tjd, tidAcc float64) float64 {
	var t, dt float64
	irec := -1
	ygreg := 2000.0 + (tjd-J2000)/365.2425

	// after the year -720 get value from spline curve
	for i := 0; i < NDTCF16; i++ {
		if tjd < dtcf16[i][0] {
			break
		}
		if tjd < dtcf16[i][1] {
			irec = i
			break
		}
	}

	if irec >= 0 {
		t = (tjd - dtcf16[irec][0]) / (dtcf16[irec][1] - dtcf16[irec][0])
		dt = dtcf16[irec][2] + dtcf16[irec][3]*t + dtcf16[irec][4]*t*t + dtcf16[irec][5]*t*t*t
	} else if ygreg < -720 {
		// for earlier epochs, use long term parabola
		t = (ygreg - 1825) / 100.0
		dt = -320 + 32.5*t*t
		dt -= 179.7337208 // to make curve continous on 1 Jan -720 (D. Koch)
	} else {
		// future
		t = (ygreg - 1825) / 100.0
		dt = -320 + 32.5*t*t
		dt += 269.4790417 // to make curve continous on 1 Jan 2016 (D. Koch)
	}

	// The parameter adjust_after_1955 must be TRUE here, because the
	// Stephenson 2016 curve is based on occultation data alone,
	// not on IERS data.
	// Note, however, the current function deltaStephensonEtc2016()
	// is called only for dates before 1 Jan 1955.
	dt = adjustForTidacc(dt, ygreg, tidAcc, SE_TIDAL_STEPHENSON_2016, true)
	dt /= 86400.0
	return dt
}

func deltaEspanakMeeus1620(tjd, tidAcc float64) float64 {
	var ans float64
	ygreg := 2000.0 + (tjd-J2000)/365.2425
	var u float64

	switch {
	case ygreg < -500:
		ans = deltatLongtermMorrisonStephenson(tjd)

	case ygreg < 500:
		u = ygreg / 100.0
		ans = (((((0.0090316521*u+0.022174192)*u-0.1798452)*u-5.952053)*u+33.78311)*u-1014.41)*u + 10583.6

	case ygreg < 1600:
		u = (ygreg - 1000) / 100.0
		ans = (((((0.0083572073*u-0.005050998)*u-0.8503463)*u+0.319781)*u+71.23472)*u-556.01)*u + 1574.2

	case ygreg < 1700:
		u = ygreg - 1600
		ans = 120 - 0.9808*u - 0.01532*u*u + u*u*u/7129.0

	case ygreg < 1800:
		u = ygreg - 1700
		ans = (((-u/1174000.0+0.00013336)*u-0.0059285)*u+0.1603)*u + 8.83

	case ygreg < 1860:
		u = ygreg - 1800
		ans = ((((((0.000000000875*u-0.0000001699)*u+0.0000121272)*u-0.00037436)*u+0.0041116)*u+0.0068612)*u-0.332447)*u + 13.72

	case ygreg < 1900:
		u = ygreg - 1860
		ans = ((((u/233174.0-0.0004473624)*u+0.01680668)*u-0.251754)*u+0.5737)*u + 7.62

	case ygreg < 1920:
		u = ygreg - 1900
		ans = (((-0.000197*u+0.0061966)*u-0.0598939)*u+1.494119)*u - 2.79

	case ygreg < 1941:
		u = ygreg - 1920
		ans = 21.20 + 0.84493*u - 0.076100*u*u + 0.0020936*u*u*u

	case ygreg < 1961:
		u = ygreg - 1950
		ans = 29.07 + 0.407*u - u*u/233.0 + u*u*u/2547.0

	case ygreg < 1986:
		u = ygreg - 1975
		ans = 45.45 + 1.067*u - u*u/260.0 - u*u*u/718.0

	case ygreg < 2005:
		u = ygreg - 2000
		ans = ((((0.00002373599*u+0.000651814)*u+0.0017275)*u-0.060374)*u+0.3345)*u + 63.86
	}

	ans = adjustForTidacc(ans, ygreg, tidAcc, SE_TIDAL_26, false)
	ans /= 86400.0
	return ans
}

const (
	AS_MAXCH = 256 // Adjust this value according to your needs
	//	TABSTART     = 1620 // Adjust according to your needs
	//	TABSIZ_SPACE = 1000 // Adjust according to your needs
)

type Swed struct {
	initDtDone         bool
	ephepath           string
	swiGuessEpheFlag   int32
	tidAcc             float64
	deltaTUserdefIsSet bool
	isTidAccManual     bool
	jplFileIsOpen      bool
	jpldenum           int32
	deltaTUserdef      float64
	deltaModel         []int
	// ... other fields as needed
}

var (
	swed Swed
	//	dt   []float64 = make([]float64, TABSIZ_SPACE)
)

// swiFopen is a helper function to emulate C's swi_fopen behavior
func swiFopen(prefix int, filename, path string, errMsg *string) (*os.File, error) {
	fullPath := filepath.Join(path, filename)
	file, err := os.Open(fullPath)
	if err != nil {
		return nil, err
	}
	return file, nil
}

func initDt() int {
	if !sweData.InitDtDone {
		sweData.InitDtDone = true

		// Try opening either swe_deltat.txt or sedeltat.txt
		file, err := swiFopen(-1, "swe_deltat.txt", sweData.EphePath, nil)
		if err != nil {
			file, err = swiFopen(-1, "sedeltat.txt", sweData.EphePath, nil)
			if err != nil {
				return TABSIZ // No error message if file is missing
			}
		}
		defer file.Close()

		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			line := scanner.Text()

			// Skip leading whitespace
			line = strings.TrimSpace(line)

			// Skip comments and empty lines
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			// Parse year
			fields := strings.Fields(line)
			if len(fields) < 2 {
				continue
			}

			year, err := strconv.Atoi(fields[0])
			if err != nil {
				continue
			}

			tabIndex := year - TABSTART

			// Table space is limited. no error msg, if exceeded
			if tabIndex >= TABSIZ_SPACE {
				continue
			}

			// Parse delta t value
			deltaT, err := strconv.ParseFloat(fields[1], 64)
			if err != nil {
				continue
			}
			dt[tabIndex] = deltaT
		}

		if err := scanner.Err(); err != nil {
			// Handle error if needed
			return TABSIZ
		}
	}

	// Find table size
	tabsiz := 2001 - TABSTART + 1
	for i := tabsiz - 1; i < TABSIZ_SPACE; i++ {
		if dt[i] == 0 {
			break
		}
		tabsiz++
	}
	tabsiz--

	return tabsiz
}

/* Original comments:
 * Astronomical Almanac table is corrected by adding the expression -0.000091 (ndot + 26)(year-1955)^2  seconds
 * to entries prior to 1955 (AA page K8), where ndot is the secular tidal term in the mean motion of the Moon.
 *
 * Entries after 1955 are referred to atomic time standards and are not affected by errors in Lunar or planetary theory.
 */

// ASBool is a type alias for boolean to match C's AS_BOOL
type ASBool bool

func adjustForTidacc(ans, y, tidAcc, tidAcc0 float64, adjustAfter1955 ASBool) float64 {
	if y < 1955.0 || adjustAfter1955 {
		b := y - 1955.0
		ans += -0.000091 * (tidAcc - tidAcc0) * b * b
	}
	return ans
}

//const (
//	SE_TIDAL_AUTOMATIC    = 0
//	SE_DELTAT_AUTOMATIC   = 0
//	SEFLG_SWIEPH         = 2
//	SEFLG_JPLEPH         = 4
//	SEFLG_MOSEPH         = 8
//	SEFLG_EPHMASK        = 0x000000FF // Adjust value as needed
//	SEI_FILE_MOON        = 0          // Adjust value as needed
//)

// SetTidAcc sets tidal acceleration of the Moon
func SetTidAcc(tAcc float64) {
	if tAcc == SE_TIDAL_AUTOMATIC {
		sweData.TidAcc = SE_TIDAL_DEFAULT
		sweData.IsTidAccManual = false
		return
	}
	sweData.TidAcc = tAcc
	sweData.IsTidAccManual = true
}

// SetDeltaTUserdef sets user-defined delta t
// Probably not in use, ignore for now
//func SetDeltaTUserdef(dt float64) {
//	if dt == SE_DELTAT_AUTOMATIC {
//		sweData.deltaTUserdefIsSet = false
//	} else {
//		sweData.deltaTUserdefIsSet = true
//		sweData.deltaTUserdef = dt
//	}
//}

// GuessEpheFlag guesses ephemeris flag
func GuessEpheFlag() int32 {
	if sweData.JplFileIsOpen {
		return SEFLG_JPLEPH
	}
	return SEFLG_SWIEPH
}

func SwiGetTidAcc(tjdUt float64, iflag int32, denum int32, denumret *int32, tidAcc *float64, serr *string) int32 {
	iflag &= SEFLG_EPHMASK

	if swed.isTidAccManual {
		*tidAcc = swed.tidAcc
		return iflag
	}

	if denum == 0 {
		if (iflag & SEFLG_MOSEPH) != 0 {
			*tidAcc = SE_TIDAL_DE404
			*denumret = 404
			return iflag
		}

		if (iflag & SEFLG_JPLEPH) != 0 {
			if swed.jplFileIsOpen {
				denum = swed.jpldenum
			}
		}

		// SEFLG_SWIEPH wanted or SEFLG_JPLEPH failed:
		if (iflag & SEFLG_SWIEPH) != 0 {
			if sweData.Fidat[SEI_FILE_MOON].Fptr != nil {
				denum = sweData.Fidat[SEI_FILE_MOON].SwephDenum
			}
		}
	}

	switch denum {
	case 200:
		*tidAcc = SE_TIDAL_DE200
	case 403:
		*tidAcc = SE_TIDAL_DE403
	case 404:
		*tidAcc = SE_TIDAL_DE404
	case 405:
		*tidAcc = SE_TIDAL_DE405
	case 406:
		*tidAcc = SE_TIDAL_DE406
	case 421:
		*tidAcc = SE_TIDAL_DE421
	case 422:
		*tidAcc = SE_TIDAL_DE422
	case 430:
		*tidAcc = SE_TIDAL_DE430
	case 431:
		*tidAcc = SE_TIDAL_DE431
	case 440:
		*tidAcc = SE_TIDAL_DE441
	case 441:
		*tidAcc = SE_TIDAL_DE441
	default:
		denum = SE_DE_NUMBER
		*tidAcc = SE_TIDAL_DEFAULT
	}

	*denumret = denum
	iflag &= SEFLG_EPHMASK
	return iflag
}

func SwiSetTidAcc(tjdUt float64, iflag int32, denum int32, serr *string) int32 {
	retc := iflag
	var denumret int32

	// manual tid_acc overrides automatic tid_acc
	if swed.isTidAccManual {
		return retc
	}

	retc = SwiGetTidAcc(tjdUt, iflag, denum, &denumret, &(swed.tidAcc), serr)

	return retc
}

// SetTidAccWithFlag sets tidal acceleration with flags
//func SetTidAccWithFlag(tjdUt float64, iflag, denum int32) (int32, error) {
//	retc := iflag
//
//	// Manual tid_acc overrides automatic tid_acc
//	if sweData.IsTidAccManual {
//		return retc, nil
//	}
//
//	var err error
//	retc, _, sweData.TidAcc, err = SwiGetTidAcc(tjdUt, iflag, denum)
//	if err != nil {
//		return retc, err
//	}
//
//	// Trace functionality could be implemented differently in Go
//	// For example, using a logger
//	//if sweData.traceEnabled {
//	//	logger.Printf("swe_set_tid_acc: %f\n", sweData.tidAcc)
//	//}
//
//	return retc, nil
//}

func SwiGuessEpheFlag() int32 {
	var iflag int32 = SEFLG_SWIEPH

	// if jpl file is open, assume SEFLG_JPLEPH
	if swed.jplFileIsOpen {
		iflag = SEFLG_JPLEPH
	} else {
		iflag = SEFLG_SWIEPH
	}

	return iflag
}
