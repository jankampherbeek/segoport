package internal

import (
	"math"
)

const (
	SE_VERSION = "2.10.03"

	J2000 = 2451545.0        // 2000 January 1.5
	B1950 = 2433282.42345905 // 1950 January 0.923
	J1900 = 2415020.0        // 1900 January 0.5
	B1850 = 2396758.2035810  // 1850 January 16:53

	MPC_CERES  = 1
	MPC_PALLAS = 2
	MPC_JUNO   = 3
	MPC_VESTA  = 4
	MPC_CHIRON = 2060
	MPC_PHOLUS = 5145

	SE_NAME_SUN       = "Sun"
	SE_NAME_MOON      = "Moon"
	SE_NAME_MERCURY   = "Mercury"
	SE_NAME_VENUS     = "Venus"
	SE_NAME_MARS      = "Mars"
	SE_NAME_JUPITER   = "Jupiter"
	SE_NAME_SATURN    = "Saturn"
	SE_NAME_URANUS    = "Uranus"
	SE_NAME_NEPTUNE   = "Neptune"
	SE_NAME_PLUTO     = "Pluto"
	SE_NAME_MEAN_NODE = "mean Node"
	SE_NAME_TRUE_NODE = "true Node"
	SE_NAME_MEAN_APOG = "mean Apogee"
	SE_NAME_OSCU_APOG = "osc. Apogee"
	SE_NAME_INTP_APOG = "intp. Apogee"
	SE_NAME_INTP_PERG = "intp. Perigee"
	SE_NAME_EARTH     = "Earth"
	SE_NAME_CERES     = "Ceres"
	SE_NAME_PALLAS    = "Pallas"
	SE_NAME_JUNO      = "Juno"
	SE_NAME_VESTA     = "Vesta"
	SE_NAME_CHIRON    = "Chiron"
	SE_NAME_PHOLUS    = "Pholus"

	SE_NAME_CUPIDO            = "Cupido"
	SE_NAME_HADES             = "Hades"
	SE_NAME_ZEUS              = "Zeus"
	SE_NAME_KRONOS            = "Kronos"
	SE_NAME_APOLLON           = "Apollon"
	SE_NAME_ADMETOS           = "Admetos"
	SE_NAME_VULKANUS          = "Vulkanus"
	SE_NAME_POSEIDON          = "Poseidon"
	SE_NAME_ISIS              = "Isis"
	SE_NAME_NIBIRU            = "Nibiru"
	SE_NAME_HARRINGTON        = "Harrington"
	SE_NAME_NEPTUNE_LEVERRIER = "Leverrier"
	SE_NAME_NEPTUNE_ADAMS     = "Adams"
	SE_NAME_PLUTO_LOWELL      = "Lowell"
	SE_NAME_PLUTO_PICKERING   = "Pickering"
	SE_NAME_VULCAN            = "Vulcan"
	SE_NAME_WHITE_MOON        = "White Moon"

	// Mathematical constants
	PI    = math.Pi
	TWOPI = 2.0 * PI

	ENDMARK = -99

	SEI_EPSILON  = -2
	SEI_NUTATION = -1
	SEI_EMB      = 0
	SEI_EARTH    = 0
	SEI_SUN      = 0
	SEI_MOON     = 1
	SEI_MERCURY  = 2
	SEI_VENUS    = 3
	SEI_MARS     = 4
	SEI_JUPITER  = 5
	SEI_SATURN   = 6
	SEI_URANUS   = 7
	SEI_NEPTUNE  = 8
	SEI_PLUTO    = 9
	SEI_SUNBARY  = 10 /* barycentric sun */
	SEI_ANYBODY  = 11 /* any asteroid */
	SEI_CHIRON   = 12
	SEI_PHOLUS   = 13
	SEI_CERES    = 14
	SEI_PALLAS   = 15
	SEI_JUNO     = 16
	SEI_VESTA    = 17

	SEI_NPLANETS = 18

	SEI_MEAN_NODE = 0
	SEI_TRUE_NODE = 1
	SEI_MEAN_APOG = 2
	SEI_OSCU_APOG = 3
	SEI_INTP_APOG = 4
	SEI_INTP_PERG = 5

	SEI_NNODE_ETC = 6

	SEI_FLG_HELIO   = 1
	SEI_FLG_ROTATE  = 2
	SEI_FLG_ELLIPSE = 4
	SEI_FLG_EMBHEL  = 8 /* TRUE, if heliocentric earth is given
	 * instead of barycentric sun
	 * i.e. bary sun is computed from
	 * barycentric and heliocentric earth */

	//#if 0
	//SEI_FILE_TEST_ENDIAN=     (97L * 65536L + 98L * 256L + 99L) /*abc*/
	//#endif
	SEI_FILE_TEST_ENDIAN = 0x616263 /* abc*/
	SEI_FILE_BIGENDIAN   = 0
	SEI_FILE_NOREORD     = 0
	SEI_FILE_LITENDIAN   = 1
	SEI_FILE_REORD       = 2

	SEI_FILE_NMAXPLAN   = 50
	SEI_FILE_EFPOSBEGIN = 500

	SE_FILE_SUFFIX = "se1"

	SEI_NEPHFILES = 7
	SEI_CURR_FPOS = -1
	SEI_NMODELS   = 8

	SEI_ECL_GEOALT_MAX = 25000.0
	SEI_ECL_GEOALT_MIN = -500.0

	/* Chiron's orbit becomes chaotic before 720 AD and after 4606 AD, because of close encounters with Saturn.
	 * Accepting a maximum error of 5 degrees, the ephemeris is good between the following dates:
	 */
	CHIRON_START = 1967601.5 /* 1.1.675 */
	CHIRON_END   = 3419437.5 /* 1.1.4650 */

	/* Pholus's orbit is unstable as well, because he sometimes approaches Saturn. Accepting a maximum error of 5 degrees,
	 * the ephemeris is good after the following date:
	 */
	PHOLUS_START = 640648.5  /* 1.1.-2958 jul */
	PHOLUS_END   = 4390617.5 /* 1.1.7309 */

	MOSHPLEPH_START = 625000.5
	MOSHPLEPH_END   = 2818000.5
	MOSHLUEPH_START = 625000.5
	MOSHLUEPH_END   = 2818000.5
	MOSHNDEPH_START = -3100015.5 /* 15 Aug -13200 00:00 ET jul.cal.*/
	MOSHNDEPH_END   = 8000016.5  /* 15 Mar 17191 00:00 ET, greg. cal */

	JPL_DE431_START = -3027215.5
	JPL_DE431_END   = 7930192.5

	//#if FALSE	/*	Alois commented out, not used anywhere  */
	//JPLEPH_START	 625307.5	/* about -3000 (DE406) */
	//JPLEPH_END	2816848.5	/* about  3000 (DE406) */
	//SWIEPH_START	 625614.927151
	//SWIEPH_END	2813641.5
	//ALLEPH_START	MOSHPLEPH_START
	//ALLEPH_END	MOSHPLEPH_END
	//BEG_YEAR       (-3000)
	//END_YEAR       3000
	//#endif

	MAXORD = 40

	NCTIES = 6.0 /* number of centuries per eph. file */

	OK                = 0
	ERR               = -1
	NOT_AVAILABLE     = -2
	BEYOND_EPH_LIMITS = -3

	J_TO_J2000 = 1
	J2000_TO_J = -1

	/* we always use Astronomical Almanac constants, if available */
	MOON_MEAN_DIST = 384400000.0 /* in m, AA 1996, F2 */
	MOON_MEAN_INCL = 5.1453964   /* AA 1996, D2 */
	MOON_MEAN_ECC  = 0.054900489 /* AA 1996, F2 */
	/* SUN_EARTH_MRAT  328900.561400           Su/(Ea+Mo) AA 2006 K7 */
	SUN_EARTH_MRAT  = 332946.050895    /* Su / (Ea only) AA 2006 K7 */
	EARTH_MOON_MRAT = 1 / 0.0123000383 /* AA 2006, K7 */
	//#if 0
	//EARTH_MOON_MRAT 81.30056907419062	/* de431 */
	//#endif
	//#if 0
	//EARTH_MOON_MRAT 81.30056		/* de406 */
	//#endif

	AUNIT  = 1.49597870700e+11 /* au in meters, DE431 */
	CLIGHT = 2.99792458e+8     /* m/s, AA 1996 K6 / DE431 */
	//#if 0
	//HELGRAVCONST    1.32712438e+20		/* G * M(sun), m^3/sec^2, AA 1996 K6 */
	//#endif
	HELGRAVCONST     = 1.32712440017987e+20     /* G * M(sun), m^3/sec^2, AA 2006 K6 */
	GEOGCONST        = 3.98600448e+14           /* G * M(earth) m^3/sec^2, AA 1996 K6 */
	KGAUSS           = 0.01720209895            /* Gaussian gravitational constant K6 */
	SUN_RADIUS       = 959.63 / 3600 * DEGTORAD /*  Meeus germ. p 391 */
	EARTH_RADIUS     = 6378136.6                /* AA 2006 K6 */
	EARTH_OBLATENESS = 1.0 / 298.25642          /* AA 2006 K6 */
	EARTH_ROT_SPEED  = 7.2921151467e-5 * 86400  /* in rad/day, expl. suppl., p 162 */

	LIGHTTIME_AUNIT = 499.0047838362 / 3600.0 / 24.0 /* 8.3167 minutes (days) */
	PARSEC_TO_AUNIT = 206264.8062471                 /* 648000/PI, according to IAU Resolution B2, 2016 */

	/* node of ecliptic measured on ecliptic 2000 */
	SSY_PLANE_NODE_E2000 = 107.582569 * DEGTORAD
	/* node of ecliptic measured on solar system rotation plane */
	SSY_PLANE_NODE = 107.58883388 * DEGTORAD
	/* inclination of ecliptic against solar system rotation plane */
	SSY_PLANE_INCL = 1.578701 * DEGTORAD

	KM_S_TO_AU_CTY       = 21.095  /* km/s to AU/century */
	MOON_SPEED_INTV      = 0.00005 /* 4.32 seconds (in days) */
	PLAN_SPEED_INTV      = 0.0001  /* 8.64 seconds (in days) */
	MEAN_NODE_SPEED_INTV = 0.001
	NODE_CALC_INTV       = 0.0001
	NODE_CALC_INTV_MOSH  = 0.1
	NUT_SPEED_INTV       = 0.0001
	DEFL_SPEED_INTV      = 0.0000005

	SE_LAPSE_RATE = 0.0065 /* deg K / m, for refraction */
)

// SquareSum returns the sum of squares of the first three elements of a slice
func SquareSum(x []float64) float64 {
	return x[0]*x[0] + x[1]*x[1] + x[2]*x[2]
}

// DotProduct returns the dot product of the first three elements of two slices
func DotProduct(x, y []float64) float64 {
	return x[0]*y[0] + x[1]*y[1] + x[2]*y[2]
}

var PNOINT2JPL = []int{
	J_EARTH,
	J_MOON,
	J_MERCURY,
	J_VENUS,
	J_MARS,
	J_JUPITER,
	J_SATURN,
	J_URANUS,
	J_NEPTUNE,
	J_PLUTO,
	J_SUN,
}

const NDIAM = SE_VESTA + 1

// PlaDiam represents planetary radii in meters
var PlaDiam = [NDIAM]float64{
	1392000000.0,   // Sun
	3475000.0,      // Moon
	2439400.0 * 2,  // Mercury
	6051800.0 * 2,  // Venus
	3389500.0 * 2,  // Mars
	69911000.0 * 2, // Jupiter
	58232000.0 * 2, // Saturn
	25362000.0 * 2, // Uranus
	24622000.0 * 2, // Neptune
	1188300.0 * 2,  // Pluto
	0, 0, 0, 0,     // nodes and apogees
	6371008.4 * 2, // Earth
	271370.0,      // Chiron
	290000.0,      // Pholus
	939400.0,      // Ceres
	545000.0,      // Pallas
	246596.0,      // Juno
	525400.0,      // Vesta
}

/* Orignal comment:
 * Ayanamsas
 * For each ayanamsa, there are the following values:
 * t0       epoch of ayanamsa, TDT (can be ET or UT)
 * ayan_t0  ayanamsa value at epoch
 * t0_is_UT true, if t0 is UT
 * prec_offset is the precession model for which the ayanamsha has to be corrected by adding/subtracting a constant offset.
 *          0, if no correction is needed
 *          -1, if correction is unclear or has not been investigated and therefore is not applied
 */

type AyaInit struct {
	T0         float64 // epoch of ayanamsa, TDT (can be ET or UT)
	AyanT0     float64 // ayanamsa value at epoch
	T0IsUT     bool    // true if t0 is UT
	PrecOffset int     // precession model offset
}

// Ayanamsa represents the array of ayanamsa initializations
var Ayanamsa = [SE_NSIDM_PREDEF]AyaInit{
	/* 0: Fagan/Bradley (Default)
	   "The American Sidereal Ephemeris, 1976-2000" (Astro Computing Services, 1981) states on S.V.P. ("Synetic
	   Vernal Point"): "The S.V.P. is the Sidereal longitude of the Vernal Equinox (the Tropical zero-point) in the
	   Fagan-Bradley school of Western Sidereal astrology. It was determined empirically, its mean value being defined
	   as 335°57'28".64 for the epoch 1950.0." Fagan/Firebrace, "Primer of Sidereal Astrology", p. 13:
	   "It was during 1957 that Garth Allen .... experimenting ... But when progressed for the dates of the calamities,
	   all were found by him to be slightly out, the mean error being equivalent to an increase of 0°06'05" in the
	   then-adopted sidereal longitude of the vernal point, determined from Spica in 29 Virgo (i.e. 29°06'05" Virgo;
	   D.K.), and the proper motion having been allowed for. In short, for the epoch 1950.0 he proposed as the mean
	   longitude of the vernal point 335°57'28.64", proper motion being disregarded."
	   If "1950.0" means the standard epoch B1950 = JD 2433282.423, and based on the then-used precession model of
	   Newcomb, this ayanamsha leads to a true position of 29°06'05.965" Virgo, based on Hipparcos position of
	   the star. */
	{2433282.42346, 24.042044444, false, SEMOD_PREC_NEWCOMB},
	/* 1: Standard Lahiri
	   according to program NOVA by Robert Hand: {J1900, 360 - 337.53953}, This corresponds to an ayanamsha 22°27'37.69
	   as given in Indian Ephemeris and Nautical Almanac" 1965, p. 459. Note, however, this value should only with a
	   precession formula where T is measured in tropical centuries. Swiss Ephemeris always uses Julian centuries.
	   The following definition is according to: Calendar Reform Committee 1956; the subtracted value is nutation:
	   {2435553.5, 23.25 - 0.00464207, FALSE},
	   Lahiri (derived from: Indian Astronomical Ephemeris 1989, p. 556; the subtracted value is nutation, according
	   to Wahr 1980) */
	{2435553.5, 23.250182778 - 0.004658035, false, SEMOD_PREC_IAU_1976},
	/* 2: Robert DeLuce (Constellational Astrology ... p. 5; birth of Jesus assumed on 1 Jan 1 BC (= 0) jul.,
	   based on Newcomb precession. {J1900, 360 - 333.58695, FALSE, 0},
	   Ayanamsha was corrected with SE 2.09 as follows: Started at zero ayanamsha epoch with value 0 and run with
	   standard precession. This makes a difference of 22" compared with previous version: */
	{1721057.5, 0, true, 0},
	/* 3: B.V. Raman (Robert Hand) See B.V. Raman, "Hindu Predictive Astrology" (1938, Introduction), pp. 279, 287.
	   This ayanamsha is apparently not based on a valid precession theory (e.g. Newcomb). We cannot reproduce precisely
	   the ayanamsha values on p. 287. */
	{J1900, 360 - 338.98556, false, SEMOD_PREC_NEWCOMB},
	/* 4: Usha/Shashi (Robert Hand) Usha and Shashi, "Hindu Astrological Calculations" (1978, Sagar Publications,
	   New Delhi). We do not have this book. */
	// 4: Usha/Shashi
	{J1900, 360 - 341.33904, false, -1},
	/* 5: Krishnamurti (Robert Hand) K.S. Krishnamurti, "Reader 1", pp. 55-59. Autor does not give precise information.
	   Zero ayanamsha year is said to be 291 CE, and there is an ayanamsha table with arc min precision for 1840 to
	   2000 on p. 58.
	   This ayanamsha reproduces the table quite well, if 1 Jan of each year is taken. (Turn off Newcomb precession in
	   order to verify.) However, D. Senthilathiban believes the values are given for the date of sidereal Aries ingress
	   of each year. ("Study of KP Ayanamsa with Modern Precession Theories", pp. 126f. */
	{J1900, 360 - 337.636111, false, SEMOD_PREC_NEWCOMB},
	/* 6: Djwhal Khool (Graham Dawson),
	   "Channeled" information: Aquarius ingress of VP on 1 July 2117
	   See Philipp Lindsay, “The Beginning of the Age of Aquarius: 2,117 A.D.”,
	   http://esotericastrologer.org/newsletters/the-age-of-aquarius-ray-and-zodiac-cycles/ */
	{J1900, 360 - 333.0369024, false, 0},
	/* 7: Shri Yukteshwar; (David Cochrane) This ayanamsha seems to be wrong. Swami Sri Yukteswar, "The Holy Science",
	   1920 (1949, 1957 and 1977, partly revised), Yogoda Satsanga Society of India.
	   Ayanamsha on the spring equinox 1893 was 20°54'36" (1894 according to the revised edition of 1977) At the same
	   time he believed that this was the distance of the spring equinox from the star Revati, which he put at the
	   initial point of Aries.  Unfortunately, this is wrong, because on that date Revati was actually 18°23' away from
	   the vernal point. The error is explained from the fact that Yukteshwar used the zero ayanamsha year 499 CE
	   and an inaccurate Suryasiddhantic precession rate of 360°/24'000 years = 54 arcsec/year. It is obvious
	   that Yukteshwar actually intended an ayanamsha that starts at the star Revati.  */
	{J1900, 360 - 338.917778, false, -1},
	/* 8: J.N. Bhasin; (David Cochrane)
	   We don't have any sources or detailed information about this ayanamsha. */
	{J1900, 360 - 338.634444, false, -1},
	/* 14 Sept. 2018: the following three ayanamshas have been wrong for many years */
	// 9-11: Babylonian, Kugler
	{1684532.5, -5.66667, true, -1},
	{1684532.5, -4.26667, true, -1},
	{1684532.5, -3.41667, true, -1},
	/* 12: Babylonian, Huber. P. Huber, "Über den Nullpunkt der babylonischen Ekliptik", in: Centaurus 1958, 5, p. 192-208.
	   This ayanamsha had a wrong initial value until 14 Sept. 2018. */
	{1684532.5, -4.46667, true, -1},
	/* 13: Babylonian, Mercier; eta Piscium culminates with zero point */
	{1673941, -5.079167, true, -1},
	/* 14: t0 is defined by Aldebaran at 15 Taurus in year -100 */
	// 14: Babylonian/Aldebaran
	{1684532.5, -4.44138598, true, 0},
	/* 15: Hipparchos */
	{1674484.0, -9.33333, true, -1},
	/* 16: Sassanian */
	{1927135.8747793, 0, true, -1},
	/* 17: Galactic Center at 0 Sagittarius */
	{0, 0, false, 0},
	/* 18: J2000 */
	{J2000, 0, false, 0},
	/* 19: J1900 */
	{J1900, 0, false, 0},
	/* 20: B1950 */
	{B1950, 0, false, 0},
	/* 21: Suryasiddhanta, assuming ingress of mean Sun into Aries at point of mean equinox of date on 21.3.499, near
	noon, Ujjain (75.7684565 E) = 7:30:31.57 UT = 12:33:36 LMT*/
	{1903396.8128654, 0, true, 0},
	/* 22: Suryasiddhanta, assuming ingress of mean Sun into Aries at true position of mean Sun at same epoch */
	{1903396.8128654, -0.21463395, true, 0},
	/* 23: Aryabhata, same date, but UT 6:56:55.57 analogous to 21 */
	{1903396.7895321, 0, true, 0},
	/* 24: Aryabhata, analogous 22 */
	{1903396.7895321, -0.23763238, true, 0},
	/* 25: Suryasiddhanta, Revati/zePsc at polar long. 359°50'*/
	{1903396.8128654, -0.79167046, true, 0},
	/* 26: Suryasiddhanta, Citra/Spica at polar long. 180° */
	{1903396.8128654, 2.11070444, true, 0},
	/* 27: True Citra (Spica exactly at 0 Libra) */
	{0, 0, false, 0},
	/* 28: True Revati (zeta Psc exactly at 29°50' Pisces) */
	{0, 0, false, 0},
	/* 29: True Pushya (delta Cnc exactly a 16 Cancer */
	{0, 0, false, 0},
	/* 30: R. Gil Brand; Galactic Center at golden section between 0 Sco and 0 Aqu;
	   note: 0° Aqu/Leo is the symmetric axis of rulerships */
	{0, 0, false, 0},
	/* 31: Galactic Equator IAU 1958, i.e. galactic/ecliptic intersection point based on galactic coordinate system */
	{0, 0, false, 0},
	/* 32: Galactic Equator True, i.e. galactic/ecliptic intersection point based on the galactic pole as given in:
	   Liu/Zhu/Zhang, „Reconsidering the galactic coordinate system“, A & A No. AA2010, Oct. 2010 */
	{0, 0, false, 0},
	/* 33: Galactic Equator Mula, i.e. galactic/ecliptic intersection point in the middle of lunar mansion Mula */
	{0, 0, false, 0},
	/* 34: Skydram/Galactic Alignment (R. Mardyks); autumn equinox aligned with Galactic Equator/Pole */
	{2451079.734892000, 30, false, 0},
	/* 35: Chandra Hari */
	{0, 0, false, 0},
	/* 36: Dhruva Galactic Centre Middle of Mula (Ernst Wilhelm) */
	{0, 0, false, 0},
	/* 37: Kali 3623 = 522 CE, Ujjain (75.7684565), based on Kali midnight and year length of Suryasiddhanta */
	{1911797.740782065, 0, true, 0},
	/* 38: Babylonian (Britton 2010)  John P. Britton, "Studies in Babylonian lunar theory: part III. The introduction
	of the uniform zodiac", in Arch. Hist. Exact. Sci. (2010)64:617-663, p. 630. */
	{1721057.5, -3.2, true, -1},
	/* 39: Sunil Sheoran ("Vedic") S. Sheoran, "The Science of Time and Timeline of World History", 2017. */
	{0, 0, false, 0},
	/* 40: Galactic Center at 0 Capricon (Cochrane) */
	{0, 0, false, 0},
	/* 41: "Galactic Equatorial" (N.A. Fiorenza) */
	{2451544.5, 25.0, true, 0},
	/* 42: Vettius Valens (Moon; derived from Holden 1995 p. 12 for epoch of Valens 1 Jan. 150 CE julian) */
	{1775845.5, -2.9422, true, -1},
	/* 43: Lahiri (1940), book "Panchanga darpan": 22°26'45".50 + 50".25748T + 0".00011115T^2 */
	{J1900, 22.44597222, false, SEMOD_PREC_NEWCOMB},
	/* 44: Lahiri (VP285), mean sun at 360° in 285CE; epoch for mean sun at 0 acc. to Simon 1994, corrected for Vondrak
	precession (Preface to Lahiri's "Indian Ephemeris" 1980) */
	{1825235.2458513028, 0.0, false, 0},
	/* 45: Krishnamurti from mean equinox 291, based on Newcomb precession, according to D. Senthilathiban, "Study of KP
	Ayanamsa with Modern Precession Theories" (2019), but using precession Vondrak 2011 and correction base on Newcomb
	precession. */
	{1827424.752255678, 0.0, false, 0},
	/* 46: Lahiri original: Calendar Reform Committee 1956, before the correction by 0.658" in IAE 1985.
	   The subtracted value is nutation according to Woolard 1953. However, nutation Woolard was used by IENA/IAE only
	   from 1960 on, so this value is not correct. In order to reproduce mean ayanamshas of IENA >=1960, we could choose
	   23.25 - 0.00464207 + 0.07 / 3600.0 as initial value in 1956. However this will not help to reproduce true
	   ayanamshas. A deviation of around 0.1" remains, for unknown reasons. The difference between Lahiri (1) and
	   Lahiri ICRC (45) amounts to 1.1". */
	{2435553.5, 23.25 - 0.00464207, false, SEMOD_PREC_NEWCOMB}, // 46: SE_SIDM_LAHIRI_ICRC
	/*************************/
}

//// Epsilon represents the obliquity of ecliptic
//type Epsilon struct {
//	Teps float64 // jd (Julian date)
//	Eps  float64 // eps (epsilon/obliquity value)
//	Seps float64 // sin(eps)
//	Ceps float64 // cos(eps)
//}

// PlanetData holds orbital and computational data for a celestial body
type PlanetData struct {
	IBdy int   // Internal body number
	IFlg int32 /* contains several bit flags describing the data:
	 * SEI_FLG_HELIO: true if helio, false if bary
	 * SEI_FLG_ROTATE: TRUE if coefficients are referred to coordinate system of orbital plane
	 * SEI_FLG_ELLIPSE: TRUE if reference ellipse */
	NCoe int /* # of coefficients of ephemeris polynomial, is polynomial order + 1  */
	/* where is the segment index on the file */
	Lndx0   int32   /* file position of begin of planet's index */
	Nndx    int32   /* number of index entries on file: computed */
	TfStart float64 /* file contains ephemeris for tfstart thru tfend */
	TfEnd   float64 /*      for this particular planet !!!            */
	DSeg    float64 /* segment size (days covered by a polynomial)  */
	/* orbital elements: */
	TElem float64 /* epoch of elements */
	Prot  float64
	Qrot  float64
	DProt float64
	DQrot float64
	Rmax  float64 /* normalisation factor of cheby coefficients */
	/* in addition, if reference ellipse is used: */
	Peri  float64
	DPeri float64
	RefEp []float64 /* pointer to cheby coeffs of reference ellipse, size of data is 2 x ncoe */
	/* unpacked segment information, only updated when a segment is read: */
	TSeg0   float64     // Start JD of current segment
	TSeg1   float64     // End JD of current segment
	SegP    []float64   /* pointer to unpacked cheby coeffs of segment; the size is 3 x ncoe */
	NEval   int         /* how many coefficients to evaluate. this may be less than ncoe */
	TEval   float64     /* time for which previous computation was made */
	IEphe   int32       /* which ephemeris was used */
	X       [6]float64  /* position and speed vectors equatorial J2000 */
	XFlgs   int32       /* hel., light-time, aberr., prec. flags etc. */
	XReturn [24]float64 /* return positions:
	 * xreturn+0	ecliptic polar coordinates
	 * xreturn+6	ecliptic cartesian coordinates
	 * xreturn+12	equatorial polar coordinates
	 * xreturn+18	equatorial cartesian coordinates
	 */
}

// STR represents radians per arc second
const STR = 4.8481368110953599359e-6

//// Nut represents nutation parameters
//type Nut struct {
//	Tnut   float64       // Time parameter for nutation
//	Nutlo  [2]float64    // Nutation in longitude and obliquity
//	Snut   float64       // Sine of nutation in obliquity
//	Cnut   float64       // Cosine of nutation in obliquity
//	Matrix [3][3]float64 // 3x3 matrix
//}

// Plantbl represents planetary table parameters
type Plantbl struct {
	MaxHarmonic [9]byte   // Maximum harmonic values
	MaxPowerOfT byte      // Maximum power of T
	ArgTbl      []int8    // Argument table (signed char -> int8)
	LonTbl      []float64 // Longitude table
	LatTbl      []float64 // Latitude table
	RadTbl      []float64 // Radius table
	Distance    float64   // Distance value
}

//// FileData holds information about ephemeris files
//type FileData struct {
//	Fnam       string
//	FVersion   int
//	AstName    string
//	SwephDenum int32
//	Fptr       *os.File
//	TfStart    float64
//	TfEnd      float64
//	IFlg       int32
//	NPl        int16
//	IPl        [SEI_FILE_NMAXPLAN]int
//}
//
//type GenConst struct {
//	Clight       float64 // Speed of light
//	Aunit        float64 // Astronomical unit
//	Helgravconst float64 // Heliocentric gravitational constant
//	Ratme        float64 // Ratio of mass of Earth
//	Sunradius    float64 // Solar radius
//}
//
//type SavePositions struct {
//	Ipl      int     // Planet index
//	Tsave    float64 // Time of save
//	Iflgsave int32   // Saved flags
//	// Position at t = tsave, in:
//	// - ecliptic polar (offset 0)
//	// - ecliptic cartesian (offset 6)
//	// - equatorial polar (offset 12)
//	// - equatorial cartesian coordinates (offset 18)
//	// 6 doubles each for position and speed coordinates
//	Xsaves [24]float64
//}

type NodeData struct {
	Teval float64    // Time for which last computation was made
	Iephe int32      // Which ephemeris was used
	X     [6]float64 // Position and speed vectors equatorial J2000
	Xflgs int32      // Heliocentric, light-time, aberration, precession flags etc.
	// Return positions:
	// - ecliptic polar coordinates (offset 0)
	// - ecliptic cartesian coordinates (offset 6)
	// - equatorial polar coordinates (offset 12)
	// - equatorial cartesian coordinates (offset 18)
	Xreturn [24]float64
}

//
//type TopoData struct {
//	Geolon float64    // Geographic longitude
//	Geolat float64    // Geographic latitude
//	Geoalt float64    // Geographic altitude
//	Teval  float64    // Evaluation time
//	TjdUt  float64    // Julian Day UT
//	Xobs   [6]float64 // Observer position
//}
//
//type SidData struct {
//	SidMode int32   // Sidereal mode
//	AyanT0  float64 // Ayanamsha at T0
//	T0      float64 // Reference time T0
//	T0IsUT  bool    // True if T0 is UT
//}
//
//const SWI_STAR_LENGTH = 40
//
//// FixedStar represents fixed star data
//type FixedStar struct {
//	Skey      [SWI_STAR_LENGTH + 2]byte // May be prefixed with comma, one char more
//	StarName  [SWI_STAR_LENGTH + 1]byte
//	StarBayer [SWI_STAR_LENGTH + 1]byte
//	StarNo    [10]byte
//	Epoch     float64
//	Ra        float64
//	De        float64
//	RaMot     float64
//	DeMot     float64
//	RadVel    float64
//	Parall    float64
//	Mag       float64
//}

// SWE_DATA_DPSI_DEPS represents the number of days for dpsi and deps data (100 years after 1962)
const SWE_DATA_DPSI_DEPS = 36525

//// Interpol represents interpolation data for nutation
//type Interpol struct {
//	TjdNut0  float64
//	TjdNut2  float64
//	NutDpsi0 float64
//	NutDpsi1 float64
//	NutDpsi2 float64
//	NutDeps0 float64
//	NutDeps1 float64
//	NutDeps2 float64
//}
//
//// SweData represents Swiss Ephemeris data
//type SweData struct {
//	EphePathIsSet      bool
//	JplFileIsOpen      bool
//	FixFp              *os.File //Fixed stars file pointer
//	EphePath           string   // [AS_MAXCH]byte
//	JplFnam            string   // [AS_MAXCH]byte
//	JplDenum           int32
//	LastEpheFlag       int32
//	GeoposIsSet        bool
//	AyanaIsSet         bool
//	IsOldStarfile      bool
//	EopTjdBeg          float64
//	EopTjdBegHorizons  float64
//	EopTjdEnd          float64
//	EopTjdEndAdd       float64
//	EopDpsiLoaded      int
//	TidAcc             float64
//	IsTidAccManual     bool
//	InitDtDone         bool
//	SwedIsInitialised  bool
//	DeltaTUserdefIsSet bool
//	DeltaTUserdef      float64
//	AstG               float64
//	AstH               float64
//	AstDiam            float64
//	Astelem            [AS_MAXCH * 10]byte
//	ISavedPlanetName   int
//	SavedPlanetName    string
//	Dpsi               []float64
//	Deps               []float64
//	Timeout            int32
//	AstroModels        [SEI_NMODELS]int32
//	DoInterpolateNut   bool
//	Interpol           Interpol
//	Fidat              [SEI_NEPHFILES]FileData
//	Gcdat              GenConst
//	//	#if 0
//	//     struct node_data nddat[SEI_NNODE_ETC];
//	//	#else
//	//     struct plan_data nddat[SEI_NNODE_ETC];
//	//	#endif
//	Pldat            [SEI_NPLANETS]PlanetData
//	Nddat            [SEI_NNODE_ETC]PlanetData
//	Savedat          [SE_NPLANETS + 1]SavePositions
//	Oec              Epsilon
//	Oec2000          Epsilon
//	Nut              Nut
//	Nut2000          Nut
//	Nutv             Nut
//	Topd             TopoData
//	Sidd             SidData
//	NFixstarsReal    bool
//	NFixstarsNamed   bool
//	NFixstarsRecords bool
//	FixedStars       *FixedStar
//}
//
//var sweData SweData
