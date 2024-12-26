package internal

import (
	"bufio"
	"errors"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

// constants swephlib.h-0061
const (
	PREC_IAU_1976_CTIES = 2.0 /* J2000 +/- two centuries */
	PREC_IAU_2000_CTIES = 2.0 /* J2000 +/- two centuries */
	// we use P03 for whole ephemeris
	PREC_IAU_2006_CTIES             = 75.0 /* J2000 +/- 75 centuries */
	DPSI_DEPS_IAU1980_TJD0_HORIZONS = 2437684.5
	HORIZONS_TJD0_DPSI_DEPS_IAU1980 = 2437684.5
)

// ===== 0104 ========= swe_degnorm swephlib.c-0104 ==================================================================

// SweDegnorm reduces x modulo 360 degrees
func SweDegnorm(x float64) float64 {
	y := math.Mod(x, 360.0)
	if math.Abs(y) < 1e-13 { // Alois fix 11-dec-1999
		y = 0
	}
	if y < 0.0 {
		y += 360.0
	}
	return y
}

// ===== 0115 ============= swe_radnorm swephlib.c-0115 =============================================================

// SweRadnorm reduces x modulo TWOPI radians
func SweRadnorm(x float64) float64 {
	y := math.Mod(x, TWOPI)
	if math.Abs(y) < 1e-13 { // Alois fix 11-dec-1999
		y = 0
	}
	if y < 0.0 {
		y += TWOPI
	}
	return y
}

// ===== 0139 ===== swi_mod2PI swephlib.c-0139

// Mod2PI reduces x modulo 2π, ensuring the result is in the range [0, 2π)
func Mod2PI(x float64) float64 {
	y := math.Mod(x, TWOPI)
	if y < 0.0 {
		y += TWOPI
	}
	return y
}

// ===== 0187 =================== swi_edcheb swephlib.c-0187 =========================================================

// swiEdcheb evaluates the derivative of a Chebyshev series.
// Parameters:
//
//	x: the point at which to evaluate the derivative
//	coef: slice containing the Chebyshev coefficients
//	ncf: number of coefficients
//
// Returns: the derivative value at point x
func swiEdcheb(x float64, coef []float64, ncf int) float64 {
	x2 := x * 2.0
	var (
		bjpl float64 // b_{j+1}
		xjpl float64 // x_{j+1}
		bf   float64 // b_{final}
		bj   float64 // b_j
		xj   float64 // x_j
		bjp2 float64 // b_{j+2}
		xjp2 float64 // x_{j+2}
	)
	// Evaluate derivative using Clenshaw's recurrence formula
	for j := ncf - 1; j >= 1; j-- {
		dj := float64(j + j)
		xj = coef[j]*dj + xjp2
		bj = x2*bjpl - bjp2 + xj
		bf = bjp2
		bjp2 = bjpl
		bjpl = bj
		xjp2 = xjpl
		xjpl = xj
	}
	return (bj - bf) * 0.5
}

// ===== 0293 ===== swi_coortrf2

// swiCoortrf handles the conversion between ecliptical and equatorial cartesian coordinates
// for ecl. to equ.  eps must be negative
// for equ. to ecl.  eps must be positive
func swiCoortrf(xpo []float64, eps float64) []float64 {
	// Create new slice for output instead of modifying in place
	xpn := make([]float64, 3)
	sineps := math.Sin(eps)
	coseps := math.Cos(eps)
	x := make([]float64, 3)
	x[0] = xpo[0]
	x[1] = xpo[1]*coseps + xpo[2]*sineps
	x[2] = -xpo[1]*sineps + xpo[2]*coseps
	xpn[0] = x[0]
	xpn[1] = x[1]
	xpn[2] = x[2]
	return xpn
}

// ===== 0310 ===== swei_cartpol swephlib.c-0310 =====================================================================

// swiCartpol converts cartesian (x[3]) to polar coordinates (l[3]).
// x = l is allowed.
// if |x| = 0, then lon, lat and rad := 0.
func swiCartpol(x []float64) []float64 {
	l := make([]float64, 3)
	if x[0] == 0 && x[1] == 0 && x[2] == 0 {
		return l
	}
	ll := make([]float64, 3)
	rxy := x[0]*x[0] + x[1]*x[1]
	ll[2] = math.Sqrt(rxy + x[2]*x[2])
	rxy = math.Sqrt(rxy)
	ll[0] = math.Atan2(x[1], x[0])
	if ll[0] < 0.0 {
		ll[0] += TWOPI
	}
	if rxy == 0 {
		if x[2] >= 0 {
			ll[1] = PI / 2
		} else {
			ll[1] = -(PI / 2)
		}
	} else {
		ll[1] = math.Atan(x[2] / rxy)
	}
	l[0] = ll[0]
	l[1] = ll[1]
	l[2] = ll[2]
	return l
}

// ===== 0340 ========================================================================================================

// swiPolcart converts polar (l[3]) to cartesian coordinates (x[3]).
// x = l is allowed.
func swiPolcart(l []float64) []float64 {
	x := make([]float64, 3)
	xx := make([]float64, 3)
	cosl1 := math.Cos(l[1])
	xx[0] = l[2] * cosl1 * math.Cos(l[0])
	xx[1] = l[2] * cosl1 * math.Sin(l[0])
	xx[2] = l[2] * math.Sin(l[1])
	x[0] = xx[0]
	x[1] = xx[1]
	x[2] = xx[2]
	return x
}

// ===== 0467 ==================== constants for swiLdPeps swephib.c-0467 ============================================
// functions for precession and ecliptic obliquity according to Vondrák et alii, 2011
const (
	AS2R      = DEGTORAD / 3600.0
	D2PI      = 2 * math.Pi
	EPS0      = 84381.406 * AS2R
	NPOL_PEPS = 4
	NPER_PEPS = 10
	NPOL_PECL = 4
	NPER_PECL = 8
	NPOL_PEQU = 4
	NPER_PEQU = 14
)

// For pre_peps:
// poynomials
var pepol = [NPOL_PEPS][2]float64{
	{+8134.017132, +84028.206305},
	{+5043.0520035, +0.3624445},
	{-0.00710733, -0.00004039},
	{+0.000000271, -0.000000110},
}

// periodics
var peper = [5][NPER_PEPS]float64{
	{+409.90, +396.15, +537.22, +402.90, +417.15, +288.92, +4043.00, +306.00, +277.00, +203.00},
	{-6908.287473, -3198.706291, +1453.674527, -857.748557, +1173.231614, -156.981465, +371.836550, -216.619040, +193.691479, +11.891524},
	{+753.872780, -247.805823, +379.471484, -53.880558, -90.109153, -353.600190, -63.115353, -28.248187, +17.703387, +38.911307},
	{-2845.175469, +449.844989, -1255.915323, +886.736783, +418.887514, +997.912441, -240.979710, +76.541307, -36.788069, -170.964086},
	{-1704.720302, -862.308358, +447.832178, -889.571909, +190.402846, -56.564991, -296.222622, -75.859952, +67.473503, +3.014055},
}

// For pre_pecl:
// polynomials
var pqpol = [NPOL_PECL][2]float64{
	{+5851.607687, -1600.886300},
	{-0.1189000, +1.1689818},
	{-0.00028913, -0.00000020},
	{+0.000000101, -0.000000437},
}

// periodics
var pqper = [5][NPER_PECL]float64{
	{708.15, 2309, 1620, 492.2, 1183, 622, 882, 547},
	{-5486.751211, -17.127623, -617.517403, 413.44294, 78.614193, -180.732815, -87.676083, 46.140315},
	// typo fixed according to A&A 541, C1 (2012)
	{-684.66156, 2446.28388, 399.671049, -356.652376, -186.387003, -316.80007, 198.296701, 101.135679},
	{667.66673, -2354.886252, -428.152441, 376.202861, 184.778874, 335.321713, -185.138669, -120.97283},
	{-5523.863691, -549.74745, -310.998056, 421.535876, -36.776172, -145.278396, -34.74445, 22.885731},
}

// For pre_pequ
// polynomials
var xypol = [NPOL_PEQU][2]float64{
	{+5453.282155, -73750.930350},
	{+0.4252841, -0.7675452},
	{-0.00037173, -0.00018725},
	{-0.000000152, +0.000000231},
}

// periodics
var xyper = [5][NPER_PEQU]float64{
	{256.75, 708.15, 274.2, 241.45, 2309, 492.2, 396.1, 288.9, 231.1, 1610, 620, 157.87, 220.3, 1200},
	{-819.940624, -8444.676815, 2600.009459, 2755.17563, -167.659835, 871.855056, 44.769698, -512.313065, -819.415595, -538.071099, -189.793622, -402.922932, 179.516345, -9.814756},
	{75004.344875, 624.033993, 1251.136893, -1102.212834, -2660.66498, 699.291817, 153.16722, -950.865637, 499.754645, -145.18821, 558.116553, -23.923029, -165.405086, 9.344131},
	{81491.287984, 787.163481, 1251.296102, -1257.950837, -2966.79973, 639.744522, 131.600209, -445.040117, 584.522874, -89.756563, 524.42963, -13.549067, -210.157124, -44.919798},
	{1558.515853, 7774.939698, -2219.534038, -2523.969396, 247.850422, -846.485643, -1393.124055, 368.526116, 749.045012, 444.704518, 235.934465, 374.049623, -171.33018, -22.899655},
}

// ===== 0535 ============= swi_ldp_peps swephlib.c-0535  see constants above =======================================

func swiLdpPeps(tjd float64) (dpre, deps float64) {
	t := (tjd - J2000) / 36525.0
	p := 0.0
	q := 0.0
	// Periodic terms
	for i := 0; i < NPER_PEPS; i++ {
		w := D2PI * t
		a := w / peper[0][i]
		s := math.Sin(a)
		c := math.Cos(a)
		p += c*peper[1][i] + s*peper[3][i]
		q += c*peper[2][i] + s*peper[4][i]
	}
	// Polynomial terms
	w := 1.0
	for i := 0; i < NPOL_PEPS; i++ {
		p += pepol[i][0] * w
		q += pepol[i][1] * w
		w *= t
	}
	// Both to radians
	p *= AS2R
	q *= AS2R
	return p, q
}

// ===== 0571 ===== pre_pecl =========================================================================================

// prePecl calclates the precession of the ecliptic
func prePecl(tjd float64) []float64 {
	vec := make([]float64, 3)
	t := (tjd - J2000) / 36525.0
	p := 0.0
	q := 0.0
	// Periodic terms
	for i := 0; i < NPER_PECL; i++ {
		w := D2PI * t
		a := w / pqper[0][i]
		s := math.Sin(a)
		c := math.Cos(a)
		p += c*pqper[1][i] + s*pqper[3][i]
		q += c*pqper[2][i] + s*pqper[4][i]
	}
	// Polynomial terms
	w := 1.0
	for i := 0; i < NPOL_PECL; i++ {
		p += pqpol[i][0] * w
		q += pqpol[i][1] * w
		w *= t
	}
	// Both to radians
	p *= AS2R
	q *= AS2R
	// Ecliptic pole vector
	z := 1 - p*p - q*q
	if z < 0 {
		z = 0
	} else {
		z = math.Sqrt(z)
	}
	s := math.Sin(EPS0)
	c := math.Cos(EPS0)
	vec[0] = p
	vec[1] = -q*c - z*s
	vec[2] = -q*s + z*c
	return vec
}

// ===== 0618 ===== pre_pequ =========================================================================================

// prePequ calculates the precession of the equator
func prePequ(tjd float64) []float64 {
	veq := make([]float64, 3)
	t := (tjd - J2000) / 36525.0
	x := 0.0
	y := 0.0
	// Periodic terms
	for i := 0; i < NPER_PEQU; i++ {
		w := D2PI * t
		a := w / xyper[0][i]
		s := math.Sin(a)
		c := math.Cos(a)
		x += c*xyper[1][i] + s*xyper[3][i]
		y += c*xyper[2][i] + s*xyper[4][i]
	}
	// Polynomial terms
	w := 1.0
	for i := 0; i < NPOL_PEQU; i++ {
		x += xypol[i][0] * w
		y += xypol[i][1] * w
		w *= t
	}
	x *= AS2R
	y *= AS2R
	// Equator pole vector
	veq[0] = x
	veq[1] = y
	w = x*x + y*y
	if w < 1 {
		veq[2] = math.Sqrt(1 - w)
	} else {
		veq[2] = 0
	}
	return veq
}

// ===== 0697 ============ constants for get_owen_t0_icof  swephlib.c-0697 ===========================================
// precession according to Owen 1990: Owen, William M., Jr., (JPL) "A Theory of the Earth's Precession Relative to the
// Invariable Plane of the Solar System", Ph.D. Dissertation, University of Florida, 1990.
// Implemented for time range -18000 to 14000.
//
// p. 177: central time Tc = -160, covering time span -200 <= T <= -120 i.e. -14000 +- 40 centuries
// p. 178: central time Tc = -80, covering time span -120 <= T <= -40 i.e. -6000 +- 40 centuries
// p. 179: central time Tc = 0, covering time span -40 <= T <= +40 i.e. 2000 +- 40 centuries
// p. 180: central time Tc = 80, covering time span 40 <= T <= 120 i.e. 10000 +- 40 centuries
// p. 181: central time Tc = 160, covering time span 120 <= T <= 200 i.e. 10000 +- 40 centuries
var owenEps0Coef = [5][10]float64{
	{23.699391439256386, 5.2330816033981775e-1, -5.6259493384864815e-2, -8.2033318431602032e-3, 6.6774163554156385e-4, 2.4931584012812606e-5, -3.1313623302407878e-6, 2.0343814827951515e-7, 2.9182026615852936e-8, -4.1118760893281951e-9},
	{24.124759551704588, -1.2094875596566286e-1, -8.3914869653015218e-2, 3.5357075322387405e-3, 6.4557467824807032e-4, -2.5092064378707704e-5, -1.7631607274450848e-6, 1.3363622791424094e-7, 1.5577817511054047e-8, -2.4613907093017122e-9},
	{23.439103144206208, -4.9386077073143590e-1, -2.3965445283267805e-4, 8.6637485629656489e-3, -5.2828151901367600e-5, -4.3951004595359217e-5, -1.1058785949914705e-6, 6.2431490022621172e-8, 3.4725376218710764e-8, 1.3658853127005757e-9},
	{22.724671295125046, -1.6041813558650337e-1, 7.0646783888132504e-2, 1.4967806745062837e-3, -6.6857270989190734e-4, 5.7578378071604775e-6, 3.3738508454638728e-6, -2.2917813537654764e-7, -2.1019907929218137e-8, 4.3139832091694682e-9},
	{22.914636050333696, 3.2123508304962416e-1, 3.6633220173792710e-2, -5.9228324767696043e-3, -1.882379107379328e-4, 3.2274552870236244e-5, 4.9052463646336507e-7, -5.9064298731578425e-8, -2.0485712675098837e-8, -6.2163304813908160e-10},
}

var owenPsiaCoef = [5][10]float64{
	{-218.57864954903122, 51.752257487741612, 1.3304715765661958e-1, 9.2048123521890745e-2, -6.0877528127241278e-3, -7.0013893644531700e-5, -4.9217728385458495e-5, -1.8578234189053723e-6, 7.4396426162029877e-7, -5.9157528981843864e-9},
	{-111.94350527506128, 55.175558131675861, 4.7366115762797613e-1, -4.7701750975398538e-2, -9.2445765329325809e-3, 7.0962838707454917e-4, 1.5140455277814658e-4, -7.7813159018954928e-7, -2.4729402281953378e-6, -1.0898887008726418e-7},
	{-2.041452011529441e-1, 55.969995858494106, -1.9295093699770936e-1, -5.6819574830421158e-3, 1.1073687302518981e-2, -9.0868489896815619e-5, -1.1999773777895820e-4, 9.9748697306154409e-6, 5.7911493603430550e-7, -2.3647526839778175e-7},
	{111.61366860604471, 56.404525305162447, 4.4403302410703782e-1, 7.1490030578883907e-2, -4.9184559079790816e-3, -1.3912698949042046e-3, -6.8490613661884005e-5, 1.2394328562905297e-6, 1.7719847841480384e-6, 2.4889095220628068e-7},
	{228.40683531269390, 60.056143904919826, 2.9583200718478960e-2, -1.5710838319490748e-1, -7.0017356811600801e-3, 3.3009615142224537e-3, 2.0318123852537664e-4, -6.5840216067828310e-5, -5.9077673352976155e-6, 1.3983942185303064e-6},
}

var owenOmaCoef = [5][10]float64{
	{25.541291140949806, 2.377889511272162e-1, -3.7337334723142133e-1, 2.4579295485161534e-2, 4.3840999514263623e-3, -3.1126873333599556e-4, -9.8443045771748915e-6, -7.9403103080496923e-7, 1.0840116743893556e-9, 9.2865105216887919e-9},
	{24.429357654237926, -9.5205745947740161e-1, 8.6738296270534816e-2, 3.0061543426062955e-2, -4.1532480523019988e-3, -3.7920928393860939e-4, 3.5117012399609737e-5, 4.6811877283079217e-6, -8.1836046585546861e-8, -6.1803706664211173e-8},
	{23.450465062489337, -9.7259278279739817e-2, 1.1082286925130981e-2, -3.1469883339372219e-2, -1.0041906996819648e-4, 5.6455168475133958e-4, -8.4403910211030209e-6, -3.8269157371098435e-6, 3.1422585261198437e-7, 9.3481729116773404e-9},
	{22.581778052947806, -8.7069701538602037e-1, -9.8140710050197307e-2, 2.6025931340678079e-2, 4.8165322168786755e-3, -1.906558772193363e-4, -4.6838759635421777e-5, -1.6608525315998471e-6, -3.2347811293516124e-8, 2.8104728109642000e-9},
	{21.518861835737142, 2.0494789509441385e-1, 3.5193604846503161e-1, 1.5305977982348925e-2, -7.5015367726336455e-3, -4.0322553186065610e-4, 1.0655320434844041e-4, 7.1792339586935752e-6, -1.603874697543020e-6, -1.613563462813512e-7},
}

var owenChiaCoef = [5][10]float64{
	{8.2378850337329404e-1, -3.7443109739678667, 4.0143936898854026e-1, 8.1822830214590811e-2, -8.5978790792656293e-3, -2.8350488448426132e-5, -4.2474671728156727e-5, -1.6214840884656678e-6, 7.8560442001953050e-7, -1.032016641696707e-8},
	{-2.1726062070318606, 7.8470515033132925e-1, 4.4044931004195718e-1, -8.0671247169971653e-2, -8.9672662444325007e-3, 9.2248978383109719e-4, 1.5143472266372874e-4, -1.6387009056475679e-6, -2.4405558979328144e-6, -1.0148113464009015e-7},
	{-4.8518673570735556e-1, 1.0016737299946743e-1, -4.7074888613099918e-1, -5.8604054305076092e-3, 1.4300208240553435e-2, -6.7127991650300028e-5, -1.3703764889645475e-4, 9.0505213684444634e-6, 6.0368690647808607e-7, -2.2135404747652171e-7},
	{-2.0950740076326087, -9.4447359463206877e-1, 4.0940512860493755e-1, 1.0261699700263508e-1, -5.3133241571955160e-3, -1.6634631550720911e-3, -5.9477519536647907e-5, 2.9651387319208926e-6, 1.6434499452070584e-6, 2.3720647656961084e-7},
	{6.3315163285678715e-1, 3.5241082918420464, 2.1223076605364606e-1, -1.5648122502767368e-1, -9.1964075390801980e-3, 3.3896161239812411e-3, 2.1485178626085787e-4, -6.6261759864793735e-5, -5.9257969712852667e-6, 1.3918759086160525e-6},
}

// ===== 0747 ======== get_owen_t0_icof swephib.c-0747, see preceding constants ======================================

func getOwenT0Icof(tjd float64) (t0 float64, icof int) {
	t0s := [5]float64{-3392455.5, -470455.5, 2451544.5, 5373544.5, 8295544.5}
	t0 = t0s[0]
	j := 0
	for i := 1; i < 5; i++ {
		if tjd >= (t0s[i-1]+t0s[i])/2 {
			t0 = t0s[i]
			j++
		}
	}
	icof = j
	return t0, icof
}

// ===== 0971 ===== precess_1 swephlib.c-0971 ========================================================================

// Precession of the equinox and ecliptic from epoch Julian date J to or from J2000.0
// Original program by Steve Moshier. Changes in program structure and implementation of IAU 2003 (P03) and Vondrak 2011
// by Dieter Koch.
// SEMOD_PREC_VONDRAK_2011 J. Vondrák, N. Capitaine, and P. Wallace, "New precession expressions, valid for long time
// intervals", A&A 534, A22 (2011)
// SEMOD_PREC_IAU_2006 N. Capitaine, P.T. Wallace, and J. Chapront, "Expressions for IAU 2000 precession quantities",
// 2003, A&A 412, 567-586 (2003). This is a "short" term model, that can be combined with other models
// SEMOD_PREC_WILLIAMS_1994 James G. Williams, "Contributions to the Earth's obliquity rate, precession, and nutation,"
// Astron. J. 108, 711-724 (1994).
// SEMOD_PREC_SIMON_1994 J. L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze', G. Francou, and J. Laskar,
// "Numerical Expressions for precession formulae and mean elements for the Moon and the planets," Astronomy and
// Astrophysics 282, 663-683 (1994).
// SEMOD_PREC_IAU_1976 IAU Coefficients are from: J. H. Lieske, T. Lederle, W. Fricke, and B. Morando, "Expressions for
// the Precession Quantities Based upon the IAU (1976) System of Astronomical Constants,"  Astronomy and Astrophysics
// 58, 1-16 (1977). This is a "short" term model, that can be combined with other models
// SEMOD_PREC_LASKAR_1986 Newer formulas that cover a much longer time span are from: J. Laskar, "Secular terms of
// classical planetary theories using the results of general theory," Astronomy and Astrophysics 157, 59070 (1986).
// See also:
// P. Bretagnon and G. Francou, "Planetary theories in rectangular and spherical variables. VSOP87 solutions,"
// Astronomy and Astrophysics 202, 309-315 (1988).
// Bretagnon and Francou's expansions for the node and inclination of the ecliptic were derived from Laskar's data but
// were truncated after the term in T**6. I have recomputed these expansions from Laskar's data, retaining powers up to
// T**10 in the result.

func precess1(R []float64, J float64, direction int, precMethod int) int {
	var T, Z, z, TH float64
	var x [3]float64

	if J == J2000 {
		return 0
	}

	T = (J - J2000) / 36525.0

	switch precMethod {
	case SEMOD_PREC_IAU_1976:
		Z = ((0.017998*T+0.30188)*T + 2306.2181) * T * DEGTORAD / 3600
		z = ((0.018203*T+1.09468)*T + 2306.2181) * T * DEGTORAD / 3600
		TH = ((-0.041833*T-0.42665)*T + 2004.3109) * T * DEGTORAD / 3600

	case SEMOD_PREC_IAU_2000:
		Z = (((((-0.0000002*T-0.0000327)*T+0.0179663)*T+0.3019015)*T+2306.0809506)*T + 2.5976176) * DEGTORAD / 3600
		z = (((((-0.0000003*T-0.000047)*T+0.0182237)*T+1.0947790)*T+2306.0803226)*T - 2.5976176) * DEGTORAD / 3600
		TH = ((((-0.0000001*T-0.0000601)*T-0.0418251)*T-0.4269353)*T + 2004.1917476) * T * DEGTORAD / 3600

	case SEMOD_PREC_IAU_2006:
		T = (J - J2000) / 36525.0
		Z = (((((-0.0000003173*T-0.000005971)*T+0.01801828)*T+0.2988499)*T+2306.083227)*T + 2.650545) * DEGTORAD / 3600
		z = (((((-0.0000002904*T-0.000028596)*T+0.01826837)*T+1.0927348)*T+2306.077181)*T - 2.650545) * DEGTORAD / 3600
		TH = ((((-0.00000011274*T-0.000007089)*T-0.04182264)*T-0.4294934)*T + 2004.191903) * T * DEGTORAD / 3600

	case SEMOD_PREC_BRETAGNON_2003:
		Z = ((((((-0.00000000013*T-0.0000003040)*T-0.000005708)*T+0.01801752)*T+0.3023262)*T+2306.080472)*T + 2.72767) * DEGTORAD / 3600
		z = ((((((-0.00000000005*T-0.0000002486)*T-0.000028276)*T+0.01826676)*T+1.0956768)*T+2306.076070)*T - 2.72767) * DEGTORAD / 3600
		TH = ((((((0.000000000009*T+0.00000000036)*T-0.0000001127)*T-0.000007291)*T-0.04182364)*T-0.4266980)*T + 2004.190936) * T * DEGTORAD / 3600

	case SEMOD_PREC_NEWCOMB:
		mills := 365242.198782 // trop. millennia
		t1 := (J2000 - B1850) / mills
		t2 := (J - B1850) / mills
		T = t2 - t1
		T2 := T * T
		T3 := T2 * T
		Z1 := 23035.5548 + 139.720*t1 + 0.069*t1*t1
		Z = Z1*T + (30.242-0.269*t1)*T2 + 17.996*T3
		z = Z1*T + (109.478-0.387*t1)*T2 + 18.324*T3
		TH = (20051.125-85.294*t1-0.365*t1*t1)*T + (-42.647-0.365*t1)*T2 - 41.802*T3
		Z *= (DEGTORAD / 3600.0)
		z *= (DEGTORAD / 3600.0)
		TH *= (DEGTORAD / 3600.0)

	default:
		return 0
	}

	sinth := math.Sin(TH)
	costh := math.Cos(TH)
	sinZ := math.Sin(Z)
	cosZ := math.Cos(Z)
	sinz := math.Sin(z)
	cosz := math.Cos(z)
	A := cosZ * costh
	B := sinZ * costh

	if direction < 0 { // From J2000.0 to J
		x[0] = (A*cosz-sinZ*sinz)*R[0] - (B*cosz+cosZ*sinz)*R[1] - sinth*cosz*R[2]
		x[1] = (A*sinz+sinZ*cosz)*R[0] - (B*sinz-cosZ*cosz)*R[1] - sinth*sinz*R[2]
		x[2] = cosZ*sinth*R[0] - sinZ*sinth*R[1] + costh*R[2]
	} else { // From J to J2000.0
		x[0] = (A*cosz-sinZ*sinz)*R[0] + (A*sinz+sinZ*cosz)*R[1] + cosZ*sinth*R[2]
		x[1] = -(B*cosz+cosZ*sinz)*R[0] - (B*sinz-cosZ*cosz)*R[1] - sinZ*sinth*R[2]
		x[2] = -sinth*cosz*R[0] - sinth*sinz*R[1] + costh*R[2]
	}
	for i := 0; i < 3; i++ {
		R[i] = x[i]
	}
	return 0
}

// ===== 0828 ==================== epsiln_owen_1986 swephlib.c-0828 ==================================================

func epsilnOwen1986(tjd float64) float64 {
	k := make([]float64, 10)
	tau := make([]float64, 10)
	t0, icof := getOwenT0Icof(tjd)
	eps := 0.0
	tau[0] = 0
	tau[1] = (tjd - t0) / 36525.0 / 40.0
	for i := 2; i <= 9; i++ {
		tau[i] = tau[1] * tau[i-1]
	}
	k[0] = 1
	k[1] = tau[1]
	k[2] = 2*tau[2] - 1
	k[3] = 4*tau[3] - 3*tau[1]
	k[4] = 8*tau[4] - 8*tau[2] + 1
	k[5] = 16*tau[5] - 20*tau[3] + 5*tau[1]
	k[6] = 32*tau[6] - 48*tau[4] + 18*tau[2] - 1
	k[7] = 64*tau[7] - 112*tau[5] + 56*tau[3] - 7*tau[1]
	k[8] = 128*tau[8] - 256*tau[6] + 160*tau[4] - 32*tau[2] + 1
	k[9] = 256*tau[9] - 576*tau[7] + 432*tau[5] - 120*tau[3] + 9*tau[1]
	for i := 0; i < 10; i++ {
		eps += k[i] * owenEps0Coef[icof][i]
	}
	return eps
}

// ===== 0856 ===== data for swi_epsiln swephlib.c-856 ===============================================================

// Obliquity of the ecliptic at Julian date J
// IAU Coefficients are from:
// J. H. Lieske, T. Lederle, W. Fricke, and B. Morando, "Expressions for the Precession Quantities Based upon the IAU
// (1976) System of Astronomical Constants,"  Astronomy and Astrophysics 58, 1-16 (1977).
// Before or after 200 years from J2000, the formula used is from:
// J. Laskar, "Secular terms of classical planetary theories using the results of general theory," Astronomy and
// Astrophysics 157, 59070 (1986).
// Bretagnon, P. et al.: 2003, "Expressions for Precession Consistent with the IAU 2000A Model". A&A 400,785
// B03  	84381.4088  	-46.836051*t  	-1667*10-7*t2  	+199911*10-8*t3  	-523*10-9*t4  	-248*10-10*t5  	-3*10-11*t6
// C03   84381.406  	-46.836769*t  	-1831*10-7*t2  	+20034*10-7*t3  	-576*10-9*t4  	-434*10-10*t5
// See precess and page B18 of the Astronomical Almanac.

// constants for calc_nutation_iau1980, swephlib.c-1432
// Nutation in longitude and obliquity computed at Julian date J.
// References:
// "Summary of 1980 IAU Theory of Nutation (Final Report of the IAU Working Group on Nutation)", P. K. Seidelmann et
// al., in Transactions of the IAU Vol. XVIII A, Reports on Astronomy, P. A. Wayman, ed.; D. Reidel Pub. Co., 1982.
// "Nutation and the Earth's Rotation", I.A.U. Symposium No. 78, May, 1977, page 256. I.A.U., 1980.
// Woolard, E.W., "A redevelopment of the theory of nutation", The Astronomical Journal, 58, 1-3 (1953).
// This program implements all of the 1980 IAU nutation series.
// Results checked at 100 points against the 1986 AA; all agreed.
// - S. L. Moshier, November 1987
//   October, 1992 - typo fixed in nutation matrix
// - D. Koch, November 1995: small changes in structure,
//   Corrections to IAU 1980 Series added from Expl. Suppl. p. 116
// Each term in the expansion has a trigonometric argument given by W = i*MM + j*MS + k*FF + l*DD + m*OM
// where the variables are defined below.
// The nutation in longitude is a sum of terms of the form (a + bT) * sin(W). The terms for nutation in obliquity
// are of the form (c + dT) * cos(W).  The coefficients are arranged in the tabulation as follows:
// Coefficient:
// i  j  k  l  m      a      b      c     d
// 0, 0, 0, 0, 1, -171996, -1742, 92025, 89,
// The first line of the table, above, is done separately since two of the values do not fit into 16 bit integers.
// The values a and c are arc seconds times 10000.  b and d are arc seconds per Julian century times 100000.
// i through m are integers.  See the program for interpretation of MM, MS, etc., which are mean orbital elements of the
// Sun and Moon.
// If terms with coefficient less than X are omitted, the peak errors will be:
//   omit	error,		  omit	error,
//   a <	longitude	  c <	obliquity
// .0005"	.0100"		.0008"	.0094"
// .0046	.0492		.0095	.0481
// .0123	.0880		.0224	.0905
// .0386	.1808		.0895	.1129

const (
	OFFSET_EPS_JPLHORIZONS = 35.95
	DCOR_EPS_JPL_TJD0      = 2437846.5
	NDCOR_EPS_JPL          = 51
)

var dcorEpsJPL = [...]float64{
	36.726, 36.627, 36.595, 36.578, 36.640, 36.659, 36.731, 36.765,
	36.662, 36.555, 36.335, 36.321, 36.354, 36.227, 36.289, 36.348, 36.257, 36.163,
	35.979, 35.896, 35.842, 35.825, 35.912, 35.950, 36.093, 36.191, 36.009, 35.943,
	35.875, 35.771, 35.788, 35.753, 35.822, 35.866, 35.771, 35.732, 35.543, 35.498,
	35.449, 35.409, 35.497, 35.556, 35.672, 35.760, 35.596, 35.565, 35.510, 35.394,
	35.385, 35.375, 35.415,
}

// ===== 0356 ===== swi_cartpol_sp swephlib.c-0356 ===================================================================

// conversion of position and speed. from cartesian (x[6]) to polar coordinates (l[6]).
// x = l is allowed.
// if position is 0, function returns direction of motion.

func swiCartpolSp(x, l []float64) {
	xx := make([]float64, 6)
	ll := make([]float64, 6)

	// zero position
	if x[0] == 0 && x[1] == 0 && x[2] == 0 {
		ll[0], ll[1], ll[3], ll[4] = 0, 0, 0, 0
		ll[5] = math.Sqrt(SquareSum(x[3:]))
		ll = swiCartpol(x[3:]) // TODO Port: check this call of swiCartpol
		ll[2] = 0
		copy(l, ll)
		return
	}

	// zero speed
	if x[3] == 0 && x[4] == 0 && x[5] == 0 {
		l[3], l[4], l[5] = 0, 0, 0
		l = swiCartpol(x) // TODO Port: check this call of swiCartpol
		return
	}

	// position
	rxy := x[0]*x[0] + x[1]*x[1]
	ll[2] = math.Sqrt(rxy + x[2]*x[2])
	rxy = math.Sqrt(rxy)
	ll[0] = math.Atan2(x[1], x[0])
	if ll[0] < 0.0 {
		ll[0] += TWOPI
	}
	ll[1] = math.Atan(x[2] / rxy)
	// speed:
	// 1. rotate coordinate system by longitude of position about z-axis, so that new x-axis = position radius projected
	//    onto x-y-plane in the new coordinate system vy'/r = dlong/dt, where r = sqrt(x^2 +y^2).
	// 2. rotate coordinate system by latitude about new y-axis. vz"/r = dlat/dt, where r = position radius. vx" = dr/dt
	coslon := x[0] / rxy   // cos(l[0])
	sinlon := x[1] / rxy   // sin(l[0])
	coslat := rxy / ll[2]  // cos(l[1])
	sinlat := x[2] / ll[2] // sin(ll[1])
	xx[3] = x[3]*coslon + x[4]*sinlon
	xx[4] = -x[3]*sinlon + x[4]*coslon
	l[3] = xx[4] / rxy // speed in longitude
	xx[4] = -sinlat*xx[3] + coslat*x[5]
	xx[5] = coslat*xx[3] + sinlat*x[5]
	l[4] = xx[4] / ll[2] // speed in latitude
	l[5] = xx[5]         // speed in radius
	// return position
	l[0] = ll[0]
	l[1] = ll[1]
	l[2] = ll[2]
}

//===== 0763 ===== owen_pre_matrix swephlib.c-0763 ==================================================================

// owenPreMatrix calculates the precession matrix using Owen 1990 method
func owenPreMatrix(tjd float64, rp []float64, iflag int32) {
	var eps0, chia, psia, oma float64
	var coseps0, sineps0, coschia, sinchia, cospsia, sinpsia, cosoma, sinoma float64
	k := make([]float64, 10)
	tau := make([]float64, 10)
	t0, icof := getOwenT0Icof(tjd)
	tau[0] = 0
	tau[1] = (tjd - t0) / 36525.0 / 40.0
	for i := 2; i <= 9; i++ {
		tau[i] = tau[1] * tau[i-1]
	}
	k[0] = 1
	k[1] = tau[1]
	k[2] = 2*tau[2] - 1
	k[3] = 4*tau[3] - 3*tau[1]
	k[4] = 8*tau[4] - 8*tau[2] + 1
	k[5] = 16*tau[5] - 20*tau[3] + 5*tau[1]
	k[6] = 32*tau[6] - 48*tau[4] + 18*tau[2] - 1
	k[7] = 64*tau[7] - 112*tau[5] + 56*tau[3] - 7*tau[1]
	k[8] = 128*tau[8] - 256*tau[6] + 160*tau[4] - 32*tau[2] + 1
	k[9] = 256*tau[9] - 576*tau[7] + 432*tau[5] - 120*tau[3] + 9*tau[1]
	// Calculate precession angles
	for i := 0; i < 10; i++ {
		psia += k[i] * owenPsiaCoef[icof][i]
		oma += k[i] * owenOmaCoef[icof][i]
		chia += k[i] * owenChiaCoef[icof][i]
	}
	// Port: removed correction for JPL horizons
	eps0 = 84381.448 / 3600.0
	eps0 *= DEGTORAD
	psia *= DEGTORAD
	chia *= DEGTORAD
	oma *= DEGTORAD
	coseps0 = math.Cos(eps0)
	sineps0 = math.Sin(eps0)
	coschia = math.Cos(chia)
	sinchia = math.Sin(chia)
	cospsia = math.Cos(psia)
	sinpsia = math.Sin(psia)
	cosoma = math.Cos(oma)
	sinoma = math.Sin(oma)
	rp[0] = coschia*cospsia + sinchia*cosoma*sinpsia
	rp[1] = (-coschia*sinpsia+sinchia*cosoma*cospsia)*coseps0 + sinchia*sinoma*sineps0
	rp[2] = (-coschia*sinpsia+sinchia*cosoma*cospsia)*sineps0 - sinchia*sinoma*coseps0
	rp[3] = -sinchia*cospsia + coschia*cosoma*sinpsia
	rp[4] = (sinchia*sinpsia+coschia*cosoma*cospsia)*coseps0 + coschia*sinoma*sineps0
	rp[5] = (sinchia*sinpsia+coschia*cosoma*cospsia)*sineps0 - coschia*sinoma*coseps0
	rp[6] = sinoma * sinpsia
	rp[7] = sinoma*cospsia*coseps0 - cosoma*sineps0
	rp[8] = sinoma*cospsia*sineps0 + cosoma*coseps0
}

// ===== 0887 ===== swi_epsiln swephlib.c-0887 =======================================================================

func swiEpsiln(J float64, iflag int32) float64 {
	var eps float64
	var T float64

	precModel := swed.AstroModels[SE_MODEL_PREC_LONGTERM]
	precModelShort := swed.AstroModels[SE_MODEL_PREC_SHORTTERM]

	// Port: removed references to JPL
	if precModel == 0 {
		precModel = SEMOD_PREC_DEFAULT
	}
	if precModelShort == 0 {
		precModelShort = SEMOD_PREC_DEFAULT_SHORT
	}
	T = (J - 2451545.0) / 36525.0

	if precModelShort == SEMOD_PREC_IAU_1976 && math.Abs(T) <= PREC_IAU_1976_CTIES {
		eps = (((1.813e-3*T-5.9e-4)*T-46.8150)*T + 84381.448) * DEGTORAD / 3600
	} else if precModel == SEMOD_PREC_IAU_1976 {
		eps = (((1.813e-3*T-5.9e-4)*T-46.8150)*T + 84381.448) * DEGTORAD / 3600
	} else if precModelShort == SEMOD_PREC_IAU_2000 && math.Abs(T) <= PREC_IAU_2000_CTIES {
		eps = (((1.813e-3*T-5.9e-4)*T-46.84024)*T + 84381.406) * DEGTORAD / 3600
	} else if precModel == SEMOD_PREC_IAU_2000 {
		eps = (((1.813e-3*T-5.9e-4)*T-46.84024)*T + 84381.406) * DEGTORAD / 3600
	} else if precModelShort == SEMOD_PREC_IAU_2006 && math.Abs(T) <= PREC_IAU_2006_CTIES {
		eps = (((((-4.34e-8*T-5.76e-7)*T+2.0034e-3)*T-1.831e-4)*T-46.836769)*T + 84381.406) * DEGTORAD / 3600.0
	} else if precModel == SEMOD_PREC_NEWCOMB {
		Tn := (J - 2396758.0) / 36525.0
		eps = (0.0017*Tn*Tn*Tn - 0.0085*Tn*Tn - 46.837*Tn + 84451.68) * DEGTORAD / 3600.0
	} else if precModel == SEMOD_PREC_IAU_2006 {
		eps = (((((-4.34e-8*T-5.76e-7)*T+2.0034e-3)*T-1.831e-4)*T-46.836769)*T + 84381.406) * DEGTORAD / 3600.0
	} else if precModel == SEMOD_PREC_BRETAGNON_2003 {
		eps = ((((((-3e-11*T-2.48e-8)*T-5.23e-7)*T+1.99911e-3)*T-1.667e-4)*T-46.836051)*T + 84381.40880) * DEGTORAD / 3600.0
	} else if precModel == SEMOD_PREC_SIMON_1994 {
		eps = (((((2.5e-8*T-5.1e-7)*T+1.9989e-3)*T-1.52e-4)*T-46.80927)*T + 84381.412) * DEGTORAD / 3600.0
	} else if precModel == SEMOD_PREC_WILLIAMS_1994 {
		eps = ((((-1.0e-6*T+2.0e-3)*T-1.74e-4)*T-46.833960)*T + 84381.409) * DEGTORAD / 3600.0
	} else if precModel == SEMOD_PREC_LASKAR_1986 || precModel == SEMOD_PREC_WILL_EPS_LASK {
		T /= 10.0
		eps = ((((((((2.45e-10*T+5.79e-9)*T+2.787e-7)*T+
			7.12e-7)*T-3.905e-5)*T-2.4967e-3)*T-
			5.138e-3)*T+1.99925)*T-0.0155)*T - 468.093*T +
			84381.448
		eps *= DEGTORAD / 3600.0
	} else if precModel == SEMOD_PREC_OWEN_1990 {
		eps = epsilnOwen1986(J)
		eps *= DEGTORAD
	} else { // SEMOD_PREC_VONDRAK_2011
		_, eps = swiLdpPeps(J)
	}
	return eps
}

// ===== 1169 ===== data for precess2 swephlib.c-1170 ================================================================

// In WILLIAMS and SIMON, Laskar's terms of order higher than t^4 have been retained, because Simon et al mention that
// the solution is the same except for the lower order terms.

// SEMOD_PREC_WILLIAMS_1994
var pAcofWilliams = []float64{
	-8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3,
	-0.235316, 0.076, 110.5407, 50287.70000,
}
var nodecofWilliams = []float64{
	6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 1.9e-10,
	-3.54e-9, -1.8103e-7, 1.26e-7, 7.436169e-5,
	-0.04207794833, 3.052115282424,
}
var inclcofWilliams = []float64{
	1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11,
	-5.4000441e-11, 1.32115526e-9, -6.012e-7, -1.62442e-5,
	0.00227850649, 0.0,
}

// SEMOD_PREC_SIMON_1994
var pAcofSimon = []float64{
	-8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3,
	-0.235316, 0.07732, 111.2022, 50288.200,
}
var nodecofSimon = []float64{
	6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 1.9e-10,
	-3.54e-9, -1.8103e-7, 2.579e-8, 7.4379679e-5,
	-0.0420782900, 3.0521126906,
}
var inclcofSimon = []float64{
	1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11,
	-5.4000441e-11, 1.32115526e-9, -5.99908e-7, -1.624383e-5,
	0.002278492868, 0.0,
}

// SEMOD_PREC_LASKAR_1986
var pAcofLaskar = []float64{
	-8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3,
	-0.235316, 0.07732, 111.1971, 50290.966,
}

// Node and inclination of the earth's orbit computed from Laskar's data as done in Bretagnon and Francou's paper.
// Units are radians.
var nodecofLaskar = []float64{
	6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 6.3190131e-10,
	-3.48388152e-9, -1.813065896e-7, 2.75036225e-8, 7.4394531426e-5,
	-0.042078604317, 3.052112654975,
}
var inclcofLaskar = []float64{
	1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11,
	-5.4000441e-11, 1.32115526e-9, -5.998737027e-7, -1.6242797091e-5,
	0.002278495537, 0.0,
}

// ===== 1219 ===== precess_2 swephlib.c-1219 ========================================================================

func precess2(R []float64, J float64, iflag int32, direction int, precMethod int) int {
	const J2000 = 2451545.0 // You'll need to define this constant

	if J == J2000 {
		return 0
	}

	var pAcof, nodecof, inclcof []float64

	switch precMethod {
	case SEMOD_PREC_LASKAR_1986:
		pAcof = pAcofLaskar
		nodecof = nodecofLaskar
		inclcof = inclcofLaskar
	case SEMOD_PREC_SIMON_1994:
		pAcof = pAcofSimon
		nodecof = nodecofSimon
		inclcof = inclcofSimon
	case SEMOD_PREC_WILLIAMS_1994:
		pAcof = pAcofWilliams
		nodecof = nodecofWilliams
		inclcof = inclcofWilliams
	default:
		pAcof = pAcofLaskar
		nodecof = nodecofLaskar
		inclcof = inclcofLaskar
	}

	T := (J - J2000) / 36525.0
	// Implementation by elementary rotations using Laskar's expansions. First rotate about the x axis from the initial
	// equator to the ecliptic. (The input is equatorial.)
	var eps float64
	if direction == 1 {
		eps = swiEpsiln(J, iflag) // To J2000
	} else {
		eps = swiEpsiln(J2000, iflag) // From J2000
	}
	sineps := math.Sin(eps)
	coseps := math.Cos(eps)
	x := make([]float64, 3)
	x[0] = R[0]
	z := coseps*R[1] + sineps*R[2]
	x[2] = -sineps*R[1] + coseps*R[2]
	x[1] = z
	// Precession in longitude
	T /= 10.0 // thousands of years
	pA := pAcof[0]
	for i := 0; i < 9; i++ {
		pA = pA*T + pAcof[i+1]
	}
	pA *= DEGTORAD / 3600 * T
	// Node of the moving ecliptic on the J2000 ecliptic
	W := nodecof[0]
	for i := 0; i < 10; i++ {
		W = W*T + nodecof[i+1]
	}
	// Rotate about z axis to the node
	var z1 float64
	if direction == 1 {
		z1 = W + pA
	} else {
		z1 = W
	}
	B := math.Cos(z1)
	A := math.Sin(z1)
	z = B*x[0] + A*x[1]
	x[1] = -A*x[0] + B*x[1]
	x[0] = z
	// Rotate about new x axis by the inclination of the moving ecliptic on the J2000 ecliptic.
	z = inclcof[0]
	for i := 0; i < 10; i++ {
		z = z*T + inclcof[i+1]
	}
	if direction == 1 {
		z = -z
	}
	B = math.Cos(z)
	A = math.Sin(z)
	z = B*x[1] + A*x[2]
	x[2] = -A*x[1] + B*x[2]
	x[1] = z
	// Rotate about new z axis back from the node
	if direction == 1 {
		z = -W
	} else {
		z = -W - pA
	}
	B = math.Cos(z)
	A = math.Sin(z)
	z = B*x[0] + A*x[1]
	x[1] = -A*x[0] + B*x[1]
	x[0] = z
	// Rotate about x axis to final equator
	if direction == 1 {
		eps = swiEpsiln(J2000, iflag)
	} else {
		eps = swiEpsiln(J, iflag)
	}
	sineps = math.Sin(eps)
	coseps = math.Cos(eps)
	z = coseps*x[1] - sineps*x[2]
	x[2] = sineps*x[1] + coseps*x[2]
	x[1] = z
	for i := 0; i < 3; i++ {
		R[i] = x[i]
	}
	return 0
}

// ===== 1487 ===== data for calcNutationIau1980 swephlib.c-1487 =====================================================

// NutationTerms represents the IAU 1980 nutation series
var nt = []int16{
	// LS and OC are units of 0.0001"
	// LS2 and OC2 are units of 0.00001"
	// MM,MS,FF,DD,OM, LS, LS2,OC, OC2
	0, 0, 0, 0, 2, 2062, 2, -895, 5,
	-2, 0, 2, 0, 1, 46, 0, -24, 0,
	2, 0, -2, 0, 0, 11, 0, 0, 0,
	-2, 0, 2, 0, 2, -3, 0, 1, 0,
	1, -1, 0, -1, 0, -3, 0, 0, 0,
	0, -2, 2, -2, 1, -2, 0, 1, 0,
	2, 0, -2, 0, 1, 1, 0, 0, 0,
	0, 0, 2, -2, 2, -13187, -16, 5736, -31,
	0, 1, 0, 0, 0, 1426, -34, 54, -1,
	0, 1, 2, -2, 2, -517, 12, 224, -6,
	0, -1, 2, -2, 2, 217, -5, -95, 3,
	0, 0, 2, -2, 1, 129, 1, -70, 0,
	2, 0, 0, -2, 0, 48, 0, 1, 0,
	0, 0, 2, -2, 0, -22, 0, 0, 0,
	0, 2, 0, 0, 0, 17, -1, 0, 0,
	0, 1, 0, 0, 1, -15, 0, 9, 0,
	0, 2, 2, -2, 2, -16, 1, 7, 0,
	0, -1, 0, 0, 1, -12, 0, 6, 0,
	-2, 0, 0, 2, 1, -6, 0, 3, 0,
	0, -1, 2, -2, 1, -5, 0, 3, 0,
	2, 0, 0, -2, 1, 4, 0, -2, 0,
	0, 1, 2, -2, 1, 4, 0, -2, 0,
	1, 0, 0, -1, 0, -4, 0, 0, 0,
	2, 1, 0, -2, 0, 1, 0, 0, 0,
	0, 0, -2, 2, 1, 1, 0, 0, 0,
	0, 1, -2, 2, 0, -1, 0, 0, 0,
	0, 1, 0, 0, 2, 1, 0, 0, 0,
	-1, 0, 0, 1, 1, 1, 0, 0, 0,
	0, 1, 2, -2, 0, -1, 0, 0, 0,
	0, 0, 2, 0, 2, -2274, -2, 977, -5,
	1, 0, 0, 0, 0, 712, 1, -7, 0,
	0, 0, 2, 0, 1, -386, -4, 200, 0,
	1, 0, 2, 0, 2, -301, 0, 129, -1,
	1, 0, 0, -2, 0, -158, 0, -1, 0,
	-1, 0, 2, 0, 2, 123, 0, -53, 0,
	0, 0, 0, 2, 0, 63, 0, -2, 0,
	1, 0, 0, 0, 1, 63, 1, -33, 0,
	-1, 0, 0, 0, 1, -58, -1, 32, 0,
	-1, 0, 2, 2, 2, -59, 0, 26, 0,
	1, 0, 2, 0, 1, -51, 0, 27, 0,
	0, 0, 2, 2, 2, -38, 0, 16, 0,
	2, 0, 0, 0, 0, 29, 0, -1, 0,
	1, 0, 2, -2, 2, 29, 0, -12, 0,
	2, 0, 2, 0, 2, -31, 0, 13, 0,
	0, 0, 2, 0, 0, 26, 0, -1, 0,
	-1, 0, 2, 0, 1, 21, 0, -10, 0,
	-1, 0, 0, 2, 1, 16, 0, -8, 0,
	1, 0, 0, -2, 1, -13, 0, 7, 0,
	-1, 0, 2, 2, 1, -10, 0, 5, 0,
	1, 1, 0, -2, 0, -7, 0, 0, 0,
	0, 1, 2, 0, 2, 7, 0, -3, 0,
	0, -1, 2, 0, 2, -7, 0, 3, 0,
	1, 0, 2, 2, 2, -8, 0, 3, 0,
	1, 0, 0, 2, 0, 6, 0, 0, 0,
	2, 0, 2, -2, 2, 6, 0, -3, 0,
	0, 0, 0, 2, 1, -6, 0, 3, 0,
	0, 0, 2, 2, 1, -7, 0, 3, 0,
	1, 0, 2, -2, 1, 6, 0, -3, 0,
	0, 0, 0, -2, 1, -5, 0, 3, 0,
	1, -1, 0, 0, 0, 5, 0, 0, 0,
	2, 0, 2, 0, 1, -5, 0, 3, 0,
	0, 1, 0, -2, 0, -4, 0, 0, 0,
	1, 0, -2, 0, 0, 4, 0, 0, 0,
	0, 0, 0, 1, 0, -4, 0, 0, 0,
	1, 1, 0, 0, 0, -3, 0, 0, 0,
	1, 0, 2, 0, 0, 3, 0, 0, 0,
	1, -1, 2, 0, 2, -3, 0, 1, 0,
	-1, -1, 2, 2, 2, -3, 0, 1, 0,
	-2, 0, 0, 0, 1, -2, 0, 1, 0,
	3, 0, 2, 0, 2, -3, 0, 1, 0,
	0, -1, 2, 2, 2, -3, 0, 1, 0,
	1, 1, 2, 0, 2, 2, 0, -1, 0,
	-1, 0, 2, -2, 1, -2, 0, 1, 0,
	2, 0, 0, 0, 1, 2, 0, -1, 0,
	1, 0, 0, 0, 2, -2, 0, 1, 0,
	3, 0, 0, 0, 0, 2, 0, 0, 0,
	0, 0, 2, 1, 2, 2, 0, -1, 0,
	-1, 0, 0, 0, 2, 1, 0, -1, 0,

	1, 0, 0, -4, 0, -1, 0, 0, 0,
	-2, 0, 2, 2, 2, 1, 0, -1, 0,
	-1, 0, 2, 4, 2, -2, 0, 1, 0,
	2, 0, 0, -4, 0, -1, 0, 0, 0,
	1, 1, 2, -2, 2, 1, 0, -1, 0,
	1, 0, 2, 2, 1, -1, 0, 1, 0,
	-2, 0, 2, 4, 2, -1, 0, 1, 0,
	-1, 0, 4, 0, 2, 1, 0, 0, 0,
	1, -1, 0, -2, 0, 1, 0, 0, 0,
	2, 0, 2, -2, 1, 1, 0, -1, 0,
	2, 0, 2, 2, 2, -1, 0, 0, 0,
	1, 0, 0, 2, 1, -1, 0, 0, 0,
	0, 0, 4, -2, 2, 1, 0, 0, 0,
	3, 0, 2, -2, 2, 1, 0, 0, 0,
	1, 0, 2, -2, 0, -1, 0, 0, 0,
	0, 1, 2, 0, 1, 1, 0, 0, 0,
	-1, -1, 0, 2, 1, 1, 0, 0, 0,
	0, 0, -2, 0, 1, -1, 0, 0, 0,
	0, 0, 2, -1, 2, -1, 0, 0, 0,
	0, 1, 0, 2, 0, -1, 0, 0, 0,
	1, 0, -2, -2, 0, -1, 0, 0, 0,
	0, -1, 2, 0, 1, -1, 0, 0, 0,
	1, 1, 0, -2, 1, -1, 0, 0, 0,
	1, 0, -2, 2, 0, -1, 0, 0, 0,
	2, 0, 0, 2, 0, 1, 0, 0, 0,
	0, 0, 2, 4, 2, -1, 0, 0, 0,
	0, 1, 0, 1, 0, 1, 0, 0, 0,
	// if NUT_CORR_1987  switch is handled in function calc_nutation_iau1980()
	// corrections to IAU 1980 nutation series by Herring 1987
	//             in 0.00001" !!!
	//              LS      OC
	101, 0, 0, 0, 1, -725, 0, 213, 0,
	101, 1, 0, 0, 0, 523, 0, 208, 0,
	101, 0, 2, -2, 2, 102, 0, -41, 0,
	101, 0, 2, 0, 2, -81, 0, 32, 0,
	// LC and OS
	102, 0, 0, 0, 1, 417, 0, 224, 0,
	102, 1, 0, 0, 0, 61, 0, -24, 0,
	102, 0, 2, -2, 2, -118, 0, -47, 0,

	ENDMARK,
}

// ===== 1615 ===== calc_nutation_iau1980 swephlib.c-1615. see also constants as defined before (slice nt[]) =========

func calcNutationIau1980(J float64, nutlo []float64) int {
	// Arrays to hold sines and cosines of multiple angles
	ss := [5][8]float64{}
	cc := [5][8]float64{}
	var args [5]float64
	var ns [5]int

	nutModel := swed.AstroModels[SE_MODEL_NUT]
	if nutModel == 0 {
		nutModel = SEMOD_NUT_DEFAULT
	}

	// Julian centuries from 2000 January 1.5, barycentric dynamical time
	T := (J - 2451545.0) / 36525.0
	T2 := T * T

	// Fundamental arguments in the FK5 reference system.
	// The coefficients, originally given to 0.001", are converted here to degrees.

	// longitude of the mean ascending node of the lunar orbit on the ecliptic, measured from the mean equinox of date
	OM := -6962890.539*T + 450160.280 + (0.008*T+7.455)*T2
	OM = SweDegnorm(OM/3600) * DEGTORAD
	// mean longitude of the Sun minus the mean longitude of the Sun's perigee
	MS := 129596581.224*T + 1287099.804 - (0.012*T+0.577)*T2
	MS = SweDegnorm(MS/3600) * DEGTORAD
	// mean longitude of the Moon minus the mean longitude of the Moon's perigee
	MM := 1717915922.633*T + 485866.733 + (0.064*T+31.310)*T2
	MM = SweDegnorm(MM/3600) * DEGTORAD
	// mean longitude of the Moon minus the mean longitude of the Moon's node
	FF := 1739527263.137*T + 335778.877 + (0.011*T-13.257)*T2
	FF = SweDegnorm(FF/3600) * DEGTORAD
	// mean elongation of the Moon from the Sun
	DD := 1602961601.328*T + 1072261.307 + (0.019*T-6.891)*T2
	DD = SweDegnorm(DD/3600) * DEGTORAD
	args[0] = MM
	ns[0] = 3
	args[1] = MS
	ns[1] = 2
	args[2] = FF
	ns[2] = 4
	args[3] = DD
	ns[3] = 4
	args[4] = OM
	ns[4] = 2

	// Calculate sin(i*MM), etc. for needed multiple angles
	for k := 0; k <= 4; k++ {
		arg := args[k]
		n := ns[k]
		su := math.Sin(arg)
		cu := math.Cos(arg)
		ss[k][0] = su // sin(L)
		cc[k][0] = cu // cos(L)
		sv := 2.0 * su * cu
		cv := cu*cu - su*su
		ss[k][1] = sv // sin(2L)
		cc[k][1] = cv

		for i := 2; i < n; i++ {
			s := su*cv + cu*sv
			cv = cu*cv - su*sv
			sv = s
			ss[k][i] = sv // sin(i+1 L)
			cc[k][i] = cv
		}
	}

	// First terms, not in table
	C := (-0.01742*T - 17.1996) * ss[4][0] // sin(OM)
	D := (0.00089*T + 9.2025) * cc[4][0]   // cos(OM)

	// Process nutation terms
	for i := 0; nt[i] != ENDMARK; i += 9 {
		if nutModel != SEMOD_NUT_IAU_CORR_1987 && (nt[i] == 101 || nt[i] == 102) {
			continue
		}
		// Argument of sine and cosine
		k1 := 0
		cv := 0.0
		sv := 0.0
		for m := 0; m < 5; m++ {
			j := nt[i+m]
			if j > 100 {
				j = 0 // nt[0] is a flag
			}
			if j != 0 {
				k := j
				if j < 0 {
					k = -k
				}
				su := ss[m][k-1] // sin(k*angle)
				if j < 0 {
					su = -su
				}
				cu := cc[m][k-1]
				if k1 == 0 { // set first angle
					sv = su
					cv = cu
					k1 = 1
				} else { // combine angles
					sw := su*cv + cu*sv
					cv = cu*cv - su*sv
					sv = sw
				}
			}
		}
		// Longitude coefficient, in 0.0001"
		f := float64(nt[i+5]) * 0.0001
		if nt[i+6] != 0 {
			f += 0.00001 * T * float64(nt[i+6])
		}
		// Obliquity coefficient, in 0.0001"
		g := float64(nt[i+7]) * 0.0001
		if nt[i+8] != 0 {
			g += 0.00001 * T * float64(nt[i+8])
		}
		if nt[i] >= 100 { // coefficients in 0.00001"
			f *= 0.1
			g *= 0.1
		}
		// Accumulate the terms
		if nt[i] != 102 {
			C += f * sv
			D += g * cv
		} else { // cos for nutl and sin for nuto
			C += f * cv
			D += g * sv
		}
	}
	// Save answers, expressed in radians
	nutlo[0] = C
	nutlo[1] = D
	return 0
}

// ===== 1766 ===== calc_nutation_iau2000ab swephlib.c-1766 ===========================================================
//
// Nutation IAU 2000A model (MHB2000 luni-solar and planetary nutation, without free core nutation)
// Function returns nutation in longitude and obliquity in radians with respect to the equinox of date.
// For the obliquity of the ecliptic the calculation of Lieske & al. (1977) must be used.
// The precision in recent years is about 0.001 arc seconds.
// The calculation includes luni-solar and planetary nutation. Free core nutation, which cannot be predicted, is
// omitted, the error being of the order of a few 0.0001 arc seconds.
// References:
// Capitaine, N., Wallace, P.T., Chapront, J., A & A 432, 366 (2005).
// Chapront, J., Chapront-Touze, M. & Francou, G., A & A 387, 700 (2002).
// Lieske, J.H., Lederle, T., Fricke, W. & Morando, B., "Expressions for the precession quantities based upon the IAU
// (1976) System of Astronomical Constants", A & A 58, 1-16 (1977).
// Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation and precession   New nutation series for nonrigid
// Earth and insights into the Earth's interior", J.Geophys.Res., 107, B4, 2002.
// Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M., Francou, G., Laskar, J., A & A 282, 663-683 (1994).
// Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M., A & A Supp. Ser. 135, 111 (1999).
// Wallace, P.T., "Software for Implementing the IAU 2000 Resolutions", in IERS Workshop 5.1 (2002).
// Nutation IAU 2000A series in: Kaplan, G.H., United States Naval Observatory Circular No. 179 (Oct. 2005)
// aa.usno.navy.mil/publications/docs/Circular_179.html
// MHB2000 code at
// - ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
// - http://www.iau-sofa.rl.ac.uk/2005_0901/Downloads.html

func calcNutationIau2000ab(J float64, nutlo []float64) int {
	var i, j, k, inls int
	var M, SM, F, D, OM float64
	var AL, ALSU, AF, AD, AOM, APA float64
	var ALME, ALVE, ALEA, ALMA, ALJU, ALSA, ALUR, ALNE float64
	var darg, sinarg, cosarg float64
	var dpsi, deps float64
	T := (J - J2000) / 36525.0
	nutModel := swed.AstroModels[SE_MODEL_NUT]
	if nutModel == 0 {
		nutModel = SEMOD_NUT_DEFAULT
	}
	// luni-solar nutation
	// Fundamental arguments, Simon & al. (1994)
	// Mean anomaly of the Moon
	M = SweDegnorm((485868.249036+
		T*(1717915923.2178+
			T*(31.8792+
				T*(0.051635+
					T*(-0.00024470)))))/3600.0) * DEGTORAD

	// Mean anomaly of the Sun
	SM = SweDegnorm((1287104.79305+
		T*(129596581.0481+
			T*(-0.5532+
				T*(0.000136+
					T*(-0.00001149)))))/3600.0) * DEGTORAD

	// Mean argument of the latitude of the Moon
	F = SweDegnorm((335779.526232+
		T*(1739527262.8478+
			T*(-12.7512+
				T*(-0.001037+
					T*(0.00000417)))))/3600.0) * DEGTORAD

	// Mean elongation of the Moon from the Sun
	D = SweDegnorm((1072260.70369+
		T*(1602961601.2090+
			T*(-6.3706+
				T*(0.006593+
					T*(-0.00003169)))))/3600.0) * DEGTORAD

	// Mean longitude of the ascending node of the Moon
	OM = SweDegnorm((450160.398036+
		T*(-6962890.5431+
			T*(7.4722+
				T*(0.007702+
					T*(-0.00005939)))))/3600.0) * DEGTORAD

	// luni-solar nutation series, in reverse order, starting with small terms
	if nutModel == SEMOD_NUT_IAU_2000B {
		inls = NLS_2000B
	} else {
		inls = NLS
	}

	// Process luni-solar nutation terms
	dpsi = 0
	deps = 0
	for i = inls - 1; i >= 0; i-- {
		j = i * 5
		darg = SweRadnorm(float64(nls[j+0])*M +
			float64(nls[j+1])*SM +
			float64(nls[j+2])*F +
			float64(nls[j+3])*D +
			float64(nls[j+4])*OM)
		sinarg = math.Sin(darg)
		cosarg = math.Cos(darg)
		k = i * 6
		dpsi += (float64(cls[k+0])+float64(cls[k+1])*T)*sinarg + float64(cls[k+2])*cosarg
		deps += (float64(cls[k+3])+float64(cls[k+4])*T)*cosarg + float64(cls[k+5])*sinarg
	}

	nutlo[0] = dpsi * O1MAS2DEG
	nutlo[1] = deps * O1MAS2DEG

	if nutModel == SEMOD_NUT_IAU_2000A {
		// Planetary nutation
		AL = SweRadnorm(2.35555598 + 8328.6914269554*T)
		ALSU = SweRadnorm(6.24006013 + 628.301955*T)
		AF = SweRadnorm(1.627905234 + 8433.466158131*T)
		AD = SweRadnorm(5.198466741 + 7771.3771468121*T)
		AOM = SweRadnorm(2.18243920 - 33.757045*T)

		// Planetary longitudes
		ALME = SweRadnorm(4.402608842 + 2608.7903141574*T)
		ALVE = SweRadnorm(3.176146697 + 1021.3285546211*T)
		ALEA = SweRadnorm(1.753470314 + 628.3075849991*T)
		ALMA = SweRadnorm(6.203480913 + 334.0612426700*T)
		ALJU = SweRadnorm(0.599546497 + 52.9690962641*T)
		ALSA = SweRadnorm(0.874016757 + 21.3299104960*T)
		ALUR = SweRadnorm(5.481293871 + 7.4781598567*T)
		ALNE = SweRadnorm(5.321159000 + 3.8127774000*T)

		// General accumulated precession in longitude
		APA = (0.02438175 + 0.00000538691*T) * T

		// Process planetary nutation terms
		dpsi = 0
		deps = 0
		for i = NPL - 1; i >= 0; i-- {
			j = i * 14
			darg = SweRadnorm(float64(npl[j+0])*AL +
				float64(npl[j+1])*ALSU +
				float64(npl[j+2])*AF +
				float64(npl[j+3])*AD +
				float64(npl[j+4])*AOM +
				float64(npl[j+5])*ALME +
				float64(npl[j+6])*ALVE +
				float64(npl[j+7])*ALEA +
				float64(npl[j+8])*ALMA +
				float64(npl[j+9])*ALJU +
				float64(npl[j+10])*ALSA +
				float64(npl[j+11])*ALUR +
				float64(npl[j+12])*ALNE +
				float64(npl[j+13])*APA)

			k = i * 4
			sinarg = math.Sin(darg)
			cosarg = math.Cos(darg)
			dpsi += float64(icpl[k+0])*sinarg + float64(icpl[k+1])*cosarg
			deps += float64(icpl[k+2])*sinarg + float64(icpl[k+3])*cosarg
		}

		nutlo[0] += dpsi * O1MAS2DEG
		nutlo[1] += deps * O1MAS2DEG

		// Changes required by adoption of P03 precession
		dpsi = -8.1*math.Sin(OM) - 0.6*math.Sin(2*F-2*D+2*OM)
		dpsi += T * (47.8*math.Sin(OM) + 3.7*math.Sin(2*F-2*D+2*OM) +
			0.6*math.Sin(2*F+2*OM) - 0.6*math.Sin(2*OM))
		deps = T * (-25.6*math.Cos(OM) - 1.6*math.Cos(2*F-2*D+2*OM))

		nutlo[0] += dpsi / (3600.0 * 1000000.0)
		nutlo[1] += deps / (3600.0 * 1000000.0)
	}

	nutlo[0] *= DEGTORAD
	nutlo[1] *= DEGTORAD

	return 0
}

// ===== 1947 ===== calc_nutation_woolard swephlib.c-1947 ============================================================

// calcNutationWoolard an incomplete implementation of nutation Woolard 1953
func calcNutationWoolard(J float64, nutlo []float64) int {
	mjd := J - J1900
	t := mjd / 36525.0
	t2 := t * t
	a := 100.0021358 * t
	b := 360.0 * (a - float64(int64(a)))
	ls := 279.697 + 0.000303*t2 + b // Sun's mean longitude
	a = 1336.855231 * t
	b = 360.0 * (a - float64(int64(a)))
	ld := 270.434 - 0.001133*t2 + b // Moon's mean longitude
	a = 99.99736056000026 * t
	b = 360.0 * (a - float64(int64(a)))
	ms := 358.476 - 0.00015*t2 + b // Sun's mean anomaly
	a = 13255523.59 * t
	b = 360.0 * (a - float64(int64(a)))
	md := 296.105 + 0.009192*t2 + b // Moon's mean anomaly
	a = 5.372616667 * t
	b = 360.0 * (a - float64(int64(a)))
	nm := 259.183 + 0.002078*t2 - b // Longitude of moon's ascending node
	// Convert to radians for trigonometric functions
	tls := 2 * ls * DEGTORAD
	nm = nm * DEGTORAD
	tnm := 2 * nm
	ms = ms * DEGTORAD
	tld := 2 * ld * DEGTORAD
	md = md * DEGTORAD
	// Calculate delta psi and eps in arcseconds
	dpsi := (-17.2327-0.01737*t)*math.Sin(nm) +
		(-1.2729-0.00013*t)*math.Sin(tls) +
		0.2088*math.Sin(tnm) -
		0.2037*math.Sin(tld) +
		(0.1261-0.00031*t)*math.Sin(ms) +
		0.0675*math.Sin(md) -
		(0.0497-0.00012*t)*math.Sin(tls+ms) -
		0.0342*math.Sin(tld-nm) -
		0.0261*math.Sin(tld+md) +
		0.0214*math.Sin(tls-ms) -
		0.0149*math.Sin(tls-tld+md) +
		0.0124*math.Sin(tls-nm) +
		0.0114*math.Sin(tld-md)
	deps := (9.21+0.00091*t)*math.Cos(nm) +
		(0.5522-0.00029*t)*math.Cos(tls) -
		0.0904*math.Cos(tnm) +
		0.0884*math.Cos(tld) +
		0.0216*math.Cos(tls+ms) +
		0.0183*math.Cos(tld-nm) +
		0.0113*math.Cos(tld+md) -
		0.0093*math.Cos(tls-ms) -
		0.0066*math.Cos(tls-nm)
	// Convert to radians
	dpsi = dpsi / 3600.0 * DEGTORAD
	deps = deps / 3600.0 * DEGTORAD
	nutlo[1] = deps
	nutlo[0] = dpsi
	return OK
}

// ===== 2069 ===== calc_nutation swephlib.c-2069 ====================================================================

func calcNutation(J float64, iflag int32, nutlo []float64) int {
	// Port: removed logic for JPL
	nutModel := swed.AstroModels[SE_MODEL_NUT]
	jplhoraModel := swed.AstroModels[SE_MODEL_JPLHORA_MODE]

	if nutModel == 0 {
		nutModel = SEMOD_NUT_DEFAULT
	}
	if nutModel == SEMOD_NUT_IAU_1980 || nutModel == SEMOD_NUT_IAU_CORR_1987 {
		calcNutationIau1980(J, nutlo)
	} else if nutModel == SEMOD_NUT_IAU_2000A || nutModel == SEMOD_NUT_IAU_2000B {
		calcNutationIau2000ab(J, nutlo)
		if (iflag&SEFLG_JPLHOR_APPROX) != 0 && jplhoraModel == SEMOD_JPLHORA_2 {
			nutlo[0] += -41.7750 / 3600.0 / 1000.0 * DEGTORAD
			nutlo[1] += -6.8192 / 3600.0 / 1000.0 * DEGTORAD
		}
	} else if nutModel == SEMOD_NUT_WOOLARD {
		calcNutationWoolard(J, nutlo)
	}
	return OK
}

// ===== 2116 ===== quadratic_intp swephlib.c-2116 ===================================================================

func quadraticIntp(ym, y0, yp, x float64) float64 {
	c := y0
	b := (yp - ym) / 2.0
	a := (yp+ym)/2.0 - c
	y := a*x*x + b*x + c
	return y
}

// ===== 2126 ===== swi_nutation swephlib.c-2126 =====================================================================

func swiNutation(tjd float64, iflag int32, nutlo []float64) int {
	var retc = OK
	dnut := make([]float64, 2)

	if !swed.DoInterpolateNut {
		retc = calcNutation(tjd, iflag, nutlo)
		// from interpolation, with three data points in 1-day steps;
		// maximum error is about 3 mas
	} else {
		// Check if precalculated data points are available
		if tjd < swed.Interpol.TjdNut2 && tjd > swed.Interpol.TjdNut0 {
			// Interpolate between existing points
			dx := (tjd - swed.Interpol.TjdNut0) - 1.0
			nutlo[0] = quadraticIntp(
				swed.Interpol.NutDpsi0,
				swed.Interpol.NutDpsi1,
				swed.Interpol.NutDpsi2,
				dx,
			)
			nutlo[1] = quadraticIntp(
				swed.Interpol.NutDeps0,
				swed.Interpol.NutDeps1,
				swed.Interpol.NutDeps2,
				dx,
			)
		} else {
			swed.Interpol.TjdNut0 = tjd - 1.0 // one day earlier
			swed.Interpol.TjdNut2 = tjd + 1.0 // one day later
			retc = calcNutation(swed.Interpol.TjdNut0, iflag, dnut)
			if retc == ERR {
				return ERR
			}
			swed.Interpol.NutDpsi0 = dnut[0]
			swed.Interpol.NutDeps0 = dnut[1]
			retc = calcNutation(swed.Interpol.TjdNut2, iflag, dnut)
			if retc == ERR {
				return ERR
			}
			swed.Interpol.NutDpsi2 = dnut[0]
			swed.Interpol.NutDeps2 = dnut[1]
			retc = calcNutation(tjd, iflag, nutlo)
			if retc == ERR {
				return ERR
			}
			swed.Interpol.NutDpsi1 = nutlo[0]
			swed.Interpol.NutDeps1 = nutlo[1]
		}
	}
	return retc
}

// ===== 2160 ===== constants for swi_approx_jlhor swephlib.c-2160 ===================================================
const (
	OFFSET_JPLHORIZONS = -52.3
	DCOR_RA_JPL_TJD0   = 2437846.5
	NDCOR_RA_JPL       = 51
)

var dcorRaJpl = []float64{
	-51.257, -51.103, -51.065, -51.503, -51.224, -50.796, -51.161, -51.181,
	-50.932, -51.064, -51.182, -51.386, -51.416, -51.428, -51.586, -51.766, -52.038, -52.370,
	-52.553, -52.397, -52.340, -52.676, -52.348, -51.964, -52.444, -52.364, -51.988, -52.212,
	-52.370, -52.523, -52.541, -52.496, -52.590, -52.629, -52.788, -53.014, -53.053, -52.902,
	-52.850, -53.087, -52.635, -52.185, -52.588, -52.292, -51.796, -51.961, -52.055, -52.134,
	-52.165, -52.141, -52.255,
}

// ===== 2172 ===== swi_approx_jplhor swephlib.c-2172 ================================================================

// swiApproxJplhor converts coordinates using JPL Horizons approximation
func swiApproxJplhor(x []float64, tjd float64, iflag int32, backward bool) {
	t := (tjd - DCOR_RA_JPL_TJD0) / 365.25
	dofs := OFFSET_JPLHORIZONS
	jplhoraModel := swed.AstroModels[SE_MODEL_JPLHORA_MODE]
	if jplhoraModel == 0 {
		jplhoraModel = SEMOD_JPLHORA_DEFAULT
	}
	if (iflag & SEFLG_JPLHOR_APPROX) == 0 {
		return
	}
	if jplhoraModel == SEMOD_JPLHORA_2 {
		return
	}
	if t < 0 {
		t = 0
		dofs = dcorRaJpl[0]
	} else if t >= float64(NDCOR_RA_JPL-1) {
		t = float64(NDCOR_RA_JPL)
		dofs = dcorRaJpl[NDCOR_RA_JPL-1]
	} else {
		t0 := int(t)
		t1 := t0 + 1
		dofs = dcorRaJpl[t0]
		dofs = (t-float64(t0))*(dcorRaJpl[t0]-dcorRaJpl[t1]) + dcorRaJpl[t0]
	}
	dofs /= (1000.0 * 3600.0)
	swiCartpol(x)
	if backward {
		x[0] -= dofs * DEGTORAD
	} else {
		x[0] += dofs * DEGTORAD
	}
	swiPolcart(x)
}

// ===== 2204 ===== swi_bias swephlib.c-2204

// SwiBias converts GCRS to J2000
func SwiBias(x []float64, tjd float64, iflag int32, backward bool) {
	xx := make([]float64, 6)
	rb := [3][3]float64{}
	biasModel := swed.AstroModels[SE_MODEL_BIAS]
	jplhoraModel := swed.AstroModels[SE_MODEL_JPLHORA_MODE]
	if biasModel == 0 {
		biasModel = SEMOD_BIAS_DEFAULT
	}
	if jplhoraModel == 0 {
		jplhoraModel = SEMOD_JPLHORA_DEFAULT
	}
	if biasModel == SEMOD_BIAS_NONE {
		return
	}
	if (iflag & SEFLG_JPLHOR_APPROX) != 0 {
		if jplhoraModel == SEMOD_JPLHORA_2 {
			return
		}
		if jplhoraModel == SEMOD_JPLHORA_3 && tjd < DPSI_DEPS_IAU1980_TJD0_HORIZONS {
			return
		}
	}
	if biasModel == SEMOD_BIAS_IAU2006 {
		rb = [3][3]float64{
			{+0.99999999999999412, +0.00000007078368695, -0.00000008056214212},
			{-0.00000007078368961, +0.99999999999999700, -0.00000003306427981},
			{+0.00000008056213978, +0.00000003306428553, +0.99999999999999634},
		}
	} else {
		rb = [3][3]float64{
			{+0.9999999999999942, +0.0000000707827948, -0.0000000805621738},
			{-0.0000000707827974, +0.9999999999999969, -0.0000000330604088},
			{+0.0000000805621715, +0.0000000330604145, +0.9999999999999962},
		}
	}

	if backward {
		swiApproxJplhor(x, tjd, iflag, true)
		for i := 0; i <= 2; i++ {
			xx[i] = x[0]*rb[i][0] +
				x[1]*rb[i][1] +
				x[2]*rb[i][2]
			if (iflag & SEFLG_SPEED) != 0 {
				xx[i+3] = x[3]*rb[i][0] +
					x[4]*rb[i][1] +
					x[5]*rb[i][2]
			}
		}
	} else {
		for i := 0; i <= 2; i++ {
			xx[i] = x[0]*rb[0][i] +
				x[1]*rb[1][i] +
				x[2]*rb[2][i]
			if (iflag & SEFLG_SPEED) != 0 {
				xx[i+3] = x[3]*rb[0][i] +
					x[4]*rb[1][i] +
					x[5]*rb[2][i]
			}
		}
		swiApproxJplhor(xx, tjd, iflag, false)
	}

	for i := 0; i <= 2; i++ {
		x[i] = xx[i]
	}
	if (iflag & SEFLG_SPEED) != 0 {
		for i := 3; i <= 5; i++ {
			x[i] = xx[i]
		}
	}
}

// ===== 2291 ===== swi_icrs2fk5 =====================================================================================

// swiIcrs2fk5 converts coordinates from GCRS (ICRS) to FK5 reference system
func swiIcrs2fk5(x *[6]float64, iflag int32, backward bool) {
	rb := [3][3]float64{
		{+0.9999999999999928, +0.0000001110223287, +0.0000000441180557},
		{-0.0000001110223330, +0.9999999999999891, +0.0000000964779176},
		{-0.0000000441180450, -0.0000000964779225, +0.9999999999999943},
	}
	var xx [6]float64
	if backward {
		// Backward transformation
		for i := 0; i <= 2; i++ {
			xx[i] = x[0]*rb[i][0] +
				x[1]*rb[i][1] +
				x[2]*rb[i][2]

			if (iflag & SEFLG_SPEED) != 0 {
				xx[i+3] = x[3]*rb[i][0] +
					x[4]*rb[i][1] +
					x[5]*rb[i][2]
			}
		}
	} else {
		// Forward transformation
		for i := 0; i <= 2; i++ {
			xx[i] = x[0]*rb[0][i] +
				x[1]*rb[1][i] +
				x[2]*rb[2][i]

			if (iflag & SEFLG_SPEED) != 0 {
				xx[i+3] = x[3]*rb[0][i] +
					x[4]*rb[1][i] +
					x[5]*rb[2][i]
			}
		}
	}
	*x = xx
}

// ===== 2335 ===== constants for delta t swphlib.c-2335 =============================================================

// DeltaT = Ephemeris Time - Universal Time, in days.
// Before 1955 we use the data developed by Stephenson, Morrison, and Hohenkerk (2016),
// 1955 - today + a couple of years:
// ---------------------------------
// The tabulated values of deltaT from the Astronomical Alamanc (AA 1997 etc. pp. K8-K9) are used. Some more recent
// values have been taken from IERS (http://maia.usno.navy.mil/ser7/deltat.data).
// Bessel's interpolation formula is implemented to obtain fourth order interpolated values at intermediate times.
// The values are adjusted depending on the ephemeris used and its inherent value of secular tidal acceleration ndot.
// future:
// ---------------------------------
// For the time after the last tabulated value, we use the formula of Stephenson (1997; p. 507), with a modification
// that avoids a jump at the end of the tabulated period. A linear term is added that makes a slow transition from the
// table to the formula over a period of 100 years. (Need not be updated, when table will be enlarged.)
// References:
// Stephenson, F. R., and L. V. Morrison, "Long-term changes in the rotation of the Earth: 700 B.C. to A.D. 1980,"
// Philosophical Transactions of the Royal Society of London Series A 313, 47-70 (1984)
// Borkowski, K. M., "ELP2000-85 and the Dynamical Time - Universal Time relation," Astronomy and Astrophysics
// 205, L8-L10 (1988)
// Borkowski's formula is derived from partly doubtful eclipses going back to 2137 BC and uses lunar position based on
// tidal coefficient of -23.9 arcsec/cy^2.
// Chapront-Touze, Michelle, and Jean Chapront, _Lunar Tables and Programs from 4000 B.C. to A.D. 8000_, Willmann-Bell
// 1991 Their table agrees with the one here, but the entries are rounded to the nearest whole second.
// Stephenson, F. R., and M. A. Houlden, _Atlas of Historical Eclipse Maps_, Cambridge U. Press (1986)
// Stephenson, F.R. & Morrison, L.V., "Long-Term Fluctuations in the Earth's Rotation: 700 BC to AD 1990", Philosophical
// Transactions of the Royal Society of London, Ser. A, 351 (1995), 165-202.
// Stephenson, F. Richard, _Historical Eclipses and Earth's Rotation_, Cambridge U. Press (1997)
// Morrison, L. V., and F.R. Stephenson, "Historical Values of the Earth's Clock Error DT and the Calculation of
// Eclipses", JHA xxxv (2004), pp.327-336
// Stephenson, F.R., Morrison, L.V., and Hohenkerk, C.Y., "Measurement of the Earth's Rotation: 720 BC to AD 2015",
// Royal Society Proceedings A  7 Dec 2016, http://rspa.royalsocietypublishing.org/lookup/doi/10.1098/rspa.2016.0404
// Table from AA for 1620 through today Note, Stephenson and Morrison's table starts at the year 1630. The Chapronts'
// table does not agree with the Almanac prior to 1630. The actual accuracy decreases rapidly prior to 1780.
// Jean Meeus, Astronomical Algorithms, 2nd edition, 1998.
// For a comprehensive collection of publications and formulae, see:
// http://www.phys.uu.nl/~vgent/deltat/deltat_modern.htm
// http://www.phys.uu.nl/~vgent/deltat/deltat_old.htm
// For future values of delta t, the following data from the Earth Orientation Department of the US Naval Observatory
// can be used:
// (TAI-UTC) from: ftp://maia.usno.navy.mil/ser7/tai-utc.dat
// (UT1-UTC) from: ftp://maia.usno.navy.mil/ser7/finals.all (cols. 59-68)
//             or: ftp://ftp.iers.org/products/eop/rapid/standard/finals.data
// file description in: ftp://maia.usno.navy.mil/ser7/readme.finals
// Delta T = TAI-UT1 + 32.184 sec = (TAI-UTC) - (UT1-UTC) + 32.184 sec
// Also, there is the following file: http://maia.usno.navy.mil/ser7/deltat.data, but it is about 3 months
// behind (on 3 feb 2009); and predictions: http://maia.usno.navy.mil/ser7/deltat.preds
// Last update of table dt[]: Dieter Koch, 18 dec 2013.
// ATTENTION: Whenever updating this table, do not forget to adjust the macros TABEND and TABSIZ !

const (
	TABSTART = 1620
	TABEND   = 2028
	TABSIZ   = TABEND - TABSTART + 1
	// we make the table greater for additional values read from external file
	TABSIZ_SPACE = TABSIZ + 100

	TAB2_SIZ   = 27
	TAB2_START = -1000
	TAB2_END   = 1600
	TAB2_STEP  = 100

	LTERM_EQUATION_YSTART = 1820
	LTERM_EQUATION_COEFF  = 32

	TAB97_SIZ   = 43
	TAB97_START = -500
	TAB97_END   = 1600
	TAB97_STEP  = 50
)

// dt contains delta T values from 1620 to 2028
var dt = [TABSIZ_SPACE]float64{
	/* 1620.0 - 1659.0 */
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

// dt2 contains delta T values from -1000 through 1600, from Morrison & Stephenson (2004)
var dt2 = [TAB2_SIZ]int16{
	// -1000  -900  -800  -700  -600  -500  -400  -300  -200  -100
	25400, 23700, 22000, 21000, 19040, 17190, 15530, 14080, 12790, 11640,
	// 0   100   200   300   400   500   600   700   800   900
	10580, 9600, 8640, 7680, 6700, 5710, 4740, 3810, 2960, 2200,
	// 1000  1100  1200  1300  1400  1500  1600
	1570, 1090, 740, 490, 320, 200, 120,
}

// dt97 contains delta T values from -500 through 1600, from Stephenson & Morrison (1995)
var dt97 = [TAB97_SIZ]int16{
	// -500  -450  -400  -350  -300  -250  -200  -150  -100   -50
	16800, 16000, 15300, 14600, 14000, 13400, 12800, 12200, 11600, 11100,
	// 0    50   100   150   200   250   300   350   400   450
	10600, 10100, 9600, 9100, 8600, 8200, 7700, 7200, 6700, 6200,
	// 500   550   600   650   700   750   800   850   900   950
	5700, 5200, 4700, 4300, 3800, 3400, 3000, 2600, 2200, 1900,
	// 1000  1050  1100  1150  1200  1250  1300  1350  1400  1450
	1600, 1350, 1100, 900, 750, 600, 470, 380, 300, 230,
	// 1500  1550  1600
	180, 140, 110,
}

// ===== 2536 ===== calc_deltat swephlib.c-2536 ======================================================================

// calcDeltat calculates the delta T value for a given Julian day.
// returns DeltaT (ET - UT) in days
// double tjd 	= 	julian day in UT
// delta t is adjusted to the tidal acceleration that is compatible with the ephemeris flag contained in iflag and with
// the ephemeris files made accessible through swe_set_ephe_path() or swe_set_jplfile().
// If iflag = -1, then the default tidal acceleration is ussed (i.e. that of DE431).
func calcDeltat(tjd float64, iflag int32) (float64, int32, error) {
	// TODO add original comments
	var err error
	var ans float64 = 0
	var B, Y, Ygreg, dd float64
	var iy int
	var retc int32
	var deltaT float64

	// Assuming these are defined elsewhere as global variables/constants
	var deltatModel int32 = swed.AstroModels[SE_MODEL_DELTAT]
	var tidAcc float64
	var denum int32

	if deltatModel == 0 {
		deltatModel = SEMOD_DELTAT_DEFAULT
	}

	epheflag := iflag & SEFLG_EPHMASK
	otherflag := iflag & ^SEFLG_EPHMASK

	// with iflag == -1, we use default tid_acc
	if iflag == -1 {
		retc, _, tidAcc = swiGetTidAcc(tjd, 0, 9999) // for default tid_acc

	} else {
		// otherwise we use tid_acc consistent with epheflag
		//Port: removed code for JPL
		denum = swed.Fidat[SEI_FILE_MOON].SwephDenum

		if swiInitSwedIfStart() == 1 && epheflag == 0 { // TODO check this statement, included originally JPL
			err = errors.New("call sweSetEphePath() before calling sweDeltatEx()")
			retc, err = swiSetTidAcc(tjd, epheflag, denum) // _set_ saves tid_acc in swed
		} else {
			retc, err = swiSetTidAcc(tjd, epheflag, denum) // _set_ saves tid_acc in swed
		}
		tidAcc = swed.TidAcc
	}

	iflag = otherflag | retc
	Y = 2000.0 + (tjd-J2000)/365.25
	Ygreg = 2000.0 + (tjd-J2000)/365.2425

	if deltatModel == SEMOD_DELTAT_STEPHENSON_ETC_2016 && tjd < 2435108.5 {
		deltaT = deltaTStephensonEtc2016(tjd, tidAcc)
		if tjd >= 2434108.5 {
			deltaT += (1.0 - (2435108.5-tjd)/1000.0) * 0.6610218 / 86400.0
		}
		return deltaT, iflag, err
	}

	if deltatModel == SEMOD_DELTAT_ESPENAK_MEEUS_2006 && tjd < 2317746.13090277789 {
		deltaT = deltaTEspenakMeeus1620(tjd, tidAcc)
		return deltaT, iflag, err
	}

	if deltatModel == SEMOD_DELTAT_STEPHENSON_MORRISON_2004 && Y < TABSTART {
		// before 1600:
		if Y < TAB2_END {
			deltaT = deltaTStephensonMorrison2004_1600(tjd, tidAcc)
			return deltaT, iflag, err
		} else {
			// between 1600 and 1620:
			// linear interpolation between end of table dt2 and start of table dt
			if Y >= TAB2_END {
				B = TABSTART - TAB2_END
				iy = int((TAB2_END - TAB2_START) / TAB2_STEP)
				dd = (Y - TAB2_END) / B
				ans = float64(dt2[iy]) + (float64(dt[0])-float64(dt2[iy]))*dd
				ans = adjustForTidacc(ans, Ygreg, tidAcc, SE_TIDAL_26, false)
				deltaT = ans / 86400.0
				return deltaT, iflag, err
			}
		}
	}

	if deltatModel == SEMOD_DELTAT_STEPHENSON_1997 && Y < TABSTART {
		// before 1600:
		if Y < TAB97_END {
			deltaT = deltaTStephensonMorrison1997_1600(tjd, tidAcc)
			return deltaT, iflag, err
		} else {
			// between 1600 and 1620:
			// linear interpolation between end of table dt97 and start of table dt
			if Y >= TAB97_END {
				B = TABSTART - TAB97_END
				iy = int((TAB97_END - TAB97_START) / TAB97_STEP)
				dd = (Y - TAB97_END) / B
				ans = float64(dt97[iy]) + dd*(float64(dt[0])-float64(dt97[iy]))
				ans = adjustForTidacc(ans, Ygreg, tidAcc, SE_TIDAL_26, false)
				deltaT = ans / 86400.0
				return deltaT, iflag, err
			}
		}
	}

	if deltatModel == SEMOD_DELTAT_STEPHENSON_MORRISON_1984 && Y < TABSTART {
		if Y >= 948.0 {
			// Stephenson and Morrison, stated domain is 948 to 1600:
			// 25.5(centuries from 1800)^2 - 1.9159(centuries from 1955)^2
			B = 0.01 * (Y - 2000.0)
			ans = (23.58*B+100.3)*B + 101.6
		} else {
			// Borkowski, before 948 and between 1600 and 1620
			B = 0.01*(Y-2000.0) + 3.75
			ans = 35.0*B*B + 40.0
		}
		deltaT = ans / 86400.0
		return deltaT, iflag, err
	}

	if Y >= TABSTART {
		deltaT = deltaTAa(tjd, tidAcc)
		return deltaT, iflag, err
	}

	deltaT = ans / 86400.0
	return deltaT, iflag, err
}

// ===== 2701 ===== swe_deltat_ex swephlib.c-2701 ====================================================================

func sweDeltatEx(tjd float64, iflag int32) (float64, error) {
	var deltat float64
	var err error
	if swed.DeltaTUserdefIsSet {
		return swed.DeltaTUserdef, nil
	}
	//func calcDeltat(tjd float64, iflag int32, err error) (float64, int32)
	// Assuming calcDeltat is defined elsewhere
	deltat, _, err = calcDeltat(tjd, iflag)
	return deltat, err
}

// ===== 2718 ===== deltat_aa swephlib.c-2718 ========================================================================
// The tabulated values of deltaT, in hundredths of a second, were taken from The Astronomical Almanac 1997etc.,
// pp. K8-K9. Some more recent values are taken from IERS http://maia.usno.navy.mil/ser7/deltat.data .
// Bessel's interpolation formula is implemented to obtain fourth order interpolated values at intermediate times.
// The values are adjusted depending on the ephemeris used and its inherent value of secular tidal acceleration ndot.
// Note by Dieter Jan. 2017: Bessel interpolation assumes equidistant sampling points. However the sampling points are
// not equidistant, because they are for first January of every year and years can have either 365 or 366 days. The
// interpolation uses a step width of 365.25 days. As a consequence, in three out of four years the interpolation does
// not reproduce the exact values of the sampling points on the days they refer to.

func deltaTAa(tjd, tidAcc float64) float64 {
	var ans, ans2, ans3 float64
	var p, B, B2, Y, dd float64
	d := make([]float64, 6)

	// Port: anonymous function to replace GOTO statements
	done := func() float64 {
		ans = adjustForTidacc(ans, Y, tidAcc, SE_TIDAL_26, false)
		return ans / 86400.0
	}

	// read additional values from swedelta.txt
	tabsiz := initDt() // This function needs to be implemented
	tabend := TABSTART + tabsiz - 1
	deltaModel := swed.AstroModels[SE_MODEL_DELTAT]
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
			return done() // No data, can't go on
		}
		// The fraction of tabulation interval
		p = Y - p
		// First order interpolated value
		ans += p * (dt[k] - dt[iy])
		if (iy-1 < 0) || (iy+2 >= tabsiz) {
			return done() // can't do second differences
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
			return done()
		}
		// Compute third differences
		for i := 0; i < 3; i++ {
			d[i] = d[i+1] - d[i]
		}
		B = 2.0 * B / 3.0
		ans += (p - 0.5) * B * d[1]
		if (iy-2 < 0) || (iy+3 > tabsiz) {
			return done()
		}
		// Compute fourth differences
		for i := 0; i < 2; i++ {
			d[i] = d[i+1] - d[i]
		}
		B = 0.125 * B * (p + 1.0) * (p - 2.0)
		ans += B * (d[0] + d[1])
	}
	// today - future: 3rd degree polynomial based on data given by Stephenson/Morrison/Hohenkerk 2016 here:
	// http://astro.ukho.gov.uk/nao/lvm/
	if deltaModel == SEMOD_DELTAT_STEPHENSON_ETC_2016 {
		B = Y - 2000
		if Y < 2500 {
			ans = B*B*B*121.0/30000000.0 + B*B/1250.0 + B*521.0/3000.0 + 64.0
			// for slow transition from tabulated data
			B2 = float64(tabend - 2000)
			ans2 = B2*B2*B2*121.0/30000000.0 + B2*B2/1250.0 + B2*521.0/3000.0 + 64.0
		} else {
			// we use a parable after 2500
			B = 0.01 * (Y - 2000)
			ans = B*B*32.5 + 42.5
		}
		// Formula Stephenson (1997; p. 507), with modification to avoid jump at end of AA table, similar to what
		// Meeus 1998 had suggested.
		// Slow transition within 100 years.
	} else {
		B = 0.01 * (Y - 1820)
		ans = -20 + 31*B*B
		// for slow transition from tabulated data
		B2 = 0.01 * (float64(tabend) - 1820)
		ans2 = -20 + 31*B2*B2
	}

	// slow transition from tabulated values to Stephenson formula
	if Y <= float64(tabend+100) {
		ans3 = dt[tabsiz-1]
		dd = (ans2 - ans3)
		ans += dd * (Y - float64(tabend+100)) * 0.01
	}

	return ans / 86400.0
}

// ===== 2841 ===== deltat_longterm_morrison_stephenson swephlib.c-2841 ==============================================

// deltaTLongTermMorrisonStephenson calculates the long-term delta T according to Morrison and Stephenson.
// Parameter: tjd: Julian Day Number
// Returns: Delta T value in seconds
func deltaTLongTermMorrisonStephenson(tjd float64) float64 {
	ygreg := 2000.0 + (tjd-J2000)/365.2425
	u := (ygreg - 1820) / 100.0
	return -20 + 32*u*u
}

// ===== 2848 ===== deltat_stephenson_morrison_1997_1600 swephlib.c-2848 =============================================

func deltaTStephensonMorrison1997_1600(tjd, tidAcc float64) float64 {
	var ans, ans2, ans3 float64
	var p, B, Y, dd float64

	Y = 2000.0 + (tjd-J2000)/365.25

	// before -500:
	// formula by Stephenson (1997; p. 508) but adjusted to fit the starting point of table dt97 (Stephenson 1997).
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

// ===== 2889 ===== deltat_stephenson_morrison_2004_1600 swephlib.c-2889 =============================================
// Stephenson & Morrison (2004)

func deltaTStephensonMorrison2004_1600(tjd, tidAcc float64) float64 {
	var ans, ans2, ans3 float64
	var B, dd float64
	var tjd0 float64

	// Calculate decimal year
	Y := 2000.0 + (tjd-J2000)/365.2425

	// Before -1000: formula by Stephenson & Morrison (2004; p. 335) but adjusted to fit the starting point of table dt2.
	if Y < TAB2_START {
		ans = deltaTLongTermMorrisonStephenson(tjd)
		ans = adjustForTidacc(ans, Y, tidAcc, SE_TIDAL_26, false)
		// Transition from formula to table over 100 years
		if Y >= TAB2_START-100 {
			// Starting value of table dt2
			ans2 = adjustForTidacc(float64(dt2[0]), TAB2_START, tidAcc, SE_TIDAL_26, false)

			// Value of formula at epoch TAB2_START
			tjd0 = (TAB2_START-2000)*365.2425 + J2000
			ans3 = deltaTLongTermMorrisonStephenson(tjd0)
			ans3 = adjustForTidacc(ans3, Y, tidAcc, SE_TIDAL_26, false)
			dd = ans3 - ans2
			B = (Y - (TAB2_START - 100)) * 0.01
			// Fit to starting point of table dt2
			ans = ans - dd*B
		}
	}
	// Between -1000 and 1600: linear interpolation between values of table dt2 (Stephenson & Morrison 2004)
	if Y >= TAB2_START && Y < TAB2_END {
		Yjul := 2000.0 + (tjd-2451557.5)/365.25
		p := math.Floor(Yjul)
		iy := int((p - TAB2_START) / TAB2_STEP)
		dd = (Yjul - (TAB2_START + TAB2_STEP*float64(iy))) / TAB2_STEP
		ans = float64(dt2[iy]) + float64(dt2[iy+1]-dt2[iy])*dd

		// Correction for tidal acceleration
		ans = adjustForTidacc(ans, Y, tidAcc, SE_TIDAL_26, false)
		ans /= 86400.0
	}
	return ans
}

// ===== 2934 ===== constants for delta t swephlib.c line 2934 =======================================================

const NDTCF16 = 54

// dtcf16 coefficients represent the spline approximation discussed in the paper "Measurement of the Earth's Rotation:
// 720 BC to AD 2015", Stephenson, F.R., Morrison, L.V., and Hohenkerk, C.Y., published by Royal Society Proceedings A
// and available from their website at http://rspa.royalsocietypublishing.org/lookup/doi/10.1098/rspa.2016.0404.
// Year numbers have been replaced by Julian day numbers by D. Koch.
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

// deltat_stephenson_etc_2016 swephlib.c-2999

// deltaTStephensonEtc2016 calculates delta T according to Stephenson et al. 2016
// tjd is the Julian Day Number
// tidAcc is the tidal acceleration value
func deltaTStephensonEtc2016(tjd, tidAcc float64) float64 {
	var dt, t float64
	Ygreg := 2000.0 + (tjd-J2000)/365.2425
	irec := -1

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
		dt = dtcf16[irec][2] +
			dtcf16[irec][3]*t +
			dtcf16[irec][4]*t*t +
			dtcf16[irec][5]*t*t*t
	} else if Ygreg < -720 {
		// for earlier epochs, use long term parabola
		t = (Ygreg - 1825) / 100.0
		dt = -320 + 32.5*t*t
		dt -= 179.7337208 // to make curve continuous on 1 Jan -720 (D. Koch)
	} else {
		// future
		t = (Ygreg - 1825) / 100.0
		dt = -320 + 32.5*t*t
		dt += 269.4790417 // to make curve continuous on 1 Jan 2016 (D. Koch)
	}

	// The parameter adjustAfter1955 must be TRUE here, because the Stephenson 2016 curve is based on occultation data
	// alone, not on IERS data.
	// Note, however, the current function DeltaTStephensonEtc2016() is called only for dates before 1 Jan 1955.
	dt = adjustForTidacc(dt, Ygreg, tidAcc, SE_TIDAL_STEPHENSON_2016, true)
	dt /= 86400.0

	return dt
}

// deltat_espenak_meeus_1620 swephlib.c line 3036

// deltaTEspenakMeeus1620 calculates delta T according to Espenak and Meeus (2006) for the period between years 1620
// and 2000.
// Parameters:
//
//	tjd: Julian Day Number
//	tidAcc: Tidal acceleration value
//
// Returns: Delta T value in days
func deltaTEspenakMeeus1620(tjd, tidAcc float64) float64 {
	var ans float64
	ygreg := 2000.0 + (tjd-J2000)/365.2425
	var u float64

	switch {
	case ygreg < -500:
		ans = deltaTLongTermMorrisonStephenson(tjd)
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

// adjust_for_tidacc swephlib.c-3133

// adjustForTidacc adjusts astronomical almanac entries for tidal acceleration.
// Astronomical Almanac table is corrected by adding the expression -0.000091 (ndot + 26)(year-1955)^2  seconds to
// entries prior to 1955 (AA page K8), where ndot is the secular tidal term in the mean motion of the Moon.
// Entries after 1955 are referred to atomic time standards and are not affected by errors in Lunar or planetary theory.
func adjustForTidacc(ans, y, tidAcc, tidAcc0 float64, adjustAfter1955 bool) float64 {
	if y < 1955.0 || adjustAfter1955 {
		b := y - 1955.0
		ans += -0.000091 * (tidAcc - tidAcc0) * b * b
	}
	return ans
}

// swe_set_tid_acc swephlib.c-3157
// SetTidAcc sets the tidal acceleration value
func sweSetTidAcc(tAcc float64) {
	if tAcc == SE_TIDAL_AUTOMATIC {
		swed.TidAcc = SE_TIDAL_DEFAULT
		swed.IsTidAccManual = false
		return
	}
	swed.TidAcc = tAcc
	swed.IsTidAccManual = true
}

// ===== 3187 ===== init_dt swephlib.c-3187 ==========================================================================

// Read delta t values from external file.
// record structure: year(whitespace)delta_t in 0.01 sec.

func initDt() int {

	if !swed.InitDtDone {
		swed.InitDtDone = true
		// no error message if file is missing
		fp, err := os.Open(filepath.Join(swed.EphePath, "swe_deltat.txt"))
		if err != nil {
			fp, err = os.Open(filepath.Join(swed.EphePath, "sedeltat.txt"))
			if err != nil {
				return TABSIZ
			}
		}
		defer fp.Close()
		scanner := bufio.NewScanner(fp)
		for scanner.Scan() {
			line := scanner.Text()
			if len(line) == 0 {
				continue
			}
			sp := strings.TrimLeft(line, " \t")
			if len(sp) == 0 || sp[0] == '#' || sp[0] == '\n' {
				continue
			}
			year, err := strconv.Atoi(sp[:4])
			if err != nil {
				continue
			}
			tabIndex := year - TABSTART
			// table space is limited. no error msg, if exceeded
			if tabIndex >= TABSIZ_SPACE {
				continue
			}
			// Skip to the value part
			sp = strings.TrimLeft(sp[4:], " \t")
			val, err := strconv.ParseFloat(sp, 64)
			if err != nil {
				continue
			}
			dt[tabIndex] = val
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

// swi_get_tid_acc swephlib.c-3196
func swiGetTidAcc(tjdUt float64, iflag int32, denum int32) (int32, int32, float64) {
	var tidAcc float64
	var denumRet int32

	iflag &= SEFLG_EPHMASK

	if swed.IsTidAccManual {
		tidAcc = swed.TidAcc
		return iflag, denumRet, tidAcc
	}

	// Port: removed handling of JPL and Moshier
	if denum == 0 {
		// Port: removed handling of JPL and Moshier

		// SEFLG_SWIEPH wanted
		if (iflag & SEFLG_SWIEPH) != 0 {
			if swed.Fidat[SEI_FILE_MOON].Fptr != nil {
				denum = swed.Fidat[SEI_FILE_MOON].SwephDenum
			}
		}
	}

	switch denum {
	case 200:
		tidAcc = SE_TIDAL_DE200
	case 403:
		tidAcc = SE_TIDAL_DE403
	case 404:
		tidAcc = SE_TIDAL_DE404
	case 405:
		tidAcc = SE_TIDAL_DE405
	case 406:
		tidAcc = SE_TIDAL_DE406
	case 421:
		tidAcc = SE_TIDAL_DE421
	case 422:
		tidAcc = SE_TIDAL_DE422
	case 430:
		tidAcc = SE_TIDAL_DE430
	case 431:
		tidAcc = SE_TIDAL_DE431
	case 440:
		tidAcc = SE_TIDAL_DE441
	case 441:
		tidAcc = SE_TIDAL_DE441
	default:
		denum = SE_DE_NUMBER
		tidAcc = SE_TIDAL_DEFAULT
	}

	denumRet = denum
	iflag &= SEFLG_EPHMASK

	return iflag, denumRet, tidAcc
}

// swi_set_tid_acc swephlib.c-3240

func swiSetTidAcc(tjdUt float64, iflag int32, denum int32) (int32, error) {
	retc := iflag

	// manual tid_acc overrides automatic tid_acc
	if swed.IsTidAccManual {
		return retc, nil
	}
	var err error
	retc, _, tidAcc := swiGetTidAcc(tjdUt, iflag, denum)
	if err != nil {
		return retc, err
	}
	// Update the global tidAcc
	swed.TidAcc = tidAcc
	return retc, nil
}

// swi_gen_filename swephlib.c-3594

// swiGenFilename generates name of ephemeris file
// * file name looks as follows:
// * swephpl.m30, where
// * "sweph"              	"swiss ephemeris"
// * "pl","mo","as"               planet, moon, or asteroid
// * "m"  or "_"                  BC or AD
// * "30"                         start century
// * tjd        	= ephemeris file for which julian day
// * ipli       	= number of planet
// * fname      	= ephemeris file name
func swiGenFilename(tjd float64, ipli int) string {
	var fname strings.Builder
	switch ipli {
	case SEI_MOON:
		fname.WriteString("semo")
	case SEI_EMB, SEI_MERCURY, SEI_VENUS, SEI_MARS, SEI_JUPITER,
		SEI_SATURN, SEI_URANUS, SEI_NEPTUNE, SEI_PLUTO, SEI_SUNBARY:
		fname.WriteString("sepl")
	case SEI_CERES, SEI_PALLAS, SEI_JUNO, SEI_VESTA, SEI_CHIRON, SEI_PHOLUS:
		fname.WriteString("seas")
	default: // asteroid or planetary moon
		if ipli > SE_PLMOON_OFFSET && ipli < SE_AST_OFFSET {
			return fmt.Sprintf("sat%ssepm%d.%s", DIR_GLUE, ipli, SE_FILE_SUFFIX)
		} else {
			sform := "ast%d%sse%05d.%s"
			if ipli-SE_AST_OFFSET > 99999 {
				sform = "ast%d%ss%06d.%s"
			}
			return fmt.Sprintf(sform, // asteroids or planetary moons: only one file 3000 bc - 3000 ad
				(ipli-SE_AST_OFFSET)/1000,
				DIR_GLUE,
				ipli-SE_AST_OFFSET,
				SE_FILE_SUFFIX)
		}
	}

	// century of tjd
	// if tjd > 1600 then gregorian calendar
	var (
		jyear    int
		gregflag bool
	)

	/* century of tjd */
	/* if tjd > 1600 then gregorian calendar */
	if tjd >= 2305447.5 {
		gregflag = true
	}

	// Call the calendar conversion function
	cal := SE_GREG_CAL
	if !gregflag {
		cal = SE_JUL_CAL
	}
	jyear, _, _, _ = SweRevJul(tjd, cal)
	// start century of file containing tjd
	var sgn int
	if jyear < 0 {
		sgn = -1
	} else {
		sgn = 1
	}
	icty := jyear / 100
	if sgn < 0 && jyear%100 != 0 {
		icty--
	}
	for icty%NCTIES != 0 {
		icty--
	}
	// Append B.C. or A.D. indicator
	if icty < 0 {
		fname.WriteString("m")
	} else {
		fname.WriteString("_")
	}
	// Append the century and file suffix
	fname.WriteString(fmt.Sprintf("%02d.%s", abs(icty), SE_FILE_SUFFIX))
	return fname.String()
}

// Helper function for absolute value of int in SwiGenFilename
func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

// ===== 3691 ===== swi_cutstr swephlib.c-3691 =======================================================================
// cut the string s at any char in cutlist; put pointers to partial strings into cpos[0..n-1], return number of partial
// strings; if less than nmax fields are found, the first empty pointer is set to NULL.
// More than one character of cutlist in direct sequence count as one separator only! cut_str_any("word,,,word2",","..)
// cuts only two parts, cpos[0] = "word" and cpos[1] = "word2".
// If more than nmax fields are found, nmax is returned and the last field nmax-1 rmains un-cut.

func swiCutstr(s string, cutlist string, cpos []string, nmax int) int {
	bytes := []byte(s)
	n := 1
	if len(bytes) == 0 {
		if nmax > 0 {
			cpos[0] = ""
		}
		return n
	}
	cpos[0] = string(bytes)
	for i := 0; i < len(bytes); i++ {
		if strings.ContainsRune(cutlist, rune(bytes[i])) && n < nmax {
			bytes[i] = 0
			for i+1 < len(bytes) && strings.ContainsRune(cutlist, rune(bytes[i+1])) {
				i++
			}
			if i+1 < len(bytes) {
				cpos[n] = string(bytes[i+1:])
				n++
			}
		}
		if bytes[i] == '\n' || bytes[i] == '\r' { // treat nl or cr like end of string
			bytes[i] = 0
			break
		}
	}
	for i := n; i < nmax; i++ {
		cpos[i] = ""
	}
	return n
}

// ===== 3722 ===== *swi_right_trim swephlib.c-3722 ==================================================================

// rightTrim removes trailing whitespace from a string.
// This version uses the standard library and is more idiomatic Go.
// Port: Changed from original version.
func rightTrim(s string) string {
	return strings.TrimRight(s, " \t\n\r\v\f")
}
