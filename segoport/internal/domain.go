package internal

import "os"

var swed SweData

// swe_data sweph.h-0790
// if this is changed, then also update initialisation in sweph.c
type SweData struct {
	EphePathIsSet bool
	JplFileIsOpen bool
	FixFp         *os.File // Fixed stars file pointer
	EphePath      string
	// Port: ignored JplFnam            [AS_MAXCH]byte
	// Port: ignored JplDenum           int32
	LastEpheFlag       int32
	GeoposIsSet        bool
	AyanaIsSet         bool
	IsOldStarfile      bool
	EopTjdBeg          float64
	EopTjdBegHorizons  float64
	EopTjdEnd          float64
	EopTjdEndAdd       float64
	EopDpsiLoaded      int
	TidAcc             float64
	IsTidAccManual     bool
	InitDtDone         bool
	SwedIsInitialised  bool
	DeltaTUserdefIsSet bool
	DeltaTUserdef      float64
	AstG               float64
	AstH               float64
	AstDiam            float64
	Astelem            [AS_MAXCH * 10]byte
	ISavedPlanetName   int
	SavedPlanetName    string
	Dpsi               []float64
	Deps               []float64
	Timeout            int32
	AstroModels        []int32 // Port: changed array[SEI_NMODELS] into a slice
	DoInterpolateNut   bool
	Interpol           Interpol
	Fidat              [SEI_NEPHFILES]FileData
	Gcdat              GenConst
	Pldat              [SEI_NPLANETS]PlanData
	Nddat              [SEI_NNODE_ETC]PlanData
	Savedat            [SE_NPLANETS + 1]SavePositions
	Oec                Epsilon
	Oec2000            Epsilon
	Nut                Nut
	Nut2000            Nut
	Nutv               Nut
	Topd               TopoData
	Sidd               SidData
	NFixstarsReal      bool // real number of fixed stars in sefstars.txt
	NFixstarsNamed     bool // number of fixed stars with tradtional name
	NFixstarsRecords   bool // number of fixed stars records in fixed_stars
	FixedStars         []FixedStar
}

var sweData SweData

// interpol sweph.h-0784
type Interpol struct {
	TjdNut0  float64
	TjdNut2  float64
	NutDpsi0 float64
	NutDpsi1 float64
	NutDpsi2 float64
	NutDeps0 float64
	NutDeps1 float64
	NutDeps2 float64
}

// file_data sweph.h-0708
type FileData struct {
	Fnam       string                   // ephemeris file name
	Fversion   int                      // version number of file
	Astnam     string                   // asteroid name, if asteroid file
	SwephDenum int32                    // DE number of JPL ephemeris, which this file is derived from
	Fptr       *os.File                 // ephemeris file pointer
	Tfstart    float64                  // file may be used from this date
	Tfend      float64                  // through this date
	Iflg       int32                    // byte reorder flag and little/bigendian flag
	Npl        int16                    // how many planets in file
	Ipl        [SEI_FILE_NMAXPLAN]int32 // planet numbers
}

// gen_const sweph.h-0722
// Port: GenConst represents general astronomical constants
type GenConst struct {
	Clight       float64 // Speed of light
	Aunit        float64 // Astronomical unit
	Helgravconst float64 // Heliocentric gravitational constant
	Ratme        float64 // Ratio of mass of Earth to mass of Moon
	Sunradius    float64 // Radius of the Sun
}

// plan_data sweph.h-0610
type PlanData struct {
	// The following data are read from file only once, immediately after file has been opened
	Ibdy int   // internal body number
	Iflg int32 // contains several bit flags describing the data:
	// SEI_FLG_HELIO: true if helio, false if bary
	// SEI_FLG_ROTATE: TRUE if coefficients are referred to coordinate system of orbital plane
	// SEI_FLG_ELLIPSE: TRUE if reference ellipse
	Ncoe int // # of coefficients of ephemeris polynomial, is polynomial order + 1
	// where is the segment index on the file
	Lndx0   int32   // file position of begin of planet's index
	Nndx    int32   // number of index entries on file: computed
	Tfstart float64 // file contains ephemeris for tfstart thru tfend
	Tfend   float64 // for this particular planet !!!
	Dseg    float64 // segment size (days covered by a polynomial)

	// orbital elements:
	Telem float64 // epoch of elements
	Prot  float64
	Qrot  float64
	Dprot float64
	Dqrot float64
	Rmax  float64 // normalisation factor of cheby coefficients
	// in addition, if reference ellipse is used:
	Peri  float64
	Dperi float64
	Refep []float64 // slice of cheby coeffs of reference ellipse, size of data is 2 x ncoe
	// unpacked segment information, only updated when a segment is read:
	Tseg0 float64   // start jd of current segment
	Tseg1 float64   // end jd of current segment
	Segp  []float64 // slice of unpacked cheby coeffs of segment; the size is 3 x ncoe
	Neval int       // how many coefficients to evaluate. this may be less than ncoe
	// result of most recent data evaluation for this body:
	Teval   float64     // time for which previous computation was made
	Iephe   int32       // which ephemeris was used
	X       [6]float64  // position and speed vectors equatorial J2000
	Xflgs   int32       // hel., light-time, aberr., prec. flags etc.
	Xreturn [24]float64 // return positions:
	// xreturn[0:6]   ecliptic polar coordinates
	// xreturn[6:12]  ecliptic cartesian coordinates
	// xreturn[12:18] equatorial polar coordinates
	// xreturn[18:24] equatorial cartesian coordinates
}

// save_positions sweph.h-730
type SavePositions struct {
	Ipl      int
	Tsave    float64
	Iflgsave int32
	// Position at t = tsave,
	// in ecliptic polar (offset 0),
	//    ecliptic cartesian (offset 6),
	//    equatorial polar (offset 12),
	//    and equatorial cartesian coordinates (offset 18).
	// 6 doubles each for position and speed coordinates.
	Xsaves [24]float64
}

// epsilon sweph.h-0600
// obliquity of ecliptic
type Epsilon struct {
	Teps float64 // jd
	Eps  float64 // eps
	Seps float64 // sin(eps)
	Ceps float64 // cos(eps)
}

// nutation sweph.h-0690
// nutation
type Nut struct {
	Tnut   float64   // time
	Nutlo  []float64 // nutation in longitude and obliquity
	Snut   float64   // sine of nutation in obliquity
	Cnut   float64   // cosine of nutation in obliquity
	Matrix [3][3]float64
}

// topo_data sweph.h-0758
type TopoData struct {
	Geolon float64    // geographical longitude
	Geolat float64    // geographical latitude
	Geoalt float64    // geographical altitude
	Teval  float64    // evaluation time
	TjdUt  float64    // Julian Day number, Universal Time
	Xobs   [6]float64 // observer position
}

// sid_data sweph.h-0765
type SidData struct {
	SidMode int32   // sidereal mode
	AyanT0  float64 // ayanamsa at t0
	T0      float64 // reference time
	T0IsUT  bool    // flag indicating if T0 is in Universal Time
}

// fixed_star sweph.h-0773
type FixedStar struct {
	Skey      string // may be prefixed with comma
	StarName  string
	StarBayer string
	StarNo    string
	Epoch     float64
	Ra        float64 // right ascension
	De        float64 // declination
	RaMot     float64 // proper motion in right ascension
	DeMot     float64 // proper motion in declination
	RadVel    float64 // radial velocity
	Parall    float64 // parallax
	Mag       float64 // magnitude
}
