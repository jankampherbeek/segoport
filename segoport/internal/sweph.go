package internal

import (
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"strings"
)

// Port: first items from sweph.h (line nrs preceded with 'h', followed by items from sweph.c =========================

// ===== h0309 ===== dot+_prod (defined function) sweph.h-0309 ========================================================

func dotProd(x, y [3]float64) float64 {
	return x[0]*y[0] + x[1]*y[1] + x[2]*y[2]
}

// free_planets sweph.c-1158
func freePlanets() {
	// Free planets data space
	for i := 0; i < SEI_NPLANETS; i++ {
		swed.Pldat[i].Segp = nil
		swed.Pldat[i].Refep = nil

		// Reset the entire structure to zero values
		swed.Pldat[i] = PlanData{} // Assuming PlanData is the Go equivalent of plan_data
	}
	for i := 0; i <= SE_NPLANETS; i++ { // "<=" is correct! see decl.
		// Reset save positions structure
		swed.Savedat[i] = SavePositions{}
	}
	// Clear node data space
	for i := 0; i < SEI_NNODE_ETC; i++ {
		swed.Nddat[i] = PlanData{}
	}
}

// swi_init_swed_if_start sweph.c-1179
// swiInitSwedIfStart initializes the swed structure if not already initialized.
// Returns 1 if initialization is performed, 0 otherwise.
func swiInitSwedIfStart() int32 {
	if !swed.SwedIsInitialised {
		swed.EphePath = SE_EPHE_PATH
		// Port: skipped JPL file
		sweSetTidAcc(SE_TIDAL_AUTOMATIC)
		swed.SwedIsInitialised = true
		return 1
	}
	return 0
}

// swi_close_keep_topo_etc sweph.c-1195
// swiCloseKeepTopoEtc closes all open files, frees space of planetary data, and deletes memory of all computed
// positions while keeping topocentric data
func swiCloseKeepTopoEtc() {
	// Close SWISSEPH files
	for i := 0; i < SEI_NEPHFILES; i++ {
		if swed.Fidat[i].Fptr != nil {
			swed.Fidat[i].Fptr.Close()
		}
		swed.Fidat[i] = FileData{} // Reset to zero value
	}

	freePlanets() // Using our previously translated function

	// Reset various structures to zero values
	swed.Oec = Epsilon{} // Using our previously defined Epsilon struct
	swed.Oec2000 = Epsilon{}
	swed.Nut = Nut{} // Using our previously defined Nut struct
	swed.Nut2000 = Nut{}
	swed.Nutv = Nut{}

	// Reset astro models
	swed.AstroModels = make([]int32, SEI_NMODELS) // Creates new slice with zero values

	// Close fixed stars
	if swed.FixFp != nil {
		swed.FixFp.Close()
		swed.FixFp = nil
	}

	// Reset other parameters
	sweSetTidAcc(SE_TIDAL_AUTOMATIC)
	swed.IsOldStarfile = false
	swed.ISavedPlanetName = 0
	swed.SavedPlanetName = "" // Assuming this is a string in Go
	swed.Timeout = 0
}

// swe_close sweph.c-1230

// SweClose closes all open files, frees space of planetary data, and deletes memory of all computed positions
func SweClose() {
	// Close SWISSEPH files
	for i := 0; i < SEI_NEPHFILES; i++ {
		if swed.Fidat[i].Fptr != nil {
			swed.Fidat[i].Fptr.Close() // Assuming fptr is *os.File
		}
		swed.Fidat[i] = FileData{} // Reset to zero value
	}

	freePlanets()

	// Reset various structures to zero values
	swed.Oec = Epsilon{}
	swed.Oec2000 = Epsilon{}
	swed.Nut = Nut{}
	swed.Nut2000 = Nut{}
	swed.Nutv = Nut{}
	swed.AstroModels = make([]int32, SEI_NMODELS)

	// Port: removed part that closes JPL files

	// Close fixed stars file
	if swed.FixFp != nil {
		swed.FixFp.Close()
		swed.FixFp = nil
	}

	sweSetTidAcc(SE_TIDAL_AUTOMATIC)
	swed.GeoposIsSet = false
	swed.AyanaIsSet = false
	swed.IsOldStarfile = false
	swed.ISavedPlanetName = 0
	swed.SavedPlanetName = ""
	swed.Topd = TopoData{}
	swed.Sidd = SidData{}
	swed.Timeout = 0
	swed.LastEpheFlag = 0
	swed.Dpsi = nil
	swed.Deps = nil

	if swed.NFixstarsRecords {
		swed.FixedStars = nil
		swed.NFixstarsReal = false
		swed.NFixstarsNamed = false
		swed.NFixstarsRecords = false
	}
}

// ===== 1310 ============ swe_set_ephe_path sweph.c-1310 ===========================================================

// SweSetEphePath sets ephemeris file path, also calls swe_close(). this makes sure that swe_calc() won't return planet
// positions previously computed from other ephemerides
func SweSetEphePath(path string) {
	const maxPathLen = AS_MAXCH - 1 - 13
	// Close all open files and delete all planetary data
	swiCloseKeepTopoEtc()
	swiInitSwedIfStart()
	swed.EphePathIsSet = true
	var s string
	// Environment variable SE_EPHE_PATH has priority
	if envPath := os.Getenv("SE_EPHE_PATH"); envPath != "" && len(envPath) <= maxPathLen {
		s = envPath
	} else if path == "" {
		s = SE_EPHE_PATH
	} else if len(path) <= maxPathLen {
		s = path
	} else {
		s = SE_EPHE_PATH
	}
	if s != "" && !strings.HasSuffix(s, string(os.PathSeparator)) {
		s += string(os.PathSeparator)
	}
	swed.EphePath = s

	// Try to open lunar ephemeris to get DE number and set tidal acceleration
	// TODO enable iflag if SweCalc is enabled again
	//iflag := int32(SEFLG_SWIEPH | SEFLG_J2000 | SEFLG_TRUEPOS | SEFLG_ICRS)
	swed.LastEpheFlag = 2

	// TODO enable xx if SweCalc is enabled again
	// var xx []float64 // TODO check if slice is OK, originally array with length 6

	// TODO enable SweCalc  . Check why this call is necessary.
	//	SweCalc(J2000, SE_MOON, iflag, xx, nil)

	if swed.Fidat[SEI_FILE_MOON].Fptr != nil {
		swiSetTidAcc(0, 0, swed.Fidat[SEI_FILE_MOON].SwephDenum)
	}

}

// ===== 1530 ========== calc_epsilon sweph.c-1530 ===================================================================

// calcEpsilon calculates obliquity of ecliptic and stores it together
// with its date, sine, and cosine
func calcEpsilon(tjd float64, iflag int32, e *Epsilon) {
	e.Teps = tjd
	e.Eps = swiEpsiln(tjd, iflag)
	e.Seps = math.Sin(e.Eps)
	e.Ceps = math.Cos(e.Eps)
}

// ----- 2359 ===========  swi_fopen sweph.c-2359 ====================================================================
// TODO Port: changed this code a lot, do need to check, maybe compare with C# version
func SwiFopen(ifno int, fname, ephepath string) (*os.File, error) {
	// var fnamp string
	var err error
	var cpos []string
	//if ifno >= 0 {
	//	fnamp = swed.Fidat[ifno].Fnam
	//} else {
	//	fnamp = ""
	//}

	// Split path using PATH_SEPARATOR
	//paths := strings.Split(ephepath, PATH_SEPARATOR)
	swiCutstr(ephepath, PATH_SEPARATOR, cpos, 20)
	for _, path := range cpos {
		s := path

		// Handle current directory case
		if s == "." {
			s = ""
		} else if s != "" && !strings.HasSuffix(s, DIR_GLUE) {
			s += DIR_GLUE
		}

		// Check path length
		if len(s)+len(fname) >= AS_MAXCH {
			return nil, fmt.Errorf("error: file path and name must be shorter than %d", AS_MAXCH)
		}
		fullPath := filepath.Join(s, fname)
		// Store the full path if ifno >= 0
		//if ifno >= 0 {
		//	fnamp = fullPath
		//}
		// Try to open the file
		fp, err := os.Open(fullPath)
		if err == nil {
			return fp, nil
		}
	}

	// File not found in any path
	err = fmt.Errorf("SwissEph file '%s' not found in PATH '%s'", fname, ephepath)
	return nil, err
}

// ===== 2406 ============ swi_get_dewnum sweph.c-2406 ==============================================================
// Port: removed a few lines for moshier and jpl

func swiGetDenum(ipli int32, iflag int32) int32 {
	var fdp *FileData
	switch {
	case ipli > SE_AST_OFFSET:
		fdp = &swed.Fidat[SEI_FILE_ANY_AST]
	case ipli > SE_PLMOON_OFFSET:
		fdp = &swed.Fidat[SEI_FILE_ANY_AST]
	case ipli == SEI_CHIRON ||
		ipli == SEI_PHOLUS ||
		ipli == SEI_CERES ||
		ipli == SEI_PALLAS ||
		ipli == SEI_JUNO ||
		ipli == SEI_VESTA:
		fdp = &swed.Fidat[SEI_FILE_MAIN_AST]
	case ipli == SEI_MOON:
		fdp = &swed.Fidat[SEI_FILE_MOON]
	default:
		fdp = &swed.Fidat[SEI_FILE_PLANET]
	}
	if fdp != nil {
		if fdp.SwephDenum != 0 {
			return fdp.SwephDenum
		}
		return SE_DE_NUMBER
	}
	return SE_DE_NUMBER
}

// ===== 2444 ===== calc_center_body sweph.c-2444 ===================================================================

func calcCenterBody(ipli int32, iflag int32, xx *[]float64, xcom *[]float64) error {
	if (iflag & SEFLG_CENTER_BODY) == 0 {
		return nil
	}
	if ipli < SEI_MARS || ipli > SEI_PLUTO {
		return nil
	}
	for i := 0; i <= 5; i++ {
		(*xx)[i] += (*xcom)[i]
	}
	return nil
}

// ===== 3588 ===== swi_nutate sweph.c-3588 =========================================================================

// SwiNutate multiplies cartesian equatorial coordinates with previously calculated nutation matrix. also corrects speed.
func SwiNutate(xx []float64, iflag int32, backward bool) {
	x := make([]float64, 6)
	xv := make([]float64, 6)
	for i := 0; i <= 2; i++ {
		if backward {
			x[i] = xx[0]*swed.Nut.Matrix[i][0] +
				xx[1]*swed.Nut.Matrix[i][1] +
				xx[2]*swed.Nut.Matrix[i][2]
		} else {
			x[i] = xx[0]*swed.Nut.Matrix[0][i] +
				xx[1]*swed.Nut.Matrix[1][i] +
				xx[2]*swed.Nut.Matrix[2][i]
		}
	}
	if (iflag & SEFLG_SPEED) != 0 {
		// Correct speed:
		// First correct rotation
		for i := 0; i <= 2; i++ {
			if backward {
				x[i+3] = xx[3]*swed.Nut.Matrix[i][0] +
					xx[4]*swed.Nut.Matrix[i][1] +
					xx[5]*swed.Nut.Matrix[i][2]
			} else {
				x[i+3] = xx[3]*swed.Nut.Matrix[0][i] +
					xx[4]*swed.Nut.Matrix[1][i] +
					xx[5]*swed.Nut.Matrix[2][i]
			}
		}
		// Then apparent motion due to change of nutation during day.
		// This makes a difference of 0.01"
		for i := 0; i <= 2; i++ {
			if backward {
				xv[i] = xx[0]*swed.Nutv.Matrix[i][0] +
					xx[1]*swed.Nutv.Matrix[i][1] +
					xx[2]*swed.Nutv.Matrix[i][2]
			} else {
				xv[i] = xx[0]*swed.Nutv.Matrix[0][i] +
					xx[1]*swed.Nutv.Matrix[1][i] +
					xx[2]*swed.Nutv.Matrix[2][i]
			}
			// New speed
			xx[i+3] = x[i+3] + (x[i]-xv[i])/NUT_SPEED_INTV
		}
	}
	// New position
	for i := 0; i <= 2; i++ {
		xx[i] = x[i]
	}
}

// ========= 4360 ======== get_new_segment sweph.c-4360 =============================================================
// fetch chebyshew coefficients from sweph file for
// tjd 		time
// ipli		planet number
// ifno		file number
// serr		error string

func getNewSegment(tjd float64, ipli, ifno int) (int, string) {
	var serr string
	var fendianValue int
	pdp := &swed.Pldat[ipli]
	fdp := &swed.Fidat[ifno]
	fp := fdp.Fptr
	freord := (fdp.Iflg & SEI_FILE_REORD) != 0
	fendian := (fdp.Iflg & SEI_FILE_LITENDIAN) != 0
	// Compute segment number
	iseg := int32((tjd - pdp.Tfstart) / pdp.Dseg)
	pdp.Tseg0 = pdp.Tfstart + float64(iseg)*pdp.Dseg
	pdp.Tseg1 = pdp.Tseg0 + pdp.Dseg
	// Get file position of coefficients
	fpos := pdp.Lndx0 + int32(iseg)*3
	var fposBytes [4]byte
	fendianValue = 0
	if fendian {
		fendianValue = 1
	}
	retc, _ := doFread(&fposBytes, 3, 1, 4, fp, fpos, freord, fendianValue, ifno)
	if retc != OK {
		return returnErrorGns(fdp)
	}
	// Seek to position
	if _, err := fp.Seek(int64(fpos), io.SeekStart); err != nil {
		return returnErrorGns(fdp)
	}
	// Allocate space for Chebyshev coefficients
	if pdp.Segp == nil {
		pdp.Segp = make([]float64, pdp.Ncoe*3)
	}
	// Read coefficients for 3 coordinates
	for icoord := 0; icoord < 3; icoord++ {
		idbl := icoord * pdp.Ncoe
		// Read header
		var c [4]byte
		fendianValue = 0
		if fendian {
			fendianValue = 1
		}
		retc, _ = doFread(&c, 1, 2, 1, fp, SEI_CURR_FPOS, freord, fendianValue, ifno)
		if retc != OK {
			return returnErrorGns(fdp)
		}
		// Process sizes
		var nsizes int
		var nsize [6]int
		var nco int
		if c[0]&128 != 0 {
			nsizes = 6
			fendianValue = 0
			if fendian {
				fendianValue = 1
			}
			retc, _ = doFread(c[2:], 1, 2, 1, fp, SEI_CURR_FPOS, freord, fendianValue, ifno)
			if retc != OK {
				return returnErrorGns(fdp)
			}
			nsize[0] = int(c[1]) / 16
			nsize[1] = int(c[1]) % 16
			nsize[2] = int(c[2]) / 16
			nsize[3] = int(c[2]) % 16
			nsize[4] = int(c[3]) / 16
			nsize[5] = int(c[3]) % 16
			nco = nsize[0] + nsize[1] + nsize[2] + nsize[3] + nsize[4] + nsize[5]
		} else {
			nsizes = 4
			nsize[0] = int(c[0]) / 16
			nsize[1] = int(c[0]) % 16
			nsize[2] = int(c[1]) / 16
			nsize[3] = int(c[1]) % 16
			nco = nsize[0] + nsize[1] + nsize[2] + nsize[3]
		}
		// Check coefficient count
		if nco > pdp.Ncoe {
			serr = fmt.Sprintf("error in ephemeris file: %d coefficients instead of %d. ",
				nco, pdp.Ncoe)
			if len(serr)+len(fdp.Fnam) < AS_MAXCH-1 {
				serr = fmt.Sprintf("error in ephemeris file %s: %d coefficients instead of %d. ",
					fdp.Fnam, nco, pdp.Ncoe)
			}
			pdp.Segp = nil
			return ERR, serr
		}
		// Unpack coefficients
		if err := unpackCoefficients(pdp, fdp, fp, nsizes, nsize, idbl, freord, fendian, ifno); err != nil {
			return returnErrorGns(fdp)
		}
	}
	return OK, ""
}

// helper function for GetNewSegments
func unpackCoefficients(pdp *PlanData, fdp *FileData, fp *os.File,
	nsizes int, nsize [6]int, idbl int,
	freord bool, fendian bool, ifno int) error {
	longs := make([]uint32, MAXORD+1)
	for i := 0; i < nsizes; i++ {
		if nsize[i] == 0 {
			continue
		}
		if i < 4 {
			// Handle full byte packing
			j := 4 - i
			k := nsize[i]
			fendianValue := 0
			if fendian {
				fendianValue = 1
			}
			retc, errStr := doFread(&longs, j, k, 4, fp, SEI_CURR_FPOS, freord, fendianValue, ifno)
			if retc != OK {
				return fmt.Errorf(errStr)
			}
			for m := 0; m < k; m++ {
				if longs[m]&1 != 0 { // negative
					pdp.Segp[idbl+m] = -float64((longs[m]+1)/2) / 1e9 * pdp.Rmax / 2
				} else {
					pdp.Segp[idbl+m] = float64(longs[m]/2) / 1e9 * pdp.Rmax / 2
				}
			}
			idbl += k
		} else if i == 4 {
			// Handle half byte packing
			if err := unpackHalfBytes(pdp, fp, longs, nsize[i], idbl, freord, fendian, ifno); err != nil {
				return err
			}
		} else if i == 5 {
			// Handle quarter byte packing
			if err := unpackQuarterBytes(pdp, fp, longs, nsize[i], idbl, freord, fendian, ifno); err != nil {
				return err
			}
		}
	}
	return nil
}

// helper function for GetNewSegments

func returnErrorGns(fdp *FileData) (int, string) {
	if fdp.Fptr != nil {
		fdp.Fptr.Close()
		fdp.Fptr = nil
	}
	freePlanets()
	return ERR, "Error in get_new_segment"
}

// TODO check this  Claude Sonnet suggested this function but it is not clear why, it is a helper function for GetNewSegments
func unpackHalfBytes(pdp *PlanData, fp *os.File, longs []uint32, nsize, idbl int,
	freord bool, fendian bool, ifno int) error {
	// Calculate how many bytes we need to read
	// Each long (4 bytes) contains 2 half-byte values
	k := (nsize + 1) / 2
	// Read the packed values
	fendianValue := 0
	if fendian {
		fendianValue = 1
	}
	retc, errStr := doFread(&longs, 1, k, 4, fp, SEI_CURR_FPOS, freord, fendianValue, ifno)
	if retc != OK {
		return fmt.Errorf(errStr)
	}
	// Current position in the output array
	currIdbl := idbl
	// Process each long value
	for m := 0; m < k && (m*2) < nsize; m++ {
		value := longs[m]
		// Process first half (16s position)
		if (m * 2) < nsize {
			if value&16 != 0 {
				// Negative value
				pdp.Segp[currIdbl] = -float64((value+16)/16/2) * pdp.Rmax / 2 / 1e9
			} else {
				// Positive value
				pdp.Segp[currIdbl] = float64(value/16/2) * pdp.Rmax / 2 / 1e9
			}
			currIdbl++
		}
		// Process second half (1s position)
		if (m*2 + 1) < nsize {
			value = value % 16 // Remove the first half
			if value&1 != 0 {
				// Negative value
				pdp.Segp[currIdbl] = -float64((value+1)/2) * pdp.Rmax / 2 / 1e9
			} else {
				// Positive value
				pdp.Segp[currIdbl] = float64(value/2) * pdp.Rmax / 2 / 1e9
			}
			currIdbl++
		}
	}
	return nil
}

// TODO check this  Claude Sonnet suggested this function but it is not clear why, it is a helper function for GetNewSegments

func unpackQuarterBytes(pdp *PlanData, fp *os.File, longs []uint32, nsize, idbl int,
	freord bool, fendian bool, ifno int) error {
	// Calculate how many bytes we need to read
	// Each long (4 bytes) contains 4 quarter-byte values
	k := (nsize + 3) / 4
	// Read the packed values
	fendianValue := 0
	if fendian {
		fendianValue = 1
	}
	retc, errStr := doFread(&longs, 1, k, 4, fp, SEI_CURR_FPOS, freord, fendianValue, ifno)
	if retc != OK {
		return fmt.Errorf(errStr)
	}
	// Current position in the output array
	currIdbl := idbl
	// Process each long value
	for m := 0; m < k && (m*4) < nsize; m++ {
		value := longs[m]
		// Process four quarters (64s, 16s, 4s, 1s positions)
		for n := 0; n < 4 && (m*4+n) < nsize; n++ {
			switch n {
			case 0:
				// First quarter (64s position)
				if value&64 != 0 {
					pdp.Segp[currIdbl] = -float64((value+64)/64/2) * pdp.Rmax / 2 / 1e9
				} else {
					pdp.Segp[currIdbl] = float64(value/64/2) * pdp.Rmax / 2 / 1e9
				}
				value = value % 64
			case 1:
				// Second quarter (16s position)
				if value&16 != 0 {
					pdp.Segp[currIdbl] = -float64((value+16)/16/2) * pdp.Rmax / 2 / 1e9
				} else {
					pdp.Segp[currIdbl] = float64(value/16/2) * pdp.Rmax / 2 / 1e9
				}
				value = value % 16
			case 2:
				// Third quarter (4s position)
				if value&4 != 0 {
					pdp.Segp[currIdbl] = -float64((value+4)/4/2) * pdp.Rmax / 2 / 1e9
				} else {
					pdp.Segp[currIdbl] = float64(value/4/2) * pdp.Rmax / 2 / 1e9
				}
				value = value % 4
			case 3:
				// Fourth quarter (1s position)
				if value&1 != 0 {
					pdp.Segp[currIdbl] = -float64((value+1)/2) * pdp.Rmax / 2 / 1e9
				} else {
					pdp.Segp[currIdbl] = float64(value/2) * pdp.Rmax / 2 / 1e9
				}
			}
			currIdbl++
		}
	}
	return nil
}

// ===== 4889 ===== do_fread sweph.c-4889 ===========================================================================

// SWISSEPH
// reads from a file and, if necessary, reorders bytes
// targ 	target pointer
// size		size of item to be read
// count	number of items
// corrsize	in what size should it be returned (e.g. 3 byte int -> 4 byte int)
// fp		file pointer
// fpos		file position: if (fpos >= 0) then fseek
// freord	reorder bytes or no
// fendian	little/bigendian
// ifno		file number
// serr		error string

func doFread(trg interface{}, size, count, corrsize int, fp *os.File, fpos int32,
	freord bool, fendian int, ifno int) (int, string) {
	totsize := size * count
	// Seek to position if specified
	if fpos >= 0 {
		if _, err := fp.Seek(int64(fpos), io.SeekStart); err != nil {
			return ERR, "Failed to seek in file"
		}
	}
	// Create buffer for reading
	space := make([]byte, totsize)
	targ := make([]byte, count*corrsize)
	// if no byte reorder has to be done, and read size == return size
	if !freord && size == corrsize {
		n, err := fp.Read(targ)
		if err != nil || n != totsize {
			serr := "Ephemeris file is damaged (1). "
			if len(serr)+len(swed.Fidat[ifno].Fnam) < AS_MAXCH-1 {
				serr = "Ephemeris file " + swed.Fidat[ifno].Fnam + " is damaged (2)."
			}
			return ERR, serr
		}
		// Copy data to target interface
		copyToInterface(trg, targ)
		return OK, ""
	}
	// Read into temporary buffer
	n, err := fp.Read(space)
	if err != nil || n != totsize {
		serr := "Ephemeris file is damaged (3). "
		if len(serr)+len(swed.Fidat[ifno].Fnam) < AS_MAXCH-1 {
			serr = "Ephemeris file " + swed.Fidat[ifno].Fnam + " is damaged (4)."
		}
		return ERR, serr
	}
	// Clear target if sizes don't match
	if size != corrsize {
		for i := range targ {
			targ[i] = 0
		}
	}
	// Reorder bytes as needed
	for i := 0; i < count; i++ {
		for j := size - 1; j >= 0; j-- {
			var k int
			if freord {
				k = size - j - 1
			} else {
				k = j
			}

			if size != corrsize {
				if (fendian == SEI_FILE_BIGENDIAN && !freord) ||
					(fendian == SEI_FILE_LITENDIAN && freord) {
					k += corrsize - size
				}
			}

			targ[i*corrsize+k] = space[i*size+j]
		}
	}
	// Copy data to target interface
	copyToInterface(trg, targ)
	return OK, ""
}

// Helper function for doFread to copy bytes to interface
func copyToInterface(dst interface{}, src []byte) {
	switch v := dst.(type) {
	case *[]byte:
		copy(*v, src)
	case *[]int32:
		data := make([]int32, len(src)/4)
		for i := range data {
			data[i] = int32(binary.LittleEndian.Uint32(src[i*4:]))
		}
		*v = data
	case *[]float64:
		data := make([]float64, len(src)/8)
		for i := range data {
			bits := binary.LittleEndian.Uint64(src[i*8:])
			data[i] = float64(bits)
		}
		*v = data
		// Add more types as needed   // TODO check, I do not see more data
	}
}

// ===== 4956 ===== rot_back sweph.c 4956 ===========================================================================
// SWISSEPH
// adds reference orbit to chebyshew series (if SEI_FLG_ELLIPSE), rotates series to mean equinox of J2000
// ipli		planet number

func rotBack(ipli int) {
	seps2000 := 0.39777715572793088 // sin(eps2000)
	ceps2000 := 0.91748206215761929 // cos(eps2000)
	x := make([][3]float64, MAXORD+1)
	pdp := &swed.Pldat[ipli]
	nco := pdp.Ncoe
	// Calculate time parameters
	t := pdp.Tseg0 + pdp.Dseg/2
	chcfx := pdp.Segp
	chcfy := chcfx[nco:]
	chcfz := chcfx[2*nco:]
	tdiff := (t - pdp.Telem) / 365250.0
	// Calculate rotation parameters
	var qav, pav float64
	if ipli == SEI_MOON {
		dn := pdp.Prot + tdiff*pdp.Dprot
		i := int(dn / TWOPI)
		dn -= float64(i) * TWOPI
		qav = (pdp.Qrot + tdiff*pdp.Dqrot) * math.Cos(dn)
		pav = (pdp.Qrot + tdiff*pdp.Dqrot) * math.Sin(dn)
	} else {
		qav = pdp.Qrot + tdiff*pdp.Dqrot
		pav = pdp.Prot + tdiff*pdp.Dprot
	}
	// Copy coefficients to working array
	for i := 0; i < nco; i++ {
		x[i][0] = chcfx[i]
		x[i][1] = chcfy[i]
		x[i][2] = chcfz[i]
	}
	// Handle elliptical case
	if pdp.Iflg&SEI_FLG_ELLIPSE != 0 {
		refepx := pdp.Refep
		refepy := refepx[nco:]

		omtild := pdp.Peri + tdiff*pdp.Dperi
		i := int(omtild / TWOPI)
		omtild -= float64(i) * TWOPI
		com := math.Cos(omtild)
		som := math.Sin(omtild)

		// Add reference orbit
		for i := 0; i < nco; i++ {
			x[i][0] = chcfx[i] + com*refepx[i] - som*refepy[i]
			x[i][1] = chcfy[i] + com*refepy[i] + som*refepx[i]
		}
	}
	// Construct right-handed orthonormal system
	cosih2 := 1.0 / (1.0 + qav*qav + pav*pav)
	// Calculate orbit pole
	uiz := [3]float64{
		2.0 * pav * cosih2,
		-2.0 * qav * cosih2,
		(1.0 - qav*qav - pav*pav) * cosih2,
	}
	// Calculate origin of longitudes vector
	uix := [3]float64{
		(1.0 + qav*qav - pav*pav) * cosih2,
		2.0 * qav * pav * cosih2,
		-2.0 * pav * cosih2,
	}
	// Calculate vector in orbital plane orthogonal to origin of longitudes
	uiy := [3]float64{
		2.0 * qav * pav * cosih2,
		(1.0 - qav*qav + pav*pav) * cosih2,
		2.0 * qav * cosih2,
	}
	// Rotate to actual orientation in space
	for i := 0; i < nco; i++ {
		xrot := x[i][0]*uix[0] + x[i][1]*uiy[0] + x[i][2]*uiz[0]
		yrot := x[i][0]*uix[1] + x[i][1]*uiy[1] + x[i][2]*uiz[1]
		zrot := x[i][0]*uix[2] + x[i][1]*uiy[2] + x[i][2]*uiz[2]
		if math.Abs(xrot)+math.Abs(yrot)+math.Abs(zrot) >= 1e-14 {
			pdp.Neval = i
		}
		x[i][0] = xrot
		x[i][1] = yrot
		x[i][2] = zrot
		if ipli == SEI_MOON {
			// Rotate to j2000 equator
			x[i][1] = ceps2000*yrot - seps2000*zrot
			x[i][2] = seps2000*yrot + ceps2000*zrot
		}
	}
	// Copy results back
	for i := 0; i < nco; i++ {
		chcfx[i] = x[i][0]
		chcfy[i] = x[i][1]
		chcfz[i] = x[i][2]
	}
}

// ===== 5055 ===== embofs sweph.c-5055 =============================================================================
// Adjust position from Earth-Moon barycenter to Earth
// xemb = hel./bar. position or velocity vectors of emb (input)
//                                                  earth (output)
// xmoon= geocentric position or velocity vector of moon

func embofs(xemb, xmoon *[3]float64) {
	for i := 0; i <= 2; i++ {
		xemb[i] -= xmoon[i] / (EARTH_MOON_MRAT + 1.0)
	}
}

// ===== 5068 ===== nut_matrix sweph.c-5068 ==========================================================================

// nutMatrix calculates the nutation matrix
// * nu		pointer to nutation data structure
// * oe		pointer to epsilon data structure
func nutMatrix(nu *Nut, oe *Epsilon) {
	psi := nu.Nutlo[0]
	eps := oe.Eps + nu.Nutlo[1]

	sinpsi := math.Sin(psi)
	cospsi := math.Cos(psi)
	sineps0 := oe.Seps
	coseps0 := oe.Ceps
	sineps := math.Sin(eps)
	coseps := math.Cos(eps)

	nu.Matrix[0][0] = cospsi
	nu.Matrix[0][1] = sinpsi * coseps
	nu.Matrix[0][2] = sinpsi * sineps

	nu.Matrix[1][0] = -sinpsi * coseps0
	nu.Matrix[1][1] = cospsi*coseps*coseps0 + sineps*sineps0
	nu.Matrix[1][2] = cospsi*sineps*coseps0 - coseps*sineps0

	nu.Matrix[2][0] = -sinpsi * sineps0
	nu.Matrix[2][1] = cospsi*coseps*sineps0 - sineps*coseps0
	nu.Matrix[2][2] = cospsi*sineps*sineps0 + coseps*coseps0
}

// ===== 5857 ===== constants for meff ==  sweph.c-5857 ==============================================================

type MeffEle struct {
	R    float64
	MEff float64
}

// EffArr contains r and m_eff values for photon passing the sun at minimum distance r (fraction of Rsun).
// Values were computed with sun_model.c using classic treatment of a photon passing a gravity field, multiplied by 2.
// The sun mass distribution m(r) is from Michael Stix, The Sun, p. 47.
var EffArr = []MeffEle{
	{1.000, 1.000000},
	{0.990, 0.999979},
	{0.980, 0.999940},
	{0.970, 0.999881},
	{0.960, 0.999811},
	{0.950, 0.999724},
	{0.940, 0.999622},
	{0.930, 0.999497},
	{0.920, 0.999354},
	{0.910, 0.999192},
	{0.900, 0.999000},
	{0.890, 0.998786},
	{0.880, 0.998535},
	{0.870, 0.998242},
	{0.860, 0.997919},
	{0.850, 0.997571},
	{0.840, 0.997198},
	{0.830, 0.996792},
	{0.820, 0.996316},
	{0.810, 0.995791},
	{0.800, 0.995226},
	{0.790, 0.994625},
	{0.780, 0.993991},
	{0.770, 0.993326},
	{0.760, 0.992598},
	{0.750, 0.991770},
	{0.740, 0.990873},
	{0.730, 0.989919},
	{0.720, 0.988912},
	{0.710, 0.987856},
	{0.700, 0.986755},
	{0.690, 0.985610},
	{0.680, 0.984398},
	{0.670, 0.982986},
	{0.660, 0.981437},
	{0.650, 0.979779},
	{0.640, 0.978024},
	{0.630, 0.976182},
	{0.620, 0.974256},
	{0.610, 0.972253},
	{0.600, 0.970174},
	{0.590, 0.968024},
	{0.580, 0.965594},
	{0.570, 0.962797},
	{0.560, 0.959758},
	{0.550, 0.956515},
	{0.540, 0.953088},
	{0.530, 0.949495},
	{0.520, 0.945741},
	{0.510, 0.941838},
	{0.500, 0.937790},
	{0.490, 0.933563},
	{0.480, 0.928668},
	{0.470, 0.923288},
	{0.460, 0.917527},
	{0.450, 0.911432},
	{0.440, 0.905035},
	{0.430, 0.898353},
	{0.420, 0.891022},
	{0.410, 0.882940},
	{0.400, 0.874312},
	{0.390, 0.865206},
	{0.380, 0.855423},
	{0.370, 0.844619},
	{0.360, 0.833074},
	{0.350, 0.820876},
	{0.340, 0.808031},
	{0.330, 0.793962},
	{0.320, 0.778931},
	{0.310, 0.763021},
	{0.300, 0.745815},
	{0.290, 0.727557},
	{0.280, 0.708234},
	{0.270, 0.687583},
	{0.260, 0.665741},
	{0.250, 0.642597},
	{0.240, 0.618252},
	{0.230, 0.592586},
	{0.220, 0.565747},
	{0.210, 0.537697},
	{0.200, 0.508554},
	{0.190, 0.478420},
	{0.180, 0.447322},
	{0.170, 0.415454},
	{0.160, 0.382892},
	{0.150, 0.349955},
	{0.140, 0.316691},
	{0.130, 0.283565},
	{0.120, 0.250431},
	{0.110, 0.218327},
	{0.100, 0.186794},
	{0.090, 0.156287},
	{0.080, 0.128421},
	{0.070, 0.102237},
	{0.060, 0.077393},
	{0.050, 0.054833},
	{0.040, 0.036361},
	{0.030, 0.020953},
	{0.020, 0.009645},
	{0.010, 0.002767},
	{0.000, 0.000000},
}

// ===== 5966 ===== meff sweph.c-5966 ================================================================================

// Meff calculates the effective mass for a given radius r.
// r is the distance as a fraction of solar radius (Rsun).
// Returns a value between 0 and 1 representing the effective mass.
func Meff(r float64) float64 {
	if r <= 0 {
		return 0.0
	} else if r >= 1 {
		return 1.0
	}
	// Find the appropriate interval in EffArr
	var i int
	for i = 0; i < len(EffArr) && EffArr[i].R > r; i++ {
	} // empty body
	// Linear interpolation between points
	f := (r - EffArr[i-1].R) / (EffArr[i].R - EffArr[i-1].R)
	m := EffArr[i-1].MEff + f*(EffArr[i].MEff-EffArr[i-1].MEff)
	return m
}

// ===== 6012 ===== swi_check_ecliptic sweph.c-6012 ==================================================================

func swiCheckEcliptic(tjd float64, iflag int32) {
	if swed.Oec2000.Teps != J2000 {
		calcEpsilon(J2000, iflag, &swed.Oec2000)
	}
	if tjd == J2000 {
		swed.Oec.Teps = swed.Oec2000.Teps
		swed.Oec.Eps = swed.Oec2000.Eps
		swed.Oec.Seps = swed.Oec2000.Seps
		swed.Oec.Ceps = swed.Oec2000.Ceps
		return
	}
	if swed.Oec.Teps != tjd || tjd == 0 {
		calcEpsilon(tjd, iflag, &swed.Oec)
	}
}

// ===== 6029 ===== swi_check_nutation sweph.c-6029 ==================================================================

// swiCheckNutation computes nutation if it is wanted and has not yet been computed.
// If speed flag has been turned on since last computation, nutation is recomputed.
func swiCheckNutation(tjd float64, iflag int32) {
	var (
		speedf1 int32
		speedf2 int32
		nutflag int32 // This would need proper initialization in the actual program
	)
	speedf1 = nutflag & SEFLG_SPEED
	speedf2 = iflag & SEFLG_SPEED
	if (iflag&SEFLG_NONUT) == 0 &&
		(tjd != swed.Nut.Tnut || tjd == 0 ||
			(speedf1 == 0 && speedf2 != 0)) {
		swiNutation(tjd, iflag, swed.Nut.Nutlo)
		swed.Nut.Tnut = tjd
		swed.Nut.Snut = math.Sin(swed.Nut.Nutlo[1])
		swed.Nut.Cnut = math.Cos(swed.Nut.Nutlo[1])
		nutflag = iflag

		nutMatrix(&swed.Nut, &swed.Oec)

		if (iflag & SEFLG_SPEED) != 0 {
			// Once more for 'speed' of nutation, which is needed for
			// planetary speeds
			t := tjd - NUT_SPEED_INTV
			swiNutation(t, iflag, swed.Nutv.Nutlo)
			swed.Nutv.Tnut = t
			swed.Nutv.Snut = math.Sin(swed.Nutv.Nutlo[1])
			swed.Nutv.Cnut = math.Cos(swed.Nutv.Nutlo[1])
			nutMatrix(&swed.Nutv, &swed.Oec)
		}
	}
}

// ===== 6061 ===== plaus_iflag - sweph.c-6061 =======================================================================

// plausIflag corrects nonsensical iflags and completes incomplete iflags
func plausIflag(iflag int32, ipl int32, tjd float64, serr *string) int32 {
	var epheflag int32 = 0

	// Port: removed parts that handle JPL and Moshier

	// if topocentric bit, turn helio- and barycentric bits off
	if (iflag & SEFLG_TOPOCTR) != 0 {
		iflag = iflag &^ (SEFLG_HELCTR | SEFLG_BARYCTR)
	}

	// if barycentric bit, turn heliocentric bit off
	if (iflag & SEFLG_BARYCTR) != 0 {
		iflag = iflag &^ SEFLG_HELCTR
	}
	if (iflag & SEFLG_HELCTR) != 0 {
		iflag = iflag &^ SEFLG_BARYCTR
	}

	// if heliocentric bit, turn aberration and deflection off
	if (iflag & (SEFLG_HELCTR | SEFLG_BARYCTR)) != 0 {
		iflag |= SEFLG_NOABERR | SEFLG_NOGDEFL // iflag |= SEFLG_TRUEPOS;
	}

	// if no_precession bit is set, set also no_nutation bit
	if (iflag & SEFLG_J2000) != 0 {
		iflag |= SEFLG_NONUT
	}

	// if truepos is set, turn off grav. defl. and aberration
	if (iflag & SEFLG_TRUEPOS) != 0 {
		iflag |= (SEFLG_NOGDEFL | SEFLG_NOABERR)
	}

	if epheflag == 0 {
		epheflag = SEFLG_DEFAULTEPH
	}

	iflag = (iflag &^ SEFLG_EPHMASK) | epheflag
	return iflag
}
