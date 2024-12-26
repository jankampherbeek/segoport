package internal

import (
	"strconv"
	"strings"
	"unicode"
)

// ===== 0522 ===== plan_oscu_elem swemplan.c-0522 ===================================================================
// Port: important, in the original code there is an ifdef pragma, with parts for SE_NEELY and another part.
// I splitted this in two ragnges of constants: planOscuElem and planOscuElemNeely

var planOscuElem = [SE_NFICT_ELEM][8]float64{
	{J1900, J1900, 104.5959, 40.99837, 0, 0, 0, 0},  /* Cupido   */
	{J1900, J1900, 337.4517, 50.667443, 0, 0, 0, 0}, /* Hades    */
	{J1900, J1900, 104.0904, 59.214362, 0, 0, 0, 0}, /* Zeus     */
	{J1900, J1900, 17.7346, 64.816896, 0, 0, 0, 0},  /* Kronos   */
	{J1900, J1900, 138.0354, 70.361652, 0, 0, 0, 0}, /* Apollon  */
	{J1900, J1900, -8.678, 73.736476, 0, 0, 0, 0},   /* Admetos  */
	{J1900, J1900, 55.9826, 77.445895, 0, 0, 0, 0},  /* Vulkanus */
	{J1900, J1900, 165.3595, 83.493733, 0, 0, 0, 0}, /* Poseidon */
	// Isis-Transpluto; elements from "Die Sterne" 3/1952, p. 70ff. Strubell does not give an equinox. 1945 is taken to best
	// reproduce ASTRON ephemeris. (This is a strange choice, though.)
	// The epoch is 1772.76. The year is understood to have 366 days. The fraction is counted from 1 Jan. 1772 */
	{2368547.66, 2431456.5, 0.0, 77.775, 0.3, 0.7, 0, 0},
	// Nibiru, elements from Christian Woeltge, Hannover
	{1856113.380954, 1856113.380954, 0.0, 234.8921, 0.981092, 103.966, -44.567, 158.708},
	// Harrington, elements from Astronomical Journal 96(4), Oct. 1988
	{2374696.5, J2000, 0.0, 101.2, 0.411, 208.5, 275.4, 32.4},
	// Leverrier's Neptune, according to W.G. Hoyt, "Planets X and Pluto", Tucson 1980, p. 63
	{2395662.5, 2395662.5, 34.05, 36.15, 0.10761, 284.75, 0, 0},
	// Adam's Neptune
	{2395662.5, 2395662.5, 24.28, 37.25, 0.12062, 299.11, 0, 0},
	// Lowell's Pluto
	{2425977.5, 2425977.5, 281, 43.0, 0.202, 204.9, 0, 0},
	// Pickering's Pluto
	{2425977.5, 2425977.5, 48.95, 55.1, 0.31, 280.1, 100, 15},
}

// ===== 0522 ===== plan_oscu_elem swemplan.c-0522 ===================================================================
// Port: important, in the original code there is an ifdef pragma, with parts for SE_NEELY and another part.
// I splitted this in two ragnges of constants: planOscuElem and planOscuElemNeely

var planOscuElemNeely = [SE_NFICT_ELEM][8]float64{
	{J1900, J1900, 163.7409, 40.99837, 0.00460, 171.4333, 129.8325, 1.0833}, /* Cupido Neely */
	{J1900, J1900, 27.6496, 50.66744, 0.00245, 148.1796, 161.3339, 1.0500},  /* Hades Neely */
	{J1900, J1900, 165.1232, 59.21436, 0.00120, 299.0440, 0.0000, 0.0000},   /* Zeus Neely */
	{J1900, J1900, 169.0193, 64.81960, 0.00305, 208.8801, 0.0000, 0.0000},   /* Kronos Neely */
	{J1900, J1900, 138.0533, 70.29949, 0.00000, 0.0000, 0.0000, 0.0000},     /* Apollon Neely */
	{J1900, J1900, 351.3350, 73.62765, 0.00000, 0.0000, 0.0000, 0.0000},     /* Admetos Neely */
	{J1900, J1900, 55.8983, 77.25568, 0.00000, 0.0000, 0.0000, 0.0000},      /* Vulcanus Neely */
	{J1900, J1900, 165.5163, 83.66907, 0.00000, 0.0000, 0.0000, 0.0000},     /* Poseidon Neely */
	// Isis-Transpluto; elements from "Die Sterne" 3/1952, p. 70ff. Strubell does not give an equinox. 1945 is taken to best
	// reproduce ASTRON ephemeris. (This is a strange choice, though.)
	// The epoch is 1772.76. The year is understood to have 366 days. The fraction is counted from 1 Jan. 1772 */
	{2368547.66, 2431456.5, 0.0, 77.775, 0.3, 0.7, 0, 0},
	// Nibiru, elements from Christian Woeltge, Hannover
	{1856113.380954, 1856113.380954, 0.0, 234.8921, 0.981092, 103.966, -44.567, 158.708},
	// Harrington, elements from Astronomical Journal 96(4), Oct. 1988
	{2374696.5, J2000, 0.0, 101.2, 0.411, 208.5, 275.4, 32.4},
	// Leverrier's Neptune, according to W.G. Hoyt, "Planets X and Pluto", Tucson 1980, p. 63
	{2395662.5, 2395662.5, 34.05, 36.15, 0.10761, 284.75, 0, 0},
	// Adam's Neptune
	{2395662.5, 2395662.5, 24.28, 37.25, 0.12062, 299.11, 0, 0},
	// Lowell's Pluto
	{2425977.5, 2425977.5, 281, 43.0, 0.202, 204.9, 0, 0},
	// Pickering's Pluto
	{2425977.5, 2425977.5, 48.95, 55.1, 0.31, 280.1, 100, 15},
}

// ===== 0916 ===== check_t_terms swemplan.c-0916 ====================================================================

func checkTTerms(t float64, sinp string) (float64, bool) {
	isgn := 1
	var doutp float64 = 0
	fac := 1.0
	z := 0
	tt := make([]float64, 5)

	tt[0] = t / 36525
	tt[1] = tt[0]
	tt[2] = tt[1] * tt[1]
	tt[3] = tt[2] * tt[1]
	tt[4] = tt[3] * tt[1]

	hasAdditionalTerms := strings.ContainsAny(sinp, "+-")

	sp := []rune(strings.TrimSpace(sinp))
	pos := 0

	for pos < len(sp) {
		for pos < len(sp) && (sp[pos] == ' ' || sp[pos] == '\t') {
			pos++
		}
		if pos >= len(sp) {
			break
		}
		if sp[pos] == '+' || sp[pos] == '-' || pos == len(sp) {
			if z > 0 {
				doutp += fac
			}
			isgn = 1
			if pos < len(sp) && sp[pos] == '-' {
				isgn = -1
			}
			fac = float64(isgn)
			if pos >= len(sp) {
				return doutp, hasAdditionalTerms
			}
			pos++
		} else {
			for pos < len(sp) && (sp[pos] == '*' || sp[pos] == ' ' || sp[pos] == '\t') {
				pos++
			}
			if pos >= len(sp) {
				break
			}
			// Handle 'T' terms
			if sp[pos] == 't' || sp[pos] == 'T' {
				pos++
				if pos < len(sp) {
					if sp[pos] == '+' || sp[pos] == '-' {
						fac *= tt[0]
					} else {
						numStr := ""
						for pos < len(sp) && unicode.IsDigit(sp[pos]) {
							numStr += string(sp[pos])
							pos++
						}
						if numStr != "" {
							i, _ := strconv.Atoi(numStr)
							if i >= 0 && i <= 4 {
								fac *= tt[i]
							}
						}
						continue
					}
				}
			} else {
				// Handle numbers
				numStr := ""
				for pos < len(sp) && (unicode.IsDigit(sp[pos]) || sp[pos] == '.') {
					numStr += string(sp[pos])
					pos++
				}
				if numStr != "" {
					num, err := strconv.ParseFloat(numStr, 64)
					if err == nil && (num != 0 || numStr == "0") {
						fac *= num
					}
				}
				continue
			}
		}
		z++
	}
	if z > 0 {
		doutp += fac
	}
	return doutp, hasAdditionalTerms
}
