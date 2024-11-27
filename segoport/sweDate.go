package segoport

import (
	"bufio"
	"errors"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)

const (
	SE_JUL_CAL  = 0
	SE_GREG_CAL = 1
	OK          = 0
	ERR         = -1
)

/*
  Transpiled from swe_date_conversion, original comments:
  This function converts some date+time input {y,m,d,uttime} into the Julian day number tjd.
  The function checks that the input is a legal combination   of dates; for illegal dates like 32 January 1993 it
  returns ERR but still converts the date correctly, i.e. like 1 Feb 1993.
  The function is usually used to convert user input of birth data into the Julian day number.
  Illegal dates should be notified to the user.

  Be aware that we always use astronomical year numbering for the years before Christ, not the historical year numbering.
  Astronomical years are done with negative numbers, historical years with indicators BC or BCE (before common era).
  Year 0 (astronomical)  	= 1 BC historical.
  year -1 (astronomical) 	= 2 BC
  etc.
  Many users of Astro programs do not know about this difference.
  Return: OK or ERR (for illegal date)
*/

// SweDateConversion converts date and time into JD and also checks for the validity of the date.
// Parameters: y, m, d: year, month and date, uttime: UT, c: 'g' for Gregorian calendar, 'j' for Julian.
func (p *Port) SweDateConversion(y, m, d int, uttime float64, c rune) (float64, error) {
	gregflag := SE_JUL_CAL
	if c == 'g' {
		gregflag = SE_GREG_CAL
	}
	jd := p.SweJulday(y, m, d, uttime, gregflag)
	ryear, rmon, rday, _ := p.SweRevjul(jd, gregflag) // ignore value for uttime (rut in SE)
	if rmon == m && rday == d && ryear == y {
		return jd, nil
	} else {
		return jd, errors.New("illegal date")
	}
}

/*
 Transpiled from swe_julday, original comments:
This function returns the absolute Julian day number (JD)
 * for a given calendar date.
 * The arguments are a calendar date: day, month, year as integers, hour as double with decimal fraction.
 * If gregflag = SE_GREG_CAL (1), Gregorian calendar is assumed,
 * if gregflag = SE_JUL_CAL (0),Julian calendar is assumed.

 The Julian day number is a system of numbering all days continously  within the time range of known human history.
 It should be familiar to every astrological or astronomical programmer. The time variable in astronomical theories is
 usually expressed in Julian days or Julian centuries (36525 days per century) relative to some start day;  the start
 day is called 'the epoch'.
 The Julian day number is a double representing the number of days since JD = 0.0 on 1 Jan -4712, 12:00 noon (in the
 Julian calendar).

 Midnight has always a JD with fraction .5, because traditionally the astronomical day started at noon.
 This was practical because then there was no change of date during a night at the telescope. From this comes also the
 fact the noon ephemerides were printed before midnight ephemerides were introduced early in the 20th century.

 NOTE: The Julian day number must not be confused with the Julian calendar system.

 Be aware the we always use astronomical year numbering for the years before Christ, not the historical year numbering.
 Astronomical years are done with negative numbers, historical years with indicators BC or BCE (before common era).
 Year 0 (astronomical)  	= 1 BC
 year -1 (astronomical) 	= 2 BC
 etc.

 Original author: Marc Pottenger, Los Angeles.
 with bug fix for year < -4711   15-aug-88 by Alois Treindl

 References: Oliver Montenbruck, Grundlagen der Ephemeridenrechnung, Verlag Sterne und Weltraum (1987), p.49 ff

 related functions: swe_revjul() reverse Julian day number: compute the calendar date from a given JD
	            swe_date_conversion() includes test for legal date values and notifies errors like 32 January.
 ****************************************************************/

// SweJulday returns the JD for a given date, time and calendar.
// hour represents the hour with a decimal fraction. gregflag is 1 for Gregorian calendar and 0 for Julian calendar
// (constants: SE_GREG_CAL and SE_JUL_CAL).
func (p *Port) SweJulday(year, month, day int, hour float64, gregflag int) float64 {
	u := float64(year)
	if month < 3 {
		u -= 1
	}
	u0 := u + 4712.0
	u1 := float64(month) + 1.0
	if u1 < 4 {
		u1 += 12.0
	}
	jd := math.Floor(u0*365.25) + math.Floor(30.6*u1+0.000001) + float64(day) + hour/24.0 - 63.5
	if gregflag == SE_GREG_CAL {
		u2 := math.Floor(math.Abs(u)/100) - math.Floor(math.Abs(u)/400)
		if u < 0.0 {
			u2 = -u2
		}
		jd = jd - u2 + 2
		if u < 0.0 && u/100 == math.Floor(u/100) && u/400 != math.Floor(u/400) {
			jd -= 1
		}
	}
	return jd
}

/*
  Transpiled from swe_revjul. Orignal comments:
  swe_revjul() is the inverse function to swe_julday(), see the description there.
  Arguments are julian day number, calendar flag (0=julian, 1=gregorian) return values are the calendar day, month,
  year and the hour of the day with decimal fraction (0 .. 23.999999).

  Be aware the we use astronomical year numbering for the years before Christ, not the historical year numbering.
  Astronomical years are done with negative numbers, historical years with indicators BC or BCE (before common era).
  Year  0 (astronomical)  	= 1 BC historical year
  year -1 (astronomical) 	= 2 BC historical year
  year -234 (astronomical) 	= 235 BC historical year
  etc.

  Original author Mark Pottenger, Los Angeles. with bug fix for year < -4711 16-aug-88 Alois Treindl
*/

// SweRevjul returns date and time for a given jd value and calendar.
// gregflag is 1 for Gregorian calendar and 0 for Julian calendar (constants: SE_GREG_CAL and SE_JUL_CAL).
// The return values are day, month, year and ut.
func (p *Port) SweRevjul(jd float64, gregflag int) (int, int, int, float64) {
	u0 := jd + 32082.5
	if gregflag == SE_GREG_CAL {
		u1 := u0 + math.Floor(u0/36525.0) - math.Floor(u0/146100.0) - 38.0
		if jd >= 1830691.5 {
			u1 += 1
		}
		u0 = u0 + math.Floor(u1/36525.0) - math.Floor(u1/146100.0) - 38.0
	}
	u2 := math.Floor(u0 + 123.0)
	u3 := math.Floor((u2 - 122.2) / 365.25)
	u4 := math.Floor((u2 - math.Floor(365.25*u3)) / 30.6001)
	jmon := int(u4 - 1.0)
	if jmon > 12 {
		jmon -= 12
	}
	jday := int(u2 - math.Floor(365.25*u3) - math.Floor(30.6001*u4))
	jyear := int(u3 + math.Floor((u4-2.0)/12.0) - 4800)
	jut := (jd - math.Floor(jd+0.5) + 0.5) * 24.0
	return jyear, jmon, jday, jut
}

/* Transpiled from swe_utc_time_zone, original comments:
   transform local time to UTC or UTC to local time
   input
   iyear ... dsec     date and time
   d_timezone		timezone offset
   output
   iyear_out ... dsec_out

   For time zones east of Greenwich, d_timezone is positive.
   For time zones west of Greenwich, d_timezone is negative.

   For conversion from local time to utc, use +d_timezone.
   For conversion from utc to local time, use -d_timezone.
*/

// SweUtcTimeZone converts between timezone and utc.
// dsec is a decimal swecond, dTimezone is the offseet for the timezone: east is positive, west is negative.
// Use +dTimeZone for conversion from local time to utc, and -dTimeZone for conversion from utc to local time.
// Output: year, month, day, hour, minute, decimal second.
func (p *Port) SweUtcTimeZone(iyear, imonth, iday, ihour, imin int, dsec, dTimezone float64) (int, int, int, int, int, float64) {
	haveLeapsec := false
	if dsec >= 60.0 {
		haveLeapsec = true
		dsec -= 1.0
	}
	dhour := float64(ihour) + float64(imin)/60.0 + dsec/3600.0
	tjd := p.SweJulday(iyear, imonth, iday, 0, SE_GREG_CAL)
	dhour -= dTimezone
	if dhour < 0.0 {
		tjd -= 1.0
		dhour += 24.0
	}
	if dhour >= 24.0 {
		tjd += 1.0
		dhour -= 24.0
	}
	iyearOut, imonthOut, idayOut, _ := p.SweRevjul(tjd+0.001, SE_GREG_CAL)
	ihourOut := int(dhour)
	d := (dhour - float64(ihourOut)) * 60
	iminOut := int(d)
	dsecOut := (d - float64(iminOut)) * 60
	if haveLeapsec {
		dsecOut += 1.0
	}
	return iyearOut, imonthOut, idayOut, ihourOut, iminOut, dsecOut
}

const (
	NLEAP_SECONDS       = 27
	NLEAP_SECONDS_SPACE = 100
)

var leapSeconds = []int{ // keep last zero value as end mark
	19720630, 19721231, 19731231, 19741231, 19751231, 19761231, 19771231,
	19781231, 19791231, 19810630, 19820630, 19830630, 19850630, 19871231,
	19891231, 19901231, 19920630, 19930630, 19940630, 19951231, 19970630,
	19981231, 20051231, 20081231, 20120630, 20150630, 20161231, 0,
}

var initLeapSecondsDone = false

/* Transpiled from int_leapsec, original comment:
   Read additional leap second dates from external file, if given.
*/
// initLeapSec Reads additional leapseconds from seleapsec.txt
func initLeapSec() int {
	if initLeapSecondsDone {
		return len(leapSeconds)
	}
	initLeapSecondsDone = true

	file, err := os.Open("seleapsec.txt")
	if err != nil {
		// File not found or other error, return current size
		return len(leapSeconds)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		ndat, err := strconv.Atoi(line)
		if err != nil {
			continue
		}
		if ndat > leapSeconds[len(leapSeconds)-1] {
			if len(leapSeconds) < NLEAP_SECONDS_SPACE {
				leapSeconds = append(leapSeconds, ndat)
			} else {
				break
			}
		}
	}

	if err := scanner.Err(); err != nil {
		fmt.Println("Error reading file:", err)
	}

	return len(leapSeconds)
}

const (
	J1972      = 2441317.5
	NLEAP_INIT = 10
)

/* Transpiled from swe_utc_to_jd, original comments:
 * Input:  Clock time UTC, year, month, day, hour, minute, second (decimal).
 *         gregflag  Calendar flag
 *         serr      error string
 * Output: An array of doubles:
 *         dret[0] = Julian day number TT (ET)
 *         dret[1] = Julian day number UT1
 *
 * Function returns OK or Error.
 *
 * - Before 1972, swe_utc_to_jd() treats its input time as UT1.
 *   Note: UTC was introduced in 1961. From 1961 - 1971, the length of the UTC second was regularly changed, so that
 *   UTC remained very close to UT1.
 * - From 1972 on, input time is treated as UTC.
 * - If delta_t - nleap - 32.184 > 1, the input time is treated as UT1.
 *   Note: Like this we avoid errors greater than 1 second in case that the leap seconds table (or the Swiss Ephemeris
 *   version) is not updated for a long time.
 */

// SweUtcToJd calculated Jd for ET and Jd for Ut for a given date, time and calendar.
// Input parameter dsec contains seconds with a decimal fraction, gregflag is 1 for Gregorian calendar and 0 for Julian
//
//	calendar (constants: SE_GREG_CAL and SE_JUL_CAL).
func (p *Port) SweUtcToJd(iyear, imonth, iday, ihour, imin int, dsec float64, gregflag int) (float64, float64, error) {
	var tjdUt1, tjdEt, tjdEt1972, dhour, d float64
	var iyear2, imonth2, iday2 int
	var nleap, ndat, tabsizNleap int

	// Error handling: invalid iyear etc.
	tjdUt1 = p.SweJulday(iyear, imonth, iday, 0, gregflag)
	iyear2, imonth2, iday2, _ = p.SweRevjul(tjdUt1, gregflag)
	if iyear != iyear2 || imonth != imonth2 || iday != iday2 {
		return 0, 0, errors.New(fmt.Sprintf("invalid date: year = %d, month = %d, day = %d", iyear, imonth, iday))
	}
	if ihour < 0 || ihour > 23 || imin < 0 || imin > 59 || dsec < 0 || dsec >= 61 || (dsec >= 60 && (imin < 59 || ihour < 23 || tjdUt1 < J1972)) {
		return 0, 0, errors.New(fmt.Sprintf("invalid time: %d:%d:%.2f", ihour, imin, dsec))
	}

	dhour = float64(ihour) + float64(imin)/60.0 + dsec/3600.0

	// Before 1972, treat input date as UT1
	if tjdUt1 < J1972 {
		tjdUt1 = p.SweJulday(iyear, imonth, iday, dhour, gregflag)
		tjdEt = tjdUt1 + sweDeltatEx(tjdUt1, -1, nil)
		return tjdEt, tjdUt1, nil
	}

	// Convert to Gregorian calendar if needed
	if gregflag == SE_JUL_CAL {
		gregflag = SE_GREG_CAL
		iyear, imonth, iday, _ = p.SweRevjul(tjdUt1, gregflag)
	}

	// Number of leap seconds since 1972
	tabsizNleap = initLeapSec()
	nleap = NLEAP_INIT
	ndat = iyear*10000 + imonth*100 + iday
	for i := 0; i < tabsizNleap; i++ {
		if ndat <= leapSeconds[i] {
			break
		}
		nleap++
	}

	// Check if delta_t - nleap - 32.184 > 0.9
	d = sweDeltatEx(tjdUt1, -1, nil) * 86400.0
	if d-float64(nleap)-32.184 >= 1.0 {
		tjdUt1 += dhour / 24.0
		tjdEt = tjdUt1 + sweDeltatEx(tjdUt1, -1, nil)
		return tjdEt, tjdUt1, nil
	}

	// Check for valid leap second
	if dsec >= 60 {
		validLeapSecond := false
		for i := 0; i < tabsizNleap; i++ {
			if ndat == leapSeconds[i] {
				validLeapSecond = true
				break
			}
		}
		if !validLeapSecond {
			return 0, 0, errors.New(fmt.Sprintf("invalid time (no leap second!): %d:%d:%.2f", ihour, imin, dsec))
		}
	}

	// Convert UTC to ET and UT1
	d = tjdUt1 - J1972
	d += float64(ihour)/24.0 + float64(imin)/1440.0 + dsec/86400.0
	tjdEt1972 = J1972 + (32.184+NLEAP_INIT)/86400.0
	tjdEt = tjdEt1972 + d + float64(nleap-NLEAP_INIT)/86400.0
	d = sweDeltatEx(tjdEt, -1, nil)
	tjdUt1 = tjdEt - sweDeltatEx(tjdEt-d, -1, nil)
	tjdUt1 = tjdEt - sweDeltatEx(tjdUt1, -1, nil)

	return tjdEt, tjdUt1, nil
}

func sweDeltatEx(jd float64, iflag int, serr *string) float64 {
	// Placeholder for delta T calculation
	// This function should return the difference between TT and UT1
	return 0.0 // TODO call actual delta T calculation: swe_deltat_ex in swiphlib.c
}
