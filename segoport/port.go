package segoport

import "github.com/jankampherbeek/segoport/internal"

// Port gives access to all the public functions of segoport module.
type Port struct{}

// Version returns the current version of segoport.
func (p *Port) Version() string {
	return "Temporary version 0.0.0"
}

// SweJulDay returns the Julian Day Number.
// Input: year, month, day as int, hour as number with fraction and gregflag to indicate the calender: 1 = Gregorian, 0 = Julian.
// Ouput: the calculated Julian Day Number.
func (p *Port) SweJulDay(year, month, day int, hour float64, gregflag int) float64 {
	d := internal.SweDate{}
	return d.SweJulday(year, month, day, hour, gregflag)
}
