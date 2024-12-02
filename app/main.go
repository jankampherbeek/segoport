/*
 *  SeGoPort.
 *  Copyright (c) Jan Kampherbeek and AstroDienst.
 *  SegoPort is open source.
 *  Please check the file copyright.txt in the root of the source for further details.
 */

package main

import (
	"fmt"
	"github.com/jankampherbeek/segoport"
)

func main() {

	Version()
}

func Version() {
	p := segoport.Port{}
	fmt.Println("Version : " + p.Version())
	jd := p.SweJulDay(2000, 1, 1, 12.0, 1)
	fmt.Printf("Julian day : %f\n", jd)
}
