package internal

// Partial set of constatns from swedef.h

const (
	M_PI     = 3.14159265358979323846
	RADTODEG = 180.0 / M_PI
	DEGTORAD = M_PI / 180.0
	AS_MAXCH = 256 // used for string declarations, allowing 255 char+\0
)
