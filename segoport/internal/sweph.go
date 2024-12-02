package internal

// SwiInitSwedIfStart initializes swed structure.
// Returns true if initialization is done, false otherwise
func SwiInitSwedIfStart() bool {

	if !sweData.SwedIsInitialised {

		sweData.EphePath = SE_EPHE_PATH
		sweData.JplFnam = SE_FNAME_DFT

		SetTidAcc(SE_TIDAL_AUTOMATIC)
		sweData.SwedIsInitialised = true
		return true
	}
	return false
}
