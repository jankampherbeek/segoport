package segoport

// Port gives access to all the public functions of segoport module.
type Port struct{}

// Version returns the current version of segoport.
func (p *Port) Version() string {
	return "Temporary version 0.0.0"
}
