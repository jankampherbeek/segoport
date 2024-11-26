package segoport

import "testing"

func TestPortVersion(t *testing.T) {
	p := Port{}
	result := p.Version()
	l := len(result)
	if l < 3 {
		t.Errorf("Port.Version should returned a string of %d characters; want 3 characters or more", l)
	}
}
