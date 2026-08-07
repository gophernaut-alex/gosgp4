package sgp4

import "errors"

// Sentinel errors returned by Propagate, matching the satrec.error codes
// assigned in sgp4() (SGP4.cpp lines 1844-2042). Checked in this exact order
// on every Propagate call: ErrMeanMotion, ErrEccentricity, (deep-space only)
// ErrPerturbedEccentricity, ErrSemiLatusRectum, then (checked last, after r/v
// are computed) ErrDecayed.
//
// Code 5 is dead code in the reference (assigned nowhere, commented out at
// SGP4.cpp:1494) and has no corresponding sentinel here.
var (
	ErrEccentricity          = errors.New("sgp4: mean eccentricity out of range")      // code 1
	ErrMeanMotion            = errors.New("sgp4: mean motion less than zero")          // code 2
	ErrPerturbedEccentricity = errors.New("sgp4: perturbed eccentricity out of range") // code 3
	ErrSemiLatusRectum       = errors.New("sgp4: semi-latus rectum less than zero")    // code 4
	ErrDecayed               = errors.New("sgp4: satellite has decayed")               // code 6

	ErrUnknownGravConst = errors.New("sgp4: unknown gravity constant model")
)
