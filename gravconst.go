package sgp4

import "math"

// pi matches the C++ SGP4.cpp's own #define pi (not math.Pi) so that every
// derived constant (twopi, etc.) reproduces the reference bit-for-bit.
const pi = 3.14159265358979323846

// GravConstType mirrors the C++ gravconsttype enum (SGP4.h).
type GravConstType int

const (
	WGS72Old GravConstType = iota
	WGS72
	WGS84
)

// GravConst holds the earth-gravity constants used throughout propagation.
// Field names match ElsetRec's existing spelling of the same quantities.
type GravConst struct {
	Tumin         float64
	Mus           float64
	Radiusearthkm float64
	Xke           float64
	J2            float64
	J3            float64
	J4            float64
	J3oj2         float64
}

// GetGravConst returns the constant set for the given gravity model,
// replacing getgravconst in SGP4.cpp (lines 2082-2136). Unlike the C++,
// which silently leaves all 8 outputs untouched for an invalid whichconst,
// this returns ErrUnknownGravConst.
func GetGravConst(which GravConstType) (GravConst, error) {
	switch which {
	case WGS72Old:
		return GravConst{
			Tumin:         1.0 / 0.0743669161,
			Mus:           398600.79964,
			Radiusearthkm: 6378.135,
			Xke:           0.0743669161,
			J2:            0.001082616,
			J3:            -0.00000253881,
			J4:            -0.00000165597,
			J3oj2:         -0.00000253881 / 0.001082616,
		}, nil
	case WGS72:
		const mu = 398600.8
		const radiusearthkm = 6378.135
		xke := 60.0 / math.Sqrt(radiusearthkm*radiusearthkm*radiusearthkm/mu)
		return GravConst{
			Tumin:         1.0 / xke,
			Mus:           mu,
			Radiusearthkm: radiusearthkm,
			Xke:           xke,
			J2:            0.001082616,
			J3:            -0.00000253881,
			J4:            -0.00000165597,
			J3oj2:         -0.00000253881 / 0.001082616,
		}, nil
	case WGS84:
		const mu = 398600.5
		const radiusearthkm = 6378.137
		xke := 60.0 / math.Sqrt(radiusearthkm*radiusearthkm*radiusearthkm/mu)
		return GravConst{
			Tumin:         1.0 / xke,
			Mus:           mu,
			Radiusearthkm: radiusearthkm,
			Xke:           xke,
			J2:            0.00108262998905,
			J3:            -0.00000253215306,
			J4:            -0.00000161098761,
			J3oj2:         -0.00000253215306 / 0.00108262998905,
		}, nil
	default:
		return GravConst{}, ErrUnknownGravConst
	}
}
