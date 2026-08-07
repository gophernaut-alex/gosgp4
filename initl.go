package sgp4

import "math"

// InitlResult holds the initl outputs. At initl's one real call site
// (SGP4.cpp:1479-1482, inside sgp4init/NewElsetRec), Method, Con41, Gsto,
// and NoUnkozai are passed as real satrec fields (and so get copied onto
// ElsetRec by NewElsetRec), while Ainv, Ao, Con42, Cosio, Cosio2, Eccsq,
// Omeosq, Posq, Rp, Rteosq, Sinio are sgp4init's own transient locals, fed
// into its subsequent near-earth coefficient computation and then
// discarded -- the same mixed field/local split as dscom's outputs, so
// initl stays a plain function returning a result struct (like
// DscomResult) rather than an (*ElsetRec) method.
type InitlResult struct {
	Method byte

	Ainv   float64
	Ao     float64
	Con41  float64
	Con42  float64
	Cosio  float64
	Cosio2 float64
	Eccsq  float64
	Omeosq float64
	Posq   float64
	Rp     float64
	Rteosq float64
	Sinio  float64
	Gsto   float64

	NoUnkozai float64
}

// initl performs one-time epoch initialization: un-Kozai's the mean motion,
// computes derived geometric quantities, sets Method='n' (NewElsetRec may
// later override this to 'd' for deep-space orbits), and computes Greenwich
// sidereal time at epoch. Replaces initl in SGP4.cpp (lines 1208-1282).
//
// opsmode is accepted for signature fidelity but unused in this body: the
// reference's opsmode-conditional "old way" sidereal-time branch is
// commented out upstream, leaving an unconditional call to gstime_SGP4 --
// so this port skips that dead branch (it computes a value, "gsto1" in the
// original, that is never actually used for anything).
func initl(xke, j2, ecco, epoch, inclo, noKozai float64, opsmode byte) InitlResult {
	const x2o3 = 2.0 / 3.0

	// ------------- calculate auxillary epoch quantities ----------
	eccsq := ecco * ecco
	omeosq := 1.0 - eccsq
	rteosq := math.Sqrt(omeosq)
	cosio := math.Cos(inclo)
	cosio2 := cosio * cosio

	// ------------------ un-kozai the mean motion -----------------
	ak := math.Pow(xke/noKozai, x2o3)
	d1 := 0.75 * j2 * (3.0*cosio2 - 1.0) / (rteosq * omeosq)
	del := d1 / (ak * ak)
	adel := ak * (1.0 - del*del - del*(1.0/3.0+134.0*del*del/81.0))
	del = d1 / (adel * adel)
	noUnkozai := noKozai / (1.0 + del)

	ao := math.Pow(xke/noUnkozai, x2o3)
	sinio := math.Sin(inclo)
	po := ao * omeosq
	con42 := 1.0 - 5.0*cosio2
	con41 := -con42 - cosio2 - cosio2
	ainv := 1.0 / ao
	posq := po * po
	rp := ao * (1.0 - ecco)

	gsto := GSTime(epoch + 2433281.5)

	return InitlResult{
		Method:    'n',
		Ainv:      ainv,
		Ao:        ao,
		Con41:     con41,
		Con42:     con42,
		Cosio:     cosio,
		Cosio2:    cosio2,
		Eccsq:     eccsq,
		Omeosq:    omeosq,
		Posq:      posq,
		Rp:        rp,
		Rteosq:    rteosq,
		Sinio:     sinio,
		Gsto:      gsto,
		NoUnkozai: noUnkozai,
	}
}
