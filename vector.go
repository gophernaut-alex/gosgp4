package sgp4

import "math"

// Undefined is the sentinel value the reference uses to flag a
// geometrically undefined result (near-zero vector in Angle, or any of the
// orbit-type-dependent elements in RVToCOE that don't apply to the given
// orbit's geometry), matching the C++ "undefined = 999999.1" constant used
// in both angle_SGP4 and rv2coe_SGP4.
const Undefined = 999999.1

// undefinedAnomaly is the sentinel NewtonNu leaves e0/m at when neither
// input is touched (matches the C++ "infinite"/initial 999999.9 value used
// in newtonnu_SGP4 and rv2coe_SGP4's semi-major-axis special case).
const undefinedAnomaly = 999999.9

// Vector3 is a 3-element Cartesian vector (km, km/s, or dimensionless
// depending on context), replacing the C++'s double[3] convention.
type Vector3 struct {
	X, Y, Z float64
}

// Dot returns the dot product of v and o, replacing dot_SGP4.
func (v Vector3) Dot(o Vector3) float64 {
	return v.X*o.X + v.Y*o.Y + v.Z*o.Z
}

// Cross returns the cross product v x o, replacing cross_SGP4.
func (v Vector3) Cross(o Vector3) Vector3 {
	return Vector3{
		X: v.Y*o.Z - v.Z*o.Y,
		Y: v.Z*o.X - v.X*o.Z,
		Z: v.X*o.Y - v.Y*o.X,
	}
}

// Magnitude returns the Euclidean norm of v, replacing mag_SGP4.
func (v Vector3) Magnitude() float64 {
	return math.Sqrt(v.X*v.X + v.Y*v.Y + v.Z*v.Z)
}

// Add returns v + o.
func (v Vector3) Add(o Vector3) Vector3 {
	return Vector3{v.X + o.X, v.Y + o.Y, v.Z + o.Z}
}

// Sub returns v - o.
func (v Vector3) Sub(o Vector3) Vector3 {
	return Vector3{v.X - o.X, v.Y - o.Y, v.Z - o.Z}
}

// Scale returns v scaled by k.
func (v Vector3) Scale(k float64) Vector3 {
	return Vector3{v.X * k, v.Y * k, v.Z * k}
}

// Sgn returns -1 if x < 0, else 1, replacing sgn_SGP4.
func Sgn(x float64) float64 {
	if x < 0.0 {
		return -1.0
	}
	return 1.0
}

// Angle returns the angle in radians (0..pi) between v1 and v2, or
// Undefined if either vector's magnitude is (near) zero, replacing
// angle_SGP4.
func Angle(v1, v2 Vector3) float64 {
	const small = 0.00000001

	magv1 := v1.Magnitude()
	magv2 := v2.Magnitude()

	if magv1*magv2 > small*small {
		temp := v1.Dot(v2) / (magv1 * magv2)
		if math.Abs(temp) > 1.0 {
			temp = Sgn(temp) * 1.0
		}
		return math.Acos(temp)
	}
	return Undefined
}

// NewtonNu solves Kepler's equation when the true anomaly nu (radians) and
// eccentricity ecc are known, returning the eccentric/parabolic/hyperbolic
// anomaly e0 and mean anomaly m, replacing newtonnu_SGP4. In the one branch
// the reference itself leaves undefined (hyperbolic orbit with |nu| too
// close to the asymptote), e0 and m are returned as undefinedAnomaly
// (999999.9), matching the reference exactly rather than introducing a Go
// error for a case upstream doesn't treat as one.
func NewtonNu(ecc, nu float64) (e0, m float64) {
	const small = 0.00000001

	e0 = undefinedAnomaly
	m = undefinedAnomaly

	switch {
	case math.Abs(ecc) < small:
		// --------------------------- circular ------------------------
		m = nu
		e0 = nu
	case ecc < 1.0-small:
		// ---------------------- elliptical -----------------------
		sine := (math.Sqrt(1.0-ecc*ecc) * math.Sin(nu)) / (1.0 + ecc*math.Cos(nu))
		cose := (ecc + math.Cos(nu)) / (1.0 + ecc*math.Cos(nu))
		e0 = math.Atan2(sine, cose)
		m = e0 - ecc*math.Sin(e0)
	case ecc > 1.0+small:
		// -------------------- hyperbolic --------------------
		if fabs := math.Abs(nu); fabs+0.00001 < pi-math.Acos(1.0/ecc) {
			sine := (math.Sqrt(ecc*ecc-1.0) * math.Sin(nu)) / (1.0 + ecc*math.Cos(nu))
			e0 = math.Asinh(sine)
			m = ecc*math.Sinh(e0) - e0
		}
	default:
		// ----------------- parabolic ---------------------
		if math.Abs(nu) < 168.0*pi/180.0 {
			e0 = math.Tan(nu * 0.5)
			m = e0 + (e0*e0*e0)/3.0
		}
	}

	if ecc < 1.0 {
		m = math.Mod(m, 2.0*pi)
		if m < 0.0 {
			m += 2.0 * pi
		}
		e0 = math.Mod(e0, 2.0*pi)
	}
	return e0, m
}

// COE holds classical orbital elements as returned by RVToCOE. Fields that
// don't apply to a given orbit's geometry (e.g. argp for an equatorial
// orbit) are set to Undefined, matching the C++'s sentinel behavior.
type COE struct {
	P, A, Ecc, Incl, Omega, Argp, Nu, M, ArgLat, TrueLon, LonPer float64
}

// RVToCOE computes classical orbital elements from position r (km) and
// velocity v (km/s) vectors and gravitational parameter mus (km^3/s^2),
// replacing rv2coe_SGP4.
func RVToCOE(r, v Vector3, mus float64) COE {
	const small = 0.00000001
	const infinite = undefinedAnomaly
	twopi := 2.0 * pi
	halfpi := 0.5 * pi

	magr := r.Magnitude()

	hbar := r.Cross(v)
	magh := hbar.Magnitude()
	if magh <= small {
		return COE{
			P: Undefined, A: Undefined, Ecc: Undefined, Incl: Undefined,
			Omega: Undefined, Argp: Undefined, Nu: Undefined, M: Undefined,
			ArgLat: Undefined, TrueLon: Undefined, LonPer: Undefined,
		}
	}

	nbar := Vector3{X: -hbar.Y, Y: hbar.X, Z: 0.0}
	magn := nbar.Magnitude()
	magv := v.Magnitude()
	c1 := magv*magv - mus/magr
	rdotv := r.Dot(v)
	ebar := Vector3{
		X: (c1*r.X - rdotv*v.X) / mus,
		Y: (c1*r.Y - rdotv*v.Y) / mus,
		Z: (c1*r.Z - rdotv*v.Z) / mus,
	}
	ecc := ebar.Magnitude()

	// ------------  find a e and semi-latus rectum   ----------
	sme := (magv * magv * 0.5) - (mus / magr)
	var a float64
	if math.Abs(sme) > small {
		a = -mus / (2.0 * sme)
	} else {
		a = infinite
	}
	p := magh * magh / mus

	// -----------------  find inclination   -------------------
	hk := hbar.Z / magh
	incl := math.Acos(hk)

	// --------  determine type of orbit for later use  --------
	// typeorbit: 1 = "ei", 2 = "ce", 3 = "ci", 4 = "ee"
	typeorbit := 1
	if ecc < small {
		if incl < small || math.Abs(incl-pi) < small {
			typeorbit = 2
		} else {
			typeorbit = 3
		}
	} else if incl < small || math.Abs(incl-pi) < small {
		typeorbit = 4
	}

	// ----------  find right ascension of the ascending node ------------
	var omega float64
	if magn > small {
		temp := nbar.X / magn
		if math.Abs(temp) > 1.0 {
			temp = Sgn(temp)
		}
		omega = math.Acos(temp)
		if nbar.Y < 0.0 {
			omega = twopi - omega
		}
	} else {
		omega = Undefined
	}

	// ---------------- find argument of perigee ---------------
	var argp float64
	if typeorbit == 1 {
		argp = Angle(nbar, ebar)
		if ebar.Z < 0.0 {
			argp = twopi - argp
		}
	} else {
		argp = Undefined
	}

	// ------------  find true anomaly at epoch    -------------
	var nu float64
	if typeorbit == 1 || typeorbit == 4 {
		nu = Angle(ebar, r)
		if rdotv < 0.0 {
			nu = twopi - nu
		}
	} else {
		nu = Undefined
	}

	// ----  find argument of latitude - circular inclined -----
	var arglat, m float64
	if typeorbit == 3 {
		arglat = Angle(nbar, r)
		if r.Z < 0.0 {
			arglat = twopi - arglat
		}
		m = arglat
	} else {
		arglat = Undefined
	}

	// -- find longitude of perigee - elliptical equatorial ----
	var lonper float64
	if ecc > small && typeorbit == 4 {
		temp := ebar.X / ecc
		if math.Abs(temp) > 1.0 {
			temp = Sgn(temp)
		}
		lonper = math.Acos(temp)
		if ebar.Y < 0.0 {
			lonper = twopi - lonper
		}
		if incl > halfpi {
			lonper = twopi - lonper
		}
	} else {
		lonper = Undefined
	}

	// -------- find true longitude - circular equatorial ------
	var truelon float64
	if magr > small && typeorbit == 2 {
		temp := r.X / magr
		if math.Abs(temp) > 1.0 {
			temp = Sgn(temp)
		}
		truelon = math.Acos(temp)
		if r.Y < 0.0 {
			truelon = twopi - truelon
		}
		if incl > halfpi {
			truelon = twopi - truelon
		}
		m = truelon
	} else {
		truelon = Undefined
	}

	// ------------ find mean anomaly for all orbits -----------
	if typeorbit == 1 || typeorbit == 4 {
		_, m = NewtonNu(ecc, nu)
	}

	return COE{
		P: p, A: a, Ecc: ecc, Incl: incl, Omega: omega, Argp: argp,
		Nu: nu, M: m, ArgLat: arglat, TrueLon: truelon, LonPer: lonper,
	}
}
