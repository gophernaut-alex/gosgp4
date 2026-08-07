package sgp4

import "math"

// OpsMode selects between the AFSPC and "improved" sidereal-time/day-number
// conventions used during initialization, replacing the C++'s bare 'a'/'i'
// opsmode char.
type OpsMode byte

const (
	OpsModeAFSPC    OpsMode = 'a'
	OpsModeImproved OpsMode = 'i'
)

// NewElsetRec builds and fully initializes a satellite record from raw mean
// elements, replacing sgp4init in SGP4.cpp (lines 1366-1664). Units match
// the reference exactly: epoch is days since 1950 Jan 0.0 UTC; inclo,
// nodeo, argpo, mo are radians; noKozai is radians/minute (Kozai mean
// motion, as found on a TLE line 2); ecco is dimensionless; bstar is
// 1/earth-radii; ndot, nddot are radians/minute^2 and radians/minute^3.
//
// Mirrors sgp4init's exact call order: GetGravConst -> initl -> (if the
// orbital period is >= 225 minutes) dscom -> dpper(init='y') -> s.dsinit ->
// near-earth d2/d3/d4/t*cof correction (only when isimp==0, which deep-space
// orbits never are) -> a final Propagate(0.0) to finish populating fields
// exactly as the reference's own "propagate to zero epoch to initialize all
// others" step does. Any error from that final Propagate call (this is how
// a malformed input actually surfaces an error, e.g. immediately-decayed or
// out-of-range elements) is returned here rather than producing a silently
// broken *ElsetRec.
func NewElsetRec(
	which GravConstType,
	opsMode OpsMode,
	satNum string,
	epoch float64,
	bstar, ndot, nddot float64,
	ecco, argpo, inclo, mo, noKozai, nodeo float64,
) (*ElsetRec, error) {
	s := &ElsetRec{}

	// ----------- set all near earth variables to zero ------------
	// (already the Go zero value on a fresh ElsetRec; Method/Init/T set below)
	s.Method = 'n'

	// ----------- set all deep space variables to zero ------------
	// (already the Go zero value on a fresh ElsetRec)

	// ------------------------ earth constants -----------------------
	grav, err := GetGravConst(which)
	if err != nil {
		return nil, err
	}
	s.Tumin, s.Mus, s.Radiusearthkm, s.Xke = grav.Tumin, grav.Mus, grav.Radiusearthkm, grav.Xke
	s.J2, s.J3, s.J4, s.J3oj2 = grav.J2, grav.J3, grav.J4, grav.J3oj2

	s.Error = 0
	s.Operationmode = byte(opsMode)
	s.SatnumStr = satNum

	// sgp4fix - note the following variables are also passed directly via
	// satrec. it is possible to streamline the sgp4init call by deleting
	// the "x" variables, but the user would need to set the satrec.*
	// values first. we include the additional assignments in case
	// twoline2rv/ParseTLE is not used.
	s.Bstar = bstar
	s.Ndot = ndot
	s.Nddot = nddot
	s.Ecco = ecco
	s.Argpo = argpo
	s.Inclo = inclo
	s.Mo = mo
	s.NoKozai = noKozai
	s.Nodeo = nodeo

	// single averaged mean elements (satrec.om/ArgpM is NOT zeroed here,
	// matching the reference's own omission -- harmless, since it's
	// already the Go zero value on a fresh struct anyway)
	s.Am, s.Em, s.Im, s.NodeM, s.Mm, s.Nm = 0, 0, 0, 0, 0, 0

	// ------------------------ earth constants -----------------------
	ss := 78.0/s.Radiusearthkm + 1.0
	qzms2ttemp := (120.0 - 78.0) / s.Radiusearthkm
	qzms2t := qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp
	const x2o3 = 2.0 / 3.0

	s.Init = 'y'
	s.T = 0.0

	initRes := initl(s.Xke, s.J2, s.Ecco, epoch, s.Inclo, s.NoKozai, s.Operationmode)
	s.Method = initRes.Method
	s.Con41 = initRes.Con41
	s.Gsto = initRes.Gsto
	s.NoUnkozai = initRes.NoUnkozai
	ainv := initRes.Ainv
	ao := initRes.Ao
	con42 := initRes.Con42
	cosio := initRes.Cosio
	cosio2 := initRes.Cosio2
	eccsq := initRes.Eccsq
	omeosq := initRes.Omeosq
	posq := initRes.Posq
	rp := initRes.Rp
	rteosq := initRes.Rteosq
	sinio := initRes.Sinio
	_ = ainv

	s.A = math.Pow(s.NoUnkozai*s.Tumin, -2.0/3.0)
	s.Alta = s.A*(1.0+s.Ecco) - 1.0
	s.Altp = s.A*(1.0-s.Ecco) - 1.0
	s.Error = 0

	// sgp4fix remove sub-orbital-at-epoch check as unnecessary: the mrt
	// check in Propagate handles decaying satellites even if the starting
	// condition is below the surface of the earth (error code 5, which
	// that removed check used, is unused dead code in the reference).

	if omeosq >= 0.0 || s.NoUnkozai >= 0.0 {
		s.Isimp = 0
		if rp < 220.0/s.Radiusearthkm+1.0 {
			s.Isimp = 1
		}
		sfour := ss
		qzms24 := qzms2t
		perige := (rp - 1.0) * s.Radiusearthkm

		// - for perigees below 156 km, s and qoms2t are altered -
		if perige < 156.0 {
			sfour = perige - 78.0
			if perige < 98.0 {
				sfour = 20.0
			}
			qzms24temp := (120.0 - sfour) / s.Radiusearthkm
			qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp
			sfour = sfour/s.Radiusearthkm + 1.0
		}
		pinvsq := 1.0 / posq

		tsi := 1.0 / (ao - sfour)
		s.Eta = ao * s.Ecco * tsi
		etasq := s.Eta * s.Eta
		eeta := s.Ecco * s.Eta
		psisq := math.Abs(1.0 - etasq)
		coef := qzms24 * math.Pow(tsi, 4.0)
		coef1 := coef / math.Pow(psisq, 3.5)
		cc2 := coef1 * s.NoUnkozai * (ao*(1.0+1.5*etasq+eeta*(4.0+etasq)) +
			0.375*s.J2*tsi/psisq*s.Con41*(8.0+3.0*etasq*(8.0+etasq)))
		s.Cc1 = s.Bstar * cc2
		cc3 := 0.0
		if s.Ecco > 1.0e-4 {
			cc3 = -2.0 * coef * tsi * s.J3oj2 * s.NoUnkozai * sinio / s.Ecco
		}
		s.X1mth2 = 1.0 - cosio2
		s.Cc4 = 2.0 * s.NoUnkozai * coef1 * ao * omeosq *
			(s.Eta*(2.0+0.5*etasq) + s.Ecco*(0.5+2.0*etasq) -
				s.J2*tsi/(ao*psisq)*(-3.0*s.Con41*(1.0-2.0*eeta+etasq*(1.5-0.5*eeta))+
					0.75*s.X1mth2*(2.0*etasq-eeta*(1.0+etasq))*math.Cos(2.0*s.Argpo)))
		s.Cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75*(etasq+eeta) + eeta*etasq)
		cosio4 := cosio2 * cosio2
		temp1 := 1.5 * s.J2 * pinvsq * s.NoUnkozai
		temp2 := 0.5 * temp1 * s.J2 * pinvsq
		temp3 := -0.46875 * s.J4 * pinvsq * pinvsq * s.NoUnkozai
		s.Mdot = s.NoUnkozai + 0.5*temp1*rteosq*s.Con41 +
			0.0625*temp2*rteosq*(13.0-78.0*cosio2+137.0*cosio4)
		s.Argpdot = -0.5*temp1*con42 + 0.0625*temp2*(7.0-114.0*cosio2+395.0*cosio4) +
			temp3*(3.0-36.0*cosio2+49.0*cosio4)
		xhdot1 := -temp1 * cosio
		s.Nodedot = xhdot1 + (0.5*temp2*(4.0-19.0*cosio2)+2.0*temp3*(3.0-7.0*cosio2))*cosio
		xpidot := s.Argpdot + s.Nodedot
		s.Omgcof = s.Bstar * cc3 * math.Cos(s.Argpo)
		s.Xmcof = 0.0
		if s.Ecco > 1.0e-4 {
			s.Xmcof = -x2o3 * coef * s.Bstar / eeta
		}
		s.Nodecf = 3.5 * omeosq * xhdot1 * s.Cc1
		s.T2cof = 1.5 * s.Cc1
		// sgp4fix for divide by zero with xinco = 180 deg
		const temp4 = 1.5e-12
		if math.Abs(cosio+1.0) > temp4 {
			s.Xlcof = -0.25 * s.J3oj2 * sinio * (3.0 + 5.0*cosio) / (1.0 + cosio)
		} else {
			s.Xlcof = -0.25 * s.J3oj2 * sinio * (3.0 + 5.0*cosio) / temp4
		}
		s.Aycof = -0.5 * s.J3oj2 * sinio
		delmotemp := 1.0 + s.Eta*math.Cos(s.Mo)
		s.Delmo = delmotemp * delmotemp * delmotemp
		s.Sinmao = math.Sin(s.Mo)
		s.X7thm1 = 7.0*cosio2 - 1.0

		// --------------- deep space initialization ---------------
		if (2*pi)/s.NoUnkozai >= 225.0 {
			s.Method = 'd'
			s.Isimp = 1
			tc := 0.0
			inclm := s.Inclo

			dsc := dscom(epoch, s.Ecco, s.Argpo, tc, s.Inclo, s.Nodeo, s.NoUnkozai)
			s.E3, s.Ee2 = dsc.E3, dsc.Ee2
			s.Peo, s.Pgho, s.Pho, s.Pinco, s.Plo = dsc.Peo, dsc.Pgho, dsc.Pho, dsc.Pinco, dsc.Plo
			s.Se2, s.Se3 = dsc.Se2, dsc.Se3
			s.Sgh2, s.Sgh3, s.Sgh4 = dsc.Sgh2, dsc.Sgh3, dsc.Sgh4
			s.Sh2, s.Sh3 = dsc.Sh2, dsc.Sh3
			s.Si2, s.Si3 = dsc.Si2, dsc.Si3
			s.Sl2, s.Sl3, s.Sl4 = dsc.Sl2, dsc.Sl3, dsc.Sl4
			s.Xgh2, s.Xgh3, s.Xgh4 = dsc.Xgh2, dsc.Xgh3, dsc.Xgh4
			s.Xh2, s.Xh3 = dsc.Xh2, dsc.Xh3
			s.Xi2, s.Xi3 = dsc.Xi2, dsc.Xi3
			s.Xl2, s.Xl3, s.Xl4 = dsc.Xl2, dsc.Xl3, dsc.Xl4
			s.Zmol, s.Zmos = dsc.Zmol, dsc.Zmos

			// dpper is a no-op at this call site: satrec.init (s.Init) is
			// 'y' here, and dpper only mutates its 5 in/out values when
			// init=='n' -- so the results are provably identical to the
			// inputs and are discarded, matching what the reference itself
			// effectively does (it calls dpper here too, for structural
			// symmetry with the per-propagate-call usage, even though nothing
			// changes).
			_, _, _, _, _ = dpper(
				s.E3, s.Ee2, s.Peo, s.Pgho, s.Pho, s.Pinco, s.Plo, s.Se2, s.Se3,
				s.Sgh2, s.Sgh3, s.Sgh4, s.Sh2, s.Sh3, s.Si2, s.Si3, s.Sl2, s.Sl3, s.Sl4,
				s.T, s.Xgh2, s.Xgh3, s.Xgh4, s.Xh2, s.Xh3, s.Xi2, s.Xi3, s.Xl2, s.Xl3, s.Xl4,
				s.Zmol, s.Zmos, inclm,
				s.Init,
				s.Ecco, s.Inclo, s.Nodeo, s.Argpo, s.Mo,
				s.Operationmode,
			)

			argpm := 0.0
			nodem := 0.0
			mm := 0.0

			// dsinit's mean-element outputs (and dndt) are sgp4init's own
			// throwaway locals here too: discarded, since the final
			// Propagate(0.0) call below recomputes everything fresh from
			// the persisted ElsetRec fields dsinit itself just wrote.
			_, _, _, _, _, _, _ = s.dsinit(
				s.Xke, dsc.Cosim, dsc.Emsq, s.Argpo, dsc.S1, dsc.S2, dsc.S3, dsc.S4, dsc.S5, dsc.Sinim,
				dsc.Ss1, dsc.Ss2, dsc.Ss3, dsc.Ss4, dsc.Ss5, dsc.Sz1, dsc.Sz3, dsc.Sz11, dsc.Sz13, dsc.Sz21, dsc.Sz23,
				dsc.Sz31, dsc.Sz33, s.T, tc, s.Gsto, s.Mo, s.Mdot, s.NoUnkozai, s.Nodeo, s.Nodedot,
				xpidot, dsc.Z1, dsc.Z3, dsc.Z11, dsc.Z13, dsc.Z21, dsc.Z23, dsc.Z31, dsc.Z33, s.Ecco,
				eccsq, dsc.Em, argpm, inclm, mm, dsc.Nm, nodem,
			)
		}

		// ----------- set variables if not deep space -----------
		if s.Isimp != 1 {
			cc1sq := s.Cc1 * s.Cc1
			s.D2 = 4.0 * ao * tsi * cc1sq
			temp := s.D2 * tsi * s.Cc1 / 3.0
			s.D3 = (17.0*ao + sfour) * temp
			s.D4 = 0.5 * temp * ao * tsi * (221.0*ao + 31.0*sfour) * s.Cc1
			s.T3cof = s.D2 + 2.0*cc1sq
			s.T4cof = 0.25 * (3.0*s.D3 + s.Cc1*(12.0*s.D2+10.0*cc1sq))
			s.T5cof = 0.2 * (3.0*s.D4 + 12.0*s.Cc1*s.D3 + 6.0*s.D2*s.D2 + 15.0*cc1sq*(2.0*s.D2+cc1sq))
		}
	}

	// finally propagate to zero epoch to initialize all others.
	// sgp4fix take out check to let satellites process until they are
	// actually below earth surface.
	if _, _, err := s.Propagate(0.0); err != nil {
		return nil, err
	}

	s.Init = 'n'

	return s, nil
}

// Propagate advances the satellite by tsince minutes from epoch and returns
// the ECI position (km) and velocity (km/s) vectors, replacing sgp4() in
// SGP4.cpp (lines 1753-2046).
func (s *ElsetRec) Propagate(tsince float64) (Vector3, Vector3, error) {
	const x2o3 = 2.0 / 3.0
	const temp4 = 1.5e-12
	twopi := 2.0 * pi

	vkmpersec := s.Radiusearthkm * s.Xke / 60.0

	// --------------------- clear sgp4 error flag -----------------
	s.T = tsince
	s.Error = 0

	// ------- update for secular gravity and atmospheric drag -----
	xmdf := s.Mo + s.Mdot*s.T
	argpdf := s.Argpo + s.Argpdot*s.T
	nodedf := s.Nodeo + s.Nodedot*s.T
	argpm := argpdf
	mm := xmdf
	t2 := s.T * s.T
	nodem := nodedf + s.Nodecf*t2
	tempa := 1.0 - s.Cc1*s.T
	tempe := s.Bstar * s.Cc4 * s.T
	templ := s.T2cof * t2

	if s.Isimp != 1 {
		delomg := s.Omgcof * s.T
		delmtemp := 1.0 + s.Eta*math.Cos(xmdf)
		delm := s.Xmcof * (delmtemp*delmtemp*delmtemp - s.Delmo)
		temp := delomg + delm
		mm = xmdf + temp
		argpm = argpdf - temp
		t3 := t2 * s.T
		t4 := t3 * s.T
		tempa = tempa - s.D2*t2 - s.D3*t3 - s.D4*t4
		tempe = tempe + s.Bstar*s.Cc5*(math.Sin(mm)-s.Sinmao)
		templ = templ + s.T3cof*t3 + t4*(s.T4cof+s.T*s.T5cof)
	}

	nm := s.NoUnkozai
	em := s.Ecco
	inclm := s.Inclo
	if s.Method == 'd' {
		em, argpm, inclm, mm, nodem, nm, _ = s.dspace(em, argpm, inclm, mm, nodem, nm)
	}

	if nm <= 0.0 {
		s.Error = 2
		return Vector3{}, Vector3{}, ErrMeanMotion
	}
	am := math.Pow(s.Xke/nm, x2o3) * tempa * tempa
	nm = s.Xke / math.Pow(am, 1.5)
	em = em - tempe

	// fix tolerance for error recognition
	if em >= 1.0 || em < -0.001 {
		s.Error = 1
		return Vector3{}, Vector3{}, ErrEccentricity
	}
	// sgp4fix fix tolerance to avoid a divide by zero
	if em < 1.0e-6 {
		em = 1.0e-6
	}
	mm += s.NoUnkozai * templ
	xlm := mm + argpm + nodem
	emsq := em * em
	_ = emsq
	temp := 1.0 - emsq

	nodem = math.Mod(nodem, twopi)
	argpm = math.Mod(argpm, twopi)
	xlm = math.Mod(xlm, twopi)
	mm = math.Mod(xlm-argpm-nodem, twopi)

	// sgp4fix recover singly averaged mean elements
	s.Am = am
	s.Em = em
	s.Im = inclm
	s.NodeM = nodem
	s.ArgpM = argpm
	s.Mm = mm
	s.Nm = nm

	// ----------------- compute extra mean quantities --------------
	sinim := math.Sin(inclm)
	cosim := math.Cos(inclm)

	// -------------------- add lunar-solar periodics ---------------
	ep := em
	xincp := inclm
	argpp := argpm
	nodep := nodem
	mp := mm
	sinip := sinim
	cosip := cosim
	if s.Method == 'd' {
		ep, xincp, nodep, argpp, mp = dpper(
			s.E3, s.Ee2, s.Peo, s.Pgho, s.Pho, s.Pinco, s.Plo, s.Se2, s.Se3,
			s.Sgh2, s.Sgh3, s.Sgh4, s.Sh2, s.Sh3, s.Si2, s.Si3, s.Sl2, s.Sl3, s.Sl4,
			s.T, s.Xgh2, s.Xgh3, s.Xgh4, s.Xh2, s.Xh3, s.Xi2, s.Xi3, s.Xl2, s.Xl3, s.Xl4,
			s.Zmol, s.Zmos, s.Inclo,
			'n',
			ep, xincp, nodep, argpp, mp,
			s.Operationmode,
		)
		if xincp < 0.0 {
			xincp = -xincp
			nodep += pi
			argpp -= pi
		}
		if ep < 0.0 || ep > 1.0 {
			s.Error = 3
			return Vector3{}, Vector3{}, ErrPerturbedEccentricity
		}
	}

	if s.Method == 'd' {
		sinip = math.Sin(xincp)
		cosip = math.Cos(xincp)
		s.Aycof = -0.5 * s.J3oj2 * sinip
		// sgp4fix for divide by zero for xincp = 180 deg
		if math.Abs(cosip+1.0) > temp4 {
			s.Xlcof = -0.25 * s.J3oj2 * sinip * (3.0 + 5.0*cosip) / (1.0 + cosip)
		} else {
			s.Xlcof = -0.25 * s.J3oj2 * sinip * (3.0 + 5.0*cosip) / temp4
		}
	}
	axnl := ep * math.Cos(argpp)
	temp = 1.0 / (am * (1.0 - ep*ep))
	aynl := ep*math.Sin(argpp) + temp*s.Aycof
	xl := mp + argpp + nodep + temp*s.Xlcof*axnl

	// --------------------- solve kepler's equation -----------------
	u := math.Mod(xl-nodep, twopi)
	eo1 := u
	tem5 := 9999.9
	ktr := 1
	var sineo1, coseo1 float64
	// sgp4fix for kepler iteration: the following iteration needs better
	// limits on corrections, per the reference's own comment.
	for math.Abs(tem5) >= 1.0e-12 && ktr <= 10 {
		sineo1 = math.Sin(eo1)
		coseo1 = math.Cos(eo1)
		tem5 = 1.0 - coseo1*axnl - sineo1*aynl
		tem5 = (u - aynl*coseo1 + axnl*sineo1 - eo1) / tem5
		if math.Abs(tem5) >= 0.95 {
			if tem5 > 0.0 {
				tem5 = 0.95
			} else {
				tem5 = -0.95
			}
		}
		eo1 += tem5
		ktr++
	}

	// ------------- short period preliminary quantities -------------
	ecose := axnl*coseo1 + aynl*sineo1
	esine := axnl*sineo1 - aynl*coseo1
	el2 := axnl*axnl + aynl*aynl
	pl := am * (1.0 - el2)
	if pl < 0.0 {
		s.Error = 4
		return Vector3{}, Vector3{}, ErrSemiLatusRectum
	}

	rl := am * (1.0 - ecose)
	rdotl := math.Sqrt(am) * esine / rl
	rvdotl := math.Sqrt(pl) / rl
	betal := math.Sqrt(1.0 - el2)
	temp = esine / (1.0 + betal)
	sinu := am / rl * (sineo1 - aynl - axnl*temp)
	cosu := am / rl * (coseo1 - axnl + aynl*temp)
	su := math.Atan2(sinu, cosu)
	sin2u := (cosu + cosu) * sinu
	cos2u := 1.0 - 2.0*sinu*sinu
	temp = 1.0 / pl
	temp1 := 0.5 * s.J2 * temp
	temp2 := temp1 * temp

	// -------------- update for short period periodics --------------
	if s.Method == 'd' {
		cosisq := cosip * cosip
		s.Con41 = 3.0*cosisq - 1.0
		s.X1mth2 = 1.0 - cosisq
		s.X7thm1 = 7.0*cosisq - 1.0
	}
	mrt := rl*(1.0-1.5*temp2*betal*s.Con41) + 0.5*temp1*s.X1mth2*cos2u
	su -= 0.25 * temp2 * s.X7thm1 * sin2u
	xnode := nodep + 1.5*temp2*cosip*sin2u
	xinc := xincp + 1.5*temp2*cosip*sinip*cos2u
	mvt := rdotl - nm*temp1*s.X1mth2*sin2u/s.Xke
	rvdot := rvdotl + nm*temp1*(s.X1mth2*cos2u+1.5*s.Con41)/s.Xke

	// --------------------- orientation vectors ----------------------
	sinsu := math.Sin(su)
	cossu := math.Cos(su)
	snod := math.Sin(xnode)
	cnod := math.Cos(xnode)
	sini := math.Sin(xinc)
	cosi := math.Cos(xinc)
	xmx := -snod * cosi
	xmy := cnod * cosi
	ux := xmx*sinsu + cnod*cossu
	uy := xmy*sinsu + snod*cossu
	uz := sini * sinsu
	vx := xmx*cossu - cnod*sinsu
	vy := xmy*cossu - snod*sinsu
	vz := sini * cossu

	// --------- position and velocity (in km and km/sec) -------------
	r := Vector3{
		X: (mrt * ux) * s.Radiusearthkm,
		Y: (mrt * uy) * s.Radiusearthkm,
		Z: (mrt * uz) * s.Radiusearthkm,
	}
	v := Vector3{
		X: (mvt*ux + rvdot*vx) * vkmpersec,
		Y: (mvt*uy + rvdot*vy) * vkmpersec,
		Z: (mvt*uz + rvdot*vz) * vkmpersec,
	}

	// sgp4fix for decaying satellites: r/v are still returned (matching
	// the reference, whose r[3]/v[3] output params are already populated
	// by the time this check runs) even though an error is also returned.
	if mrt < 1.0 {
		s.Error = 6
		return r, v, ErrDecayed
	}

	return r, v, nil
}
