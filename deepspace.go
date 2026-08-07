package sgp4

import "math"

// DscomResult holds the deep-space common items produced by dscom, shared
// by dsinit and dpper. Field names mirror the C++ output parameter names
// in SGP4.cpp (capitalized) for easy side-by-side comparison.
type DscomResult struct {
	Snodm, Cnodm, Sinim, Cosim, Sinomm float64
	Cosomm, Day, E3, Ee2, Em           float64
	Emsq, Gam, Peo, Pgho, Pho          float64
	Pinco, Plo, Rtemsq, Se2, Se3       float64
	Sgh2, Sgh3, Sgh4, Sh2, Sh3         float64
	Si2, Si3, Sl2, Sl3, Sl4            float64
	S1, S2, S3, S4, S5                 float64
	S6, S7, Ss1, Ss2, Ss3              float64
	Ss4, Ss5, Ss6, Ss7, Sz1            float64
	Sz2, Sz3, Sz11, Sz12, Sz13         float64
	Sz21, Sz22, Sz23, Sz31, Sz32       float64
	Sz33, Xgh2, Xgh3, Xgh4, Xh2        float64
	Xh3, Xi2, Xi3, Xl2, Xl3            float64
	Xl4, Nm, Z1, Z2, Z3                float64
	Z11, Z12, Z13, Z21, Z22            float64
	Z23, Z31, Z32, Z33, Zmol           float64
	Zmos                               float64
}

// dscom computes deep-space common items (solar/lunar perturbation
// coefficients) shared by dsinit and dpper. Mirrors dscom in SGP4.cpp,
// lines 426-624. All ~78 C++ output parameters are write-only at its one
// call site (inside sgp4init), so this stays a pure function returning a
// struct rather than an (*ElsetRec) method.
func dscom(epoch, ep, argpp, tc, inclp, nodep, np float64) DscomResult {
	// -------------------------- constants -------------------------
	const (
		zes    = 0.01675
		zel    = 0.05490
		c1ss   = 2.9864797e-6
		c1l    = 4.7968065e-7
		zsinis = 0.39785416
		zcosis = 0.91744867
		zcosgs = 0.1945905
		zsings = -0.98088458
	)
	twopi := 2.0 * pi

	// --------------------- local variables ------------------------
	var a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, betasq, cc, ctem, stem float64
	var x1, x2, x3, x4, x5, x6, x7, x8, xnodce, xnoi float64
	var zcosg, zcosgl, zcosh, zcoshl, zcosi, zcosil float64
	var zsing, zsingl, zsinh, zsinhl, zsini, zsinil, zx, zy float64

	nm := np
	em := ep
	snodm := math.Sin(nodep)
	cnodm := math.Cos(nodep)
	sinomm := math.Sin(argpp)
	cosomm := math.Cos(argpp)
	sinim := math.Sin(inclp)
	cosim := math.Cos(inclp)
	emsq := em * em
	betasq = 1.0 - emsq
	rtemsq := math.Sqrt(betasq)

	// ----------------- initialize lunar solar terms ---------------
	peo := 0.0
	pinco := 0.0
	plo := 0.0
	pgho := 0.0
	pho := 0.0
	day := epoch + 18261.5 + tc/1440.0
	xnodce = math.Mod(4.5236020-9.2422029e-4*day, twopi)
	stem = math.Sin(xnodce)
	ctem = math.Cos(xnodce)
	zcosil = 0.91375164 - 0.03568096*ctem
	zsinil = math.Sqrt(1.0 - zcosil*zcosil)
	zsinhl = 0.089683511 * stem / zsinil
	zcoshl = math.Sqrt(1.0 - zsinhl*zsinhl)
	gam := 5.8351514 + 0.0019443680*day
	zx = 0.39785416 * stem / zsinil
	zy = zcoshl*ctem + 0.91744867*zsinhl*stem
	zx = math.Atan2(zx, zy)
	zx = gam + zx - xnodce
	zcosgl = math.Cos(zx)
	zsingl = math.Sin(zx)

	// ------------------------- do solar terms ----------------------
	zcosg = zcosgs
	zsing = zsings
	zcosi = zcosis
	zsini = zsinis
	zcosh = cnodm
	zsinh = snodm
	cc = c1ss
	xnoi = 1.0 / nm

	var z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33 float64
	var s1, s2, s3, s4, s5, s6, s7 float64
	var ss1, ss2, ss3, ss4, ss5, ss6, ss7 float64
	var sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33 float64

	for lsflg := 1; lsflg <= 2; lsflg++ {
		a1 = zcosg*zcosh + zsing*zcosi*zsinh
		a3 = -zsing*zcosh + zcosg*zcosi*zsinh
		a7 = -zcosg*zsinh + zsing*zcosi*zcosh
		a8 = zsing * zsini
		a9 = zsing*zsinh + zcosg*zcosi*zcosh
		a10 = zcosg * zsini
		a2 = cosim*a7 + sinim*a8
		a4 = cosim*a9 + sinim*a10
		a5 = -sinim*a7 + cosim*a8
		a6 = -sinim*a9 + cosim*a10

		x1 = a1*cosomm + a2*sinomm
		x2 = a3*cosomm + a4*sinomm
		x3 = -a1*sinomm + a2*cosomm
		x4 = -a3*sinomm + a4*cosomm
		x5 = a5 * sinomm
		x6 = a6 * sinomm
		x7 = a5 * cosomm
		x8 = a6 * cosomm

		z31 = 12.0*x1*x1 - 3.0*x3*x3
		z32 = 24.0*x1*x2 - 6.0*x3*x4
		z33 = 12.0*x2*x2 - 3.0*x4*x4
		z1 = 3.0*(a1*a1+a2*a2) + z31*emsq
		z2 = 6.0*(a1*a3+a2*a4) + z32*emsq
		z3 = 3.0*(a3*a3+a4*a4) + z33*emsq
		z11 = -6.0*a1*a5 + emsq*(-24.0*x1*x7-6.0*x3*x5)
		z12 = -6.0*(a1*a6+a3*a5) + emsq*(-24.0*(x2*x7+x1*x8)-6.0*(x3*x6+x4*x5))
		z13 = -6.0*a3*a6 + emsq*(-24.0*x2*x8-6.0*x4*x6)
		z21 = 6.0*a2*a5 + emsq*(24.0*x1*x5-6.0*x3*x7)
		z22 = 6.0*(a4*a5+a2*a6) + emsq*(24.0*(x2*x5+x1*x6)-6.0*(x4*x7+x3*x8))
		z23 = 6.0*a4*a6 + emsq*(24.0*x2*x6-6.0*x4*x8)
		z1 = z1 + z1 + betasq*z31
		z2 = z2 + z2 + betasq*z32
		z3 = z3 + z3 + betasq*z33
		s3 = cc * xnoi
		s2 = -0.5 * s3 / rtemsq
		s4 = s3 * rtemsq
		s1 = -15.0 * em * s4
		s5 = x1*x3 + x2*x4
		s6 = x2*x3 + x1*x4
		s7 = x2*x4 - x1*x3

		// ----------------------- do lunar terms -------------------
		if lsflg == 1 {
			ss1 = s1
			ss2 = s2
			ss3 = s3
			ss4 = s4
			ss5 = s5
			ss6 = s6
			ss7 = s7
			sz1 = z1
			sz2 = z2
			sz3 = z3
			sz11 = z11
			sz12 = z12
			sz13 = z13
			sz21 = z21
			sz22 = z22
			sz23 = z23
			sz31 = z31
			sz32 = z32
			sz33 = z33
			zcosg = zcosgl
			zsing = zsingl
			zcosi = zcosil
			zsini = zsinil
			zcosh = zcoshl*cnodm + zsinhl*snodm
			zsinh = snodm*zcoshl - cnodm*zsinhl
			cc = c1l
		}
	}

	zmol := math.Mod(4.7199672+0.22997150*day-gam, twopi)
	zmos := math.Mod(6.2565837+0.017201977*day, twopi)

	// ------------------------ do solar terms -------------------------
	se2 := 2.0 * ss1 * ss6
	se3 := 2.0 * ss1 * ss7
	si2 := 2.0 * ss2 * sz12
	si3 := 2.0 * ss2 * (sz13 - sz11)
	sl2 := -2.0 * ss3 * sz2
	sl3 := -2.0 * ss3 * (sz3 - sz1)
	sl4 := -2.0 * ss3 * (-21.0 - 9.0*emsq) * zes
	sgh2 := 2.0 * ss4 * sz32
	sgh3 := 2.0 * ss4 * (sz33 - sz31)
	sgh4 := -18.0 * ss4 * zes
	sh2 := -2.0 * ss2 * sz22
	sh3 := -2.0 * ss2 * (sz23 - sz21)

	// ------------------------ do lunar terms -------------------------
	ee2 := 2.0 * s1 * s6
	e3 := 2.0 * s1 * s7
	xi2 := 2.0 * s2 * z12
	xi3 := 2.0 * s2 * (z13 - z11)
	xl2 := -2.0 * s3 * z2
	xl3 := -2.0 * s3 * (z3 - z1)
	xl4 := -2.0 * s3 * (-21.0 - 9.0*emsq) * zel
	xgh2 := 2.0 * s4 * z32
	xgh3 := 2.0 * s4 * (z33 - z31)
	xgh4 := -18.0 * s4 * zel
	xh2 := -2.0 * s2 * z22
	xh3 := -2.0 * s2 * (z23 - z21)

	return DscomResult{
		Snodm: snodm, Cnodm: cnodm, Sinim: sinim, Cosim: cosim, Sinomm: sinomm,
		Cosomm: cosomm, Day: day, E3: e3, Ee2: ee2, Em: em,
		Emsq: emsq, Gam: gam, Peo: peo, Pgho: pgho, Pho: pho,
		Pinco: pinco, Plo: plo, Rtemsq: rtemsq, Se2: se2, Se3: se3,
		Sgh2: sgh2, Sgh3: sgh3, Sgh4: sgh4, Sh2: sh2, Sh3: sh3,
		Si2: si2, Si3: si3, Sl2: sl2, Sl3: sl3, Sl4: sl4,
		S1: s1, S2: s2, S3: s3, S4: s4, S5: s5,
		S6: s6, S7: s7, Ss1: ss1, Ss2: ss2, Ss3: ss3,
		Ss4: ss4, Ss5: ss5, Ss6: ss6, Ss7: ss7, Sz1: sz1,
		Sz2: sz2, Sz3: sz3, Sz11: sz11, Sz12: sz12, Sz13: sz13,
		Sz21: sz21, Sz22: sz22, Sz23: sz23, Sz31: sz31, Sz32: sz32,
		Sz33: sz33, Xgh2: xgh2, Xgh3: xgh3, Xgh4: xgh4, Xh2: xh2,
		Xh3: xh3, Xi2: xi2, Xi3: xi3, Xl2: xl2, Xl3: xl3,
		Xl4: xl4, Nm: nm, Z1: z1, Z2: z2, Z3: z3,
		Z11: z11, Z12: z12, Z13: z13, Z21: z21, Z22: z22,
		Z23: z23, Z31: z31, Z32: z32, Z33: z33, Zmol: zmol,
		Zmos: zmos,
	}
}

// dpper applies deep-space lunar/solar long-period periodic corrections to
// the mean elements ep, inclp, nodep, argpp, mp, replacing dpper in
// SGP4.cpp (lines 231-356). At both of its real call sites (inside
// NewElsetRec and inside Propagate) these 5 values are transient locals of
// the caller, never ElsetRec fields directly, so dpper stays a plain
// function (value in, value out) rather than an (*ElsetRec) method. When
// init=='y' the 5 elements are returned unchanged -- periodics are zero at
// epoch by design. inclo is accepted for signature fidelity with the
// reference but is unused in this body (the original strn3 branch that used
// it is commented out upstream in favor of the gsfc perturbed-inclination
// check on inclp).
func dpper(e3, ee2, peo, pgho, pho, pinco, plo, se2, se3, sgh2, sgh3, sgh4,
	sh2, sh3, si2, si3, sl2, sl3, sl4, t, xgh2, xgh3, xgh4, xh2, xh3,
	xi2, xi3, xl2, xl3, xl4, zmol, zmos, inclo float64,
	init byte,
	ep, inclp, nodep, argpp, mp float64,
	opsmode byte,
) (epOut, inclpOut, nodepOut, argppOut, mpOut float64) {
	twopi := 2.0 * pi

	const (
		zns = 1.19459e-5
		zes = 0.01675
		znl = 1.5835218e-4
		zel = 0.05490
	)

	// --------------- calculate time varying periodics -----------
	zm := zmos + zns*t
	if init == 'y' {
		zm = zmos
	}
	zf := zm + 2.0*zes*math.Sin(zm)
	sinzf := math.Sin(zf)
	f2 := 0.5*sinzf*sinzf - 0.25
	f3 := -0.5 * sinzf * math.Cos(zf)
	ses := se2*f2 + se3*f3
	sis := si2*f2 + si3*f3
	sls := sl2*f2 + sl3*f3 + sl4*sinzf
	sghs := sgh2*f2 + sgh3*f3 + sgh4*sinzf
	shs := sh2*f2 + sh3*f3
	zm = zmol + znl*t
	if init == 'y' {
		zm = zmol
	}
	zf = zm + 2.0*zel*math.Sin(zm)
	sinzf = math.Sin(zf)
	f2 = 0.5*sinzf*sinzf - 0.25
	f3 = -0.5 * sinzf * math.Cos(zf)
	sel := ee2*f2 + e3*f3
	sil := xi2*f2 + xi3*f3
	sll := xl2*f2 + xl3*f3 + xl4*sinzf
	sghl := xgh2*f2 + xgh3*f3 + xgh4*sinzf
	shll := xh2*f2 + xh3*f3
	pe := ses + sel
	pinc := sis + sil
	pl := sls + sll
	pgh := sghs + sghl
	ph := shs + shll

	if init == 'n' {
		pe -= peo
		pinc -= pinco
		pl -= plo
		pgh -= pgho
		ph -= pho
		inclp += pinc
		ep += pe
		sinip := math.Sin(inclp)
		cosip := math.Cos(inclp)

		// ----------------- apply periodics directly ------------
		//  sgp4fix for lyddane choice: strn3 used original inclination
		//  (technically feasible); gsfc used perturbed inclination (also
		//  technically feasible) -- this uses the gsfc/perturbed check.
		if inclp >= 0.2 {
			ph /= sinip
			pgh -= cosip * ph
			argpp += pgh
			nodep += ph
			mp += pl
		} else {
			// ---- apply periodics with lyddane modification ----
			sinop := math.Sin(nodep)
			cosop := math.Cos(nodep)
			alfdp := sinip * sinop
			betdp := sinip * cosop
			dalf := ph*cosop + pinc*cosip*sinop
			dbet := -ph*sinop + pinc*cosip*cosop
			alfdp += dalf
			betdp += dbet
			nodep = math.Mod(nodep, twopi)
			if nodep < 0.0 && opsmode == 'a' {
				nodep += twopi
			}
			xls := mp + argpp + cosip*nodep
			dls := pl + pgh - pinc*nodep*sinip
			xls += dls
			xnoh := nodep
			nodep = math.Atan2(alfdp, betdp)
			if nodep < 0.0 && opsmode == 'a' {
				nodep += twopi
			}
			if math.Abs(xnoh-nodep) > pi {
				if nodep < xnoh {
					nodep += twopi
				} else {
					nodep -= twopi
				}
			}
			mp += pl
			argpp = xls - mp - cosip*nodep
		}
	}

	return ep, inclp, nodep, argpp, mp
}

// dsinit computes deep-space secular rates and (conditionally, when a 12hr
// or 24hr geopotential resonance is detected) initializes the resonance
// integrator terms, replacing dsinit in SGP4.cpp (lines 707-931).
//
// Confirmed against its real call site (SGP4.cpp:1617-1632, inside
// sgp4init/NewElsetRec): em, argpm, inclm, mm, nm, nodem, and dndt are all
// transient locals of the caller there (reset/discarded around the call,
// never ElsetRec fields), so they're ordinary value-in/value-out here. But
// irez, atime, d2201..d5433, dedt, didt, dmdt, dnodt, domdt, del1-3, xfact,
// xlamo, xli, xni ARE passed as real satrec fields at that call site, so
// dsinit mutates them directly on the receiver.
func (s *ElsetRec) dsinit(
	xke, cosim, emsq, argpo, s1, s2, s3, s4, s5, sinim,
	ss1, ss2, ss3, ss4, ss5, sz1, sz3, sz11, sz13, sz21, sz23,
	sz31, sz33, t, tc, gsto, mo, mdot, no, nodeo, nodedot,
	xpidot, z1, z3, z11, z13, z21, z23, z31, z33, ecco,
	eccsq, em, argpm, inclm, mm, nm, nodem float64,
) (emOut, argpmOut, inclmOut, mmOut, nmOut, nodemOut, dndt float64) {
	twopi := 2.0 * pi

	const (
		q22    = 1.7891679e-6
		q31    = 2.1460748e-6
		q33    = 2.2123015e-7
		root22 = 1.7891679e-6
		root44 = 7.3636953e-9
		root54 = 2.1765803e-9
		rptim  = 4.37526908801129966e-3 // equates to 7.29211514668855e-5 rad/sec
		root32 = 3.7393792e-7
		root52 = 1.1428639e-7
		x2o3   = 2.0 / 3.0
		znl    = 1.5835218e-4
		zns    = 1.19459e-5
	)

	// -------------------- deep space initialization ------------
	irez := 0
	if nm < 0.0052359877 && nm > 0.0034906585 {
		irez = 1
	}
	if nm >= 8.26e-3 && nm <= 9.24e-3 && em >= 0.5 {
		irez = 2
	}

	// ------------------------ do solar terms -------------------
	ses := ss1 * zns * ss5
	sis := ss2 * zns * (sz11 + sz13)
	sls := -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0*emsq)
	sghs := ss4 * zns * (sz31 + sz33 - 6.0)
	shs := -zns * ss2 * (sz21 + sz23)
	// sgp4fix for 180 deg incl
	if inclm < 5.2359877e-2 || inclm > pi-5.2359877e-2 {
		shs = 0.0
	}
	if sinim != 0.0 {
		shs /= sinim
	}
	sgs := sghs - cosim*shs

	// ------------------------- do lunar terms ------------------
	dedt := ses + s1*znl*s5
	didt := sis + s2*znl*(z11+z13)
	dmdt := sls - znl*s3*(z1+z3-14.0-6.0*emsq)
	sghl := s4 * znl * (z31 + z33 - 6.0)
	shll := -znl * s2 * (z21 + z23)
	// sgp4fix for 180 deg incl
	if inclm < 5.2359877e-2 || inclm > pi-5.2359877e-2 {
		shll = 0.0
	}
	domdt := sgs + sghl
	dnodt := shs
	if sinim != 0.0 {
		domdt -= cosim / sinim * shll
		dnodt += shll / sinim
	}

	// ----------- calculate deep space resonance effects --------
	dndt = 0.0
	theta := math.Mod(gsto+tc*rptim, twopi)
	em += dedt * t
	inclm += didt * t
	argpm += domdt * t
	nodem += dnodt * t
	mm += dmdt * t
	// sgp4fix for negative inclinations -- the corresponding if statement
	// is commented out upstream too; preserved here only as a comment for
	// fidelity, not applied.
	// if inclm < 0.0 { inclm = -inclm; argpm -= pi; nodem += pi }

	var atime, d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 float64
	var del1, del2, del3, xfact, xlamo, xli, xni float64

	// -------------- initialize the resonance terms -------------
	if irez != 0 {
		aonv := math.Pow(nm/xke, x2o3)

		// ---------- geopotential resonance for 12 hour orbits ------
		if irez == 2 {
			cosisq := cosim * cosim
			emo := em
			em = ecco
			emsqo := emsq
			emsq = eccsq
			eoc := em * emsq
			g201 := -0.306 - (em-0.64)*0.440

			var g211, g310, g322, g410, g422, g520 float64
			if em <= 0.65 {
				g211 = 3.616 - 13.2470*em + 16.2900*emsq
				g310 = -19.302 + 117.3900*em - 228.4190*emsq + 156.5910*eoc
				g322 = -18.9068 + 109.7927*em - 214.6334*emsq + 146.5816*eoc
				g410 = -41.122 + 242.6940*em - 471.0940*emsq + 313.9530*eoc
				g422 = -146.407 + 841.8800*em - 1629.014*emsq + 1083.4350*eoc
				g520 = -532.114 + 3017.977*em - 5740.032*emsq + 3708.2760*eoc
			} else {
				g211 = -72.099 + 331.819*em - 508.738*emsq + 266.724*eoc
				g310 = -346.844 + 1582.851*em - 2415.925*emsq + 1246.113*eoc
				g322 = -342.585 + 1554.908*em - 2366.899*emsq + 1215.972*eoc
				g410 = -1052.797 + 4758.686*em - 7193.992*emsq + 3651.957*eoc
				g422 = -3581.690 + 16178.110*em - 24462.770*emsq + 12422.520*eoc
				if em > 0.715 {
					g520 = -5149.66 + 29936.92*em - 54087.36*emsq + 31324.56*eoc
				} else {
					g520 = 1464.74 - 4664.75*em + 3763.64*emsq
				}
			}

			var g533, g521, g532 float64
			if em < 0.7 {
				g533 = -919.22770 + 4988.6100*em - 9064.7700*emsq + 5542.21*eoc
				g521 = -822.71072 + 4568.6173*em - 8491.4146*emsq + 5337.524*eoc
				g532 = -853.66600 + 4690.2500*em - 8624.7700*emsq + 5341.4*eoc
			} else {
				g533 = -37995.780 + 161616.52*em - 229838.20*emsq + 109377.94*eoc
				g521 = -51752.104 + 218913.95*em - 309468.16*emsq + 146349.42*eoc
				g532 = -40023.880 + 170470.89*em - 242699.48*emsq + 115605.82*eoc
			}

			sini2 := sinim * sinim
			f220 := 0.75 * (1.0 + 2.0*cosim + cosisq)
			f221 := 1.5 * sini2
			f321 := 1.875 * sinim * (1.0 - 2.0*cosim - 3.0*cosisq)
			f322 := -1.875 * sinim * (1.0 + 2.0*cosim - 3.0*cosisq)
			f441 := 35.0 * sini2 * f220
			f442 := 39.3750 * sini2 * sini2
			f522 := 9.84375 * sinim * (sini2*(1.0-2.0*cosim-5.0*cosisq) +
				0.33333333*(-2.0+4.0*cosim+6.0*cosisq))
			f523 := sinim * (4.92187512*sini2*(-2.0-4.0*cosim+10.0*cosisq) +
				6.56250012*(1.0+2.0*cosim-3.0*cosisq))
			f542 := 29.53125 * sinim * (2.0 - 8.0*cosim + cosisq*(-12.0+8.0*cosim+10.0*cosisq))
			f543 := 29.53125 * sinim * (-2.0 - 8.0*cosim + cosisq*(12.0+8.0*cosim-10.0*cosisq))
			xno2 := nm * nm
			ainv2 := aonv * aonv
			temp1 := 3.0 * xno2 * ainv2
			temp := temp1 * root22
			d2201 = temp * f220 * g201
			d2211 = temp * f221 * g211
			temp1 *= aonv
			temp = temp1 * root32
			d3210 = temp * f321 * g310
			d3222 = temp * f322 * g322
			temp1 *= aonv
			temp = 2.0 * temp1 * root44
			d4410 = temp * f441 * g410
			d4422 = temp * f442 * g422
			temp1 *= aonv
			temp = temp1 * root52
			d5220 = temp * f522 * g520
			d5232 = temp * f523 * g532
			temp = 2.0 * temp1 * root54
			d5421 = temp * f542 * g521
			d5433 = temp * f543 * g533
			xlamo = math.Mod(mo+nodeo+nodeo-theta-theta, twopi)
			xfact = mdot + dmdt + 2.0*(nodedot+dnodt-rptim) - no
			em = emo
			emsq = emsqo
		}

		// ---------------- synchronous resonance terms --------------
		if irez == 1 {
			g200 := 1.0 + emsq*(-2.5+0.8125*emsq)
			g310 := 1.0 + 2.0*emsq
			g300 := 1.0 + emsq*(-6.0+6.60937*emsq)
			f220 := 0.75 * (1.0 + cosim) * (1.0 + cosim)
			f311 := 0.9375*sinim*sinim*(1.0+3.0*cosim) - 0.75*(1.0+cosim)
			f330 := 1.0 + cosim
			f330 = 1.875 * f330 * f330 * f330
			del1 = 3.0 * nm * nm * aonv * aonv
			del2 = 2.0 * del1 * f220 * g200 * q22
			del3 = 3.0 * del1 * f330 * g300 * q33 * aonv
			del1 = del1 * f311 * g310 * q31 * aonv
			xlamo = math.Mod(mo+nodeo+argpo-theta, twopi)
			xfact = mdot + xpidot - rptim + dmdt + domdt + dnodt - no
		}

		// ------------ for sgp4, initialize the integrator ----------
		xli = xlamo
		xni = no
		atime = 0.0
		nm = no + dndt
	}

	s.Irez = irez
	s.Atime = atime
	s.D2201, s.D2211 = d2201, d2211
	s.D3210, s.D3222 = d3210, d3222
	s.D4410, s.D4422 = d4410, d4422
	s.D5220, s.D5232 = d5220, d5232
	s.D5421, s.D5433 = d5421, d5433
	s.Dedt, s.Didt, s.Dmdt = dedt, didt, dmdt
	s.Dnodt, s.Domdt = dnodt, domdt
	s.Del1, s.Del2, s.Del3 = del1, del2, del3
	s.Xfact, s.Xlamo, s.Xli, s.Xni = xfact, xlamo, xli, xni

	return em, argpm, inclm, mm, nm, nodem, dndt
}

// dspace advances deep-space secular/resonance terms to the current time,
// replacing dspace in SGP4.cpp (lines 1006-1153).
//
// Confirmed against its one real call site (SGP4.cpp:1827-1841, inside
// sgp4()/Propagate): atime, xli, xni ARE satrec fields there (genuine
// cross-call persistent deep-space integrator state), so they're mutated
// directly on the receiver; em, argpm, inclm, mm, nodem, nm, and dndt are
// transient locals of the caller (seeded from satrec fields moments
// earlier, written back to ElsetRec only much later after the Kepler
// solve), so they stay ordinary value in/out.
//
// The other C++ inputs to dspace (irez, d2201..d5433, dedt, del1-3, didt,
// dmdt, dnodt, domdt, argpo, argpdot, gsto, xfact, xlamo, no) are all
// ElsetRec fields dspace only reads, so this reads them directly off the
// receiver instead of duplicating them as parameters. tc is dropped too:
// at dspace's one call site tc is always computed as satrec.t immediately
// beforehand, and Propagate sets s.T = tsince before calling this, so tc
// and s.T are always the same value here.
func (s *ElsetRec) dspace(em, argpm, inclm, mm, nodem, nm float64) (emOut, argpmOut, inclmOut, mmOut, nodemOut, nmOut, dndt float64) {
	twopi := 2.0 * pi
	tc := s.T

	const (
		fasx2 = 0.13130908
		fasx4 = 2.8843198
		fasx6 = 0.37448087
		g22   = 5.7686396
		g32   = 0.95240898
		g44   = 1.8014998
		g52   = 1.0508330
		g54   = 4.4108898
		rptim = 4.37526908801129966e-3 // equates to 7.29211514668855e-5 rad/sec
		stepp = 720.0
		stepn = -720.0
		step2 = 259200.0
	)

	// ----------- calculate deep space resonance effects -----------
	dndt = 0.0
	theta := math.Mod(s.Gsto+tc*rptim, twopi)
	em += s.Dedt * s.T
	inclm += s.Didt * s.T
	argpm += s.Domdt * s.T
	nodem += s.Dnodt * s.T
	mm += s.Dmdt * s.T

	// - update resonances: numerical (euler-maclaurin) integration -
	ft := 0.0
	if s.Irez != 0 {
		// sgp4fix streamline check
		if s.Atime == 0.0 || s.T*s.Atime <= 0.0 || math.Abs(s.T) < math.Abs(s.Atime) {
			s.Atime = 0.0
			s.Xni = s.NoUnkozai
			s.Xli = s.Xlamo
		}
		var delt float64
		if s.T > 0.0 {
			delt = stepp
		} else {
			delt = stepn
		}

		var xndt, xldot, xnddt float64
		for {
			// ------------------- dot terms calculated -------------
			if s.Irez != 2 {
				// ----------- near-synchronous resonance terms -------
				xndt = s.Del1*math.Sin(s.Xli-fasx2) + s.Del2*math.Sin(2.0*(s.Xli-fasx4)) +
					s.Del3*math.Sin(3.0*(s.Xli-fasx6))
				xldot = s.Xni + s.Xfact
				xnddt = s.Del1*math.Cos(s.Xli-fasx2) +
					2.0*s.Del2*math.Cos(2.0*(s.Xli-fasx4)) +
					3.0*s.Del3*math.Cos(3.0*(s.Xli-fasx6))
				xnddt *= xldot
			} else {
				// --------- near-half-day resonance terms --------
				xomi := s.Argpo + s.Argpdot*s.Atime
				x2omi := xomi + xomi
				x2li := s.Xli + s.Xli
				xndt = s.D2201*math.Sin(x2omi+s.Xli-g22) + s.D2211*math.Sin(s.Xli-g22) +
					s.D3210*math.Sin(xomi+s.Xli-g32) + s.D3222*math.Sin(-xomi+s.Xli-g32) +
					s.D4410*math.Sin(x2omi+x2li-g44) + s.D4422*math.Sin(x2li-g44) +
					s.D5220*math.Sin(xomi+s.Xli-g52) + s.D5232*math.Sin(-xomi+s.Xli-g52) +
					s.D5421*math.Sin(xomi+x2li-g54) + s.D5433*math.Sin(-xomi+x2li-g54)
				xldot = s.Xni + s.Xfact
				xnddt = s.D2201*math.Cos(x2omi+s.Xli-g22) + s.D2211*math.Cos(s.Xli-g22) +
					s.D3210*math.Cos(xomi+s.Xli-g32) + s.D3222*math.Cos(-xomi+s.Xli-g32) +
					s.D5220*math.Cos(xomi+s.Xli-g52) + s.D5232*math.Cos(-xomi+s.Xli-g52) +
					2.0*(s.D4410*math.Cos(x2omi+x2li-g44)+
						s.D4422*math.Cos(x2li-g44)+s.D5421*math.Cos(xomi+x2li-g54)+
						s.D5433*math.Cos(-xomi+x2li-g54))
				xnddt *= xldot
			}

			// ----------------------- integrator -------------------
			if math.Abs(s.T-s.Atime) >= stepp {
				// full step: keep integrating
				s.Xli += xldot*delt + xndt*step2
				s.Xni += xndt*delt + xnddt*step2
				s.Atime += delt
				continue
			}
			// exit here: final partial step
			ft = s.T - s.Atime
			break
		}

		nm = s.Xni + xndt*ft + xnddt*ft*ft*0.5
		xl := s.Xli + xldot*ft + xndt*ft*ft*0.5
		if s.Irez != 1 {
			mm = xl - 2.0*nodem + 2.0*theta
			dndt = nm - s.NoUnkozai
		} else {
			mm = xl - nodem - argpm + theta
			dndt = nm - s.NoUnkozai
		}
		nm = s.NoUnkozai + dndt
	}

	return em, argpm, inclm, mm, nodem, nm, dndt
}
