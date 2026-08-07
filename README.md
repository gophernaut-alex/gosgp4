# gosgp4

A Go port of the reference **SGP4/SDP4** analytical satellite propagator.

## Provenance

This is a line-for-line port of the C++ reference implementation from
[CelesTrak/fundamentals-of-astrodynamics](https://github.com/CelesTrak/fundamentals-of-astrodynamics)
(`software/cpp/SGP4/SGP4/SGP4.cpp` + `SGP4.h`), itself the canonical
distribution of David Vallado's algorithm described in:

> Vallado, D., Crawford, P., Hujsak, R., Kelso, T.S., "Revisiting
> Spacetrack Report #3," AIAA 2006-6753, presented at the AIAA/AAS
> Astrodynamics Specialist Conference, August 2006.

A copy of that reference source is vendored under `reference/vallado/` in
this repository for side-by-side comparison during development; it is not
part of the Go module itself.

An independently-distributed, older Vallado release (`SGP4.cpp`/`SGP4.h`
from the 2020-era SGP4DC orbit-determination package, vendored for
comparison under `reference2/sgp4dc/`) was diffed against the above and
found byte-identical in the core propagator (`sgp4init`, `sgp4`, and all
deep-space routines) — corroborating that this algorithm has been the same,
consistently-authored Vallado code for years.

Every internal function here carries a doc comment pointing back to its
originating C++ function name and line range, and every exported struct
field mirrors the C++ variable name it replaces (capitalized) — the goal
throughout was traceability against the original, not just a working
result.

The vector-math and orbital-element utilities in `vector.go`/`timeutil.go`
(`RVToCOE`, `Angle`, `NewtonNu`, `JDay`, ...) aren't used by the propagator
itself, but are included for parity with the reference's own public API
surface.

## License

**MIT** (see `LICENSE`) for this package's own Go source. No copyleft
obligations — you can use this in closed-source or commercial projects
without having to open-source anything of your own.

The underlying SGP4/SDP4 algorithm this package ports is David Vallado's;
CelesTrak's own FAQ for the paper it's based on states the algorithm/code
may be used "for any purpose — personal or commercial — as you wish,"
conditioned on citing the source (see
https://celestrak.org/publications/AIAA/2006-6753/faq.php). See `NOTICE`
for that citation.

## Install

```
go get github.com/gophernaut-alex/gosgp4
```

## Usage

### From a TLE

```go
package main

import (
	"fmt"
	sgp4 "github.com/gophernaut-alex/gosgp4"
)

func main() {
	line1 := "1 25544U 98067A   24025.51782528  .00016717  00000-0  30412-3 0  9990"
	line2 := "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.49560856433624"

	sat, err := sgp4.ParseTLE(line1, line2, sgp4.WGS72, sgp4.OpsModeImproved)
	if err != nil {
		panic(err)
	}

	// Position (km) and velocity (km/s) in the TEME frame, 90 minutes
	// after epoch.
	r, v, err := sat.Propagate(90.0)
	if err != nil {
		panic(err)
	}
	fmt.Printf("r = %+v km\n", r)
	fmt.Printf("v = %+v km/s\n", v)
}
```

### From raw mean elements

```go
sat, err := sgp4.NewElsetRec(
	sgp4.WGS72, sgp4.OpsModeImproved,
	"25544",     // satellite number
	epochDays,   // days since 1950 Jan 0.0 UTC
	bstar, ndot, nddot,
	ecco, argpo, inclo, mo, noKozai, nodeo, // all angles in radians
)
```

### Errors

`Propagate` returns one of the sentinel errors in `errors.go`
(`ErrMeanMotion`, `ErrEccentricity`, `ErrPerturbedEccentricity`,
`ErrSemiLatusRectum`, `ErrDecayed`) when the propagator hits one of the
reference's own failure conditions — check with `errors.Is`.

### `OpsMode`

Pass `OpsModeAFSPC` or `OpsModeImproved` depending on which sidereal-time/
node-wrap convention you need to match. They only produce different results
for deep-space, low-inclination orbits whose ascending node crosses near
zero during the propagated window — see the note atop `regression_test.go`
for how this was tracked down. If you don't have a specific reason to
match AFSPC's legacy behavior, use `OpsModeImproved`.

## Architecture

Each file maps to specific function(s) in the C++ reference
(`reference/vallado/software/cpp/SGP4/SGP4/SGP4.cpp` unless noted):

| Go file | C++ source | Notes |
|---|---|---|
| `elsetrec.go` | `elsetrec` struct (`SGP4.h`) | Field-for-field port; exported names capitalize the C++ ones |
| `gravconst.go` | `gravconsttype`, `getgravconst` | `GetGravConst` returns an error instead of the reference's silent no-op on an invalid model |
| `initl.go` | `initl` | Plain function returning `InitlResult` (mixed field/local outputs, like `dscom`) |
| `deepspace.go` | `dscom`, `dpper`, `dsinit`, `dspace` | `dscom`/`dpper` are pure functions; `dsinit`/`dspace` are `(*ElsetRec)` methods — see the in/out translation rule below |
| `propagate.go` | `sgp4init`, `sgp4` | `NewElsetRec` / `(*ElsetRec) Propagate` |
| `tle.go` | `twoline2rv` (non-interactive path only) | `ParseTLE`; never mutates caller strings; interactive manual-entry modes are out of scope |
| `vector.go` | `mag/cross/dot/angle/sgn/newtonnu/rv2coe_SGP4` | `Vector3` type + methods; `asinh_SGP4` replaced by stdlib `math.Asinh` |
| `timeutil.go` | `jday/days2mdhms/invjday/gstime_SGP4` | Ported faithfully (not rebuilt on `time.Time`) to preserve the split whole/fractional-Julian-date precision the propagator depends on |
| `errors.go` | `satrec.error` codes | Sentinel errors for codes 1-4 and 6 (code 5 is dead code upstream) |

**The in/out parameter translation rule:** the C++ leans on `double&`
reference parameters to fake multi-value returns, which don't map 1:1 onto
idiomatic Go. Each function here was translated according to what its
*real call site* in the reference actually does with each reference
parameter (not just what the signature suggests):

- If a parameter is passed as `&satrec.xxx` at its real call site (genuinely
  persistent satellite state, e.g. `dsinit`'s resonance coefficients or
  `dspace`'s `atime`/`xli`/`xni` integrator state), the Go function is an
  unexported method on `*ElsetRec` that mutates that field directly.
- If it's a transient scratch value never stored in `elsetrec` (e.g.
  `dpper`'s 5 in/out mean-element params, or nearly all of `dscom`'s ~80
  outputs), the Go function stays a plain function using ordinary
  value-in/value-out returns, packaged into a named result struct
  (`DscomResult`, `InitlResult`) once there are more than a handful of
  return values.

## Testing

```
go test ./...
```

- **Unit tests** (`initl_test.go`, `deepspace_test.go`) check `initl`,
  `dscom`, `dpper`, `dsinit`, and `dspace` directly against hand-derived
  reference vectors, covering all three deep-space resonance
  classifications (`irez` = 0/1/2). Those vectors were generated by
  compiling a patched copy of the upstream `SGP4.cpp` with a small scratch
  harness that calls each function directly and prints full-precision
  output — the harness itself isn't part of this module; regenerate it
  from `reference/vallado/software/cpp/SGP4/SGP4/SGP4.cpp` if these
  functions ever change and the fixtures need updating.
- **Constant tests** (`gravconst_test.go`) check every field of all three
  gravity models (`WGS72Old`, `WGS72`, `WGS84`) against the literal values in
  `getgravconst` (`SGP4.cpp`). These are the only tests that exercise
  `WGS72Old`/`WGS84` at all — every other test in this repo uses `WGS72` —
  which is exactly how a real bug (`WGS72Old`'s `Tumin` was transcribed as
  `1.0/60.0` instead of `1.0/xke`) went undetected until this file was added.
- **Regression tests** (`regression_test.go`) propagate every case in
  `testdata/SGP4-VER.TLE` (the reference's own verification TLE set — 34
  satellites covering near-earth, all three deep-space resonance types,
  deliberately malformed error-code cases, and a ~9.6-year long-duration
  edge case) and compare against `testdata/java_sgp4_ver.out`, a golden
  output from an independent port of the same reference algorithm.
  `testdata/testsgp4.out` (Matlab-derived) is also included for anyone who
  wants to cross-check by hand.
