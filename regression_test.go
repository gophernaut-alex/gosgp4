package sgp4_test

// Regression suite: propagates every test case in testdata/SGP4-VER.TLE
// (the reference verification TLE set) and checks the results against
// testdata/java_sgp4_ver.out, a golden output produced by an independent
// port of the same reference algorithm. That golden file agrees with
// testdata/testsgp4.out (a second, Matlab-derived golden file, also copied
// into testdata/ for anyone who wants to cross-check by hand) to ~9-10
// significant digits on ordinary-duration cases, so the tolerances below
// are set about an order of magnitude looser than that mutual agreement --
// tight enough to catch real bugs, loose enough to absorb ordinary
// cross-language float rounding. (One test case -- satellite 20413's
// second block, a ~9.6-year extrapolation used to probe the Lyddane-choice
// discontinuity -- diverges between the two golden files by more than that
// at the very end of the run; this is expected float-drift accumulation
// over an enormous extrapolation, not a correctness bug, which is why this
// suite is pinned to just one golden file rather than requiring agreement
// with both.)
//
// A handful of cases in SGP4-VER.TLE are deliberately malformed to exercise
// specific error paths (see the per-satellite comments in the TLE file
// itself); those are recognized by how far the golden file's own data
// extends for that satellite (a single all-zero row means even
// initialization failed; data that stops partway through the requested
// time range means Propagate is expected to start erroring at that point)
// rather than being hardcoded here by satellite number.
//
// Note on OpsModeImproved: this suite uses OpsModeImproved (not
// OpsModeAFSPC) because that's what the golden files were generated with --
// confirmed by isolating a real (if small) discrepancy for satellite 23599
// down to dpper's opsmode-gated Lyddane node-wrap fix, then verifying
// against a fresh reference-C++ compile that OpsModeAFSPC reproduces this
// Go port's result exactly while OpsModeImproved reproduces the golden
// file's. Nearly all other cases are insensitive to this choice (it only
// matters for deep-space, low-inclination orbits whose node crosses near
// zero during the propagated window).

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"regexp"
	"strconv"
	"strings"
	"testing"

	"github.com/gophernaut-alex/gosgp4"
)

type tleCase struct {
	satnumStr             string
	ephtype               int
	line1, line2          string
	startmfe, stopmfe, dt float64
}

func loadTLECases(t *testing.T, path string) []tleCase {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("opening %s: %v", path, err)
	}
	defer f.Close()

	var cases []tleCase
	var pendingLine1 string
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := sc.Text()
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		switch line[0] {
		case '1':
			pendingLine1 = line
		case '2':
			if pendingLine1 == "" {
				t.Fatalf("line 2 with no preceding line 1: %q", line)
			}
			if len(pendingLine1) < 68 {
				t.Fatalf("line 1 too short: %q", pendingLine1)
			}
			satnumStr := strings.TrimSpace(pendingLine1[2:7])
			ephtype := 0
			if s := strings.TrimSpace(pendingLine1[62:63]); s != "" {
				ephtype, err = strconv.Atoi(s)
				if err != nil {
					t.Fatalf("parsing ephtype in %q: %v", pendingLine1, err)
				}
			}

			fields := strings.Fields(line[69:])
			if len(fields) != 3 {
				t.Fatalf("expected 3 trailing verification fields on line 2, got %d: %q", len(fields), line)
			}
			startmfe, err := strconv.ParseFloat(fields[0], 64)
			if err != nil {
				t.Fatalf("parsing startmfe: %v", err)
			}
			stopmfe, err := strconv.ParseFloat(fields[1], 64)
			if err != nil {
				t.Fatalf("parsing stopmfe: %v", err)
			}
			deltamin, err := strconv.ParseFloat(fields[2], 64)
			if err != nil {
				t.Fatalf("parsing deltamin: %v", err)
			}

			cases = append(cases, tleCase{
				satnumStr: satnumStr,
				ephtype:   ephtype,
				line1:     pendingLine1,
				line2:     line[:69],
				startmfe:  startmfe,
				stopmfe:   stopmfe,
				dt:        deltamin,
			})
			pendingLine1 = ""
		}
	}
	if err := sc.Err(); err != nil {
		t.Fatalf("reading %s: %v", path, err)
	}
	return cases
}

type goldenRow struct {
	tsince float64
	r, v   sgp4.Vector3
}

var goldenHeaderRE = regexp.MustCompile(`^(\d+) xx\s*$`)

// loadGolden returns, per satellite number, the ordered list of data blocks
// for that satellite -- a plain []goldenRow per occurrence, in file order.
// A satellite number can legitimately appear more than once (e.g. 20413
// appears three times in SGP4-VER.TLE, with different start/stop/delta
// ranges each time), so blocks must be consumed in order rather than merged.
func loadGolden(t *testing.T, path string) map[string][][]goldenRow {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("opening %s: %v", path, err)
	}
	defer f.Close()

	golden := make(map[string][][]goldenRow)
	var current string
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := sc.Text()
		if m := goldenHeaderRE.FindStringSubmatch(line); m != nil {
			current = m[1]
			golden[current] = append(golden[current], nil)
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 7 || current == "" {
			continue
		}
		vals := make([]float64, 7)
		ok := true
		for i := 0; i < 7; i++ {
			vals[i], err = strconv.ParseFloat(fields[i], 64)
			if err != nil {
				ok = false
				break
			}
		}
		if !ok {
			continue
		}
		row := goldenRow{
			tsince: vals[0],
			r:      sgp4.Vector3{X: vals[1], Y: vals[2], Z: vals[3]},
			v:      sgp4.Vector3{X: vals[4], Y: vals[5], Z: vals[6]},
		}
		last := len(golden[current]) - 1
		golden[current][last] = append(golden[current][last], row)
	}
	if err := sc.Err(); err != nil {
		t.Fatalf("reading %s: %v", path, err)
	}
	return golden
}

// closeEnough compares two values with a relative tolerance, falling back
// to an absolute floor near zero (where a purely relative check is
// meaningless).
func closeEnough(got, want, relTol, absFloor float64) bool {
	diff := math.Abs(got - want)
	if diff <= absFloor {
		return true
	}
	return diff <= relTol*math.Abs(want)
}

func checkVector(t *testing.T, label string, tsince float64, got, want sgp4.Vector3, relTol, absFloor float64) {
	t.Helper()
	if !closeEnough(got.X, want.X, relTol, absFloor) ||
		!closeEnough(got.Y, want.Y, relTol, absFloor) ||
		!closeEnough(got.Z, want.Z, relTol, absFloor) {
		t.Errorf("t=%.2f %s mismatch: got (%.8f, %.8f, %.8f), want (%.8f, %.8f, %.8f)",
			tsince, label, got.X, got.Y, got.Z, want.X, want.Y, want.Z)
	}
}

// expectedTsinceSeq builds the sequence of tsince values the reference test
// driver actually propagates at: regular steps of dt from startmfe up to
// (and possibly past) stopmfe, followed by one final call at exactly
// stopmfe if the regular steps didn't already land on it exactly -- this
// "always hit the exact endpoint" tail call is a real behavior of the
// reference's own test driver, confirmed against the golden data (e.g.
// satellite 26900's 60-minute steps from 9300 land on 9300, 9360, then
// jump straight to 9400 instead of overshooting to 9420).
func expectedTsinceSeq(startmfe, stopmfe, dt float64) []float64 {
	var seq []float64
	for tsince := startmfe; tsince <= stopmfe+1e-9; tsince += dt {
		seq = append(seq, tsince)
	}
	if len(seq) == 0 || math.Abs(seq[len(seq)-1]-stopmfe) > 1e-6 {
		seq = append(seq, stopmfe)
	}
	return seq
}

func TestRegressionSGP4VER(t *testing.T) {
	const relTol = 1e-8
	const rAbsFloor = 1e-6 // km
	const vAbsFloor = 1e-9 // km/s

	cases := loadTLECases(t, "testdata/SGP4-VER.TLE")
	golden := loadGolden(t, "testdata/java_sgp4_ver.out")
	// A satellite number can appear more than once in SGP4-VER.TLE with
	// different start/stop/delta ranges each time (e.g. 20413 appears
	// three times); blockIdx tracks, per satellite number, which golden
	// block to consume next, matching file order on both sides.
	blockIdx := make(map[string]int)

	for _, tc := range cases {
		tc := tc
		if tc.ephtype == 4 {
			// matches the reference's own "only process if ephtype == 0" guard
			continue
		}

		// All SGP4-VER.TLE satellite numbers are plain 5-digit numerics (no
		// alpha-5 catalog numbers), so a direct int parse is enough to
		// match the golden file's "<satnum> xx" block headers.
		satnum, err := strconv.Atoi(tc.satnumStr)
		if err != nil {
			t.Fatalf("parsing satellite number %q: %v", tc.satnumStr, err)
		}
		key := strconv.Itoa(satnum)
		blocks := golden[key]
		idx := blockIdx[key]
		blockIdx[key] = idx + 1

		t.Run(fmt.Sprintf("sat%s/%d", key, idx), func(t *testing.T) {
			if idx >= len(blocks) {
				t.Fatalf("no golden data block #%d found for satellite %s", idx, key)
			}
			rows := blocks[idx]
			if len(rows) == 0 {
				t.Fatalf("golden data block #%d for satellite %s is empty", idx, key)
			}

			rec, err := sgp4.ParseTLE(tc.line1, tc.line2, sgp4.WGS72, sgp4.OpsModeImproved)
			if err != nil {
				if len(rows) == 1 {
					// a single all-zero golden row means the reference's own
					// init-time propagate-to-epoch failed too -- expected.
					t.Logf("expected initialization failure: %v", err)
					return
				}
				t.Fatalf("ParseTLE: unexpected error: %v", err)
			}

			// The golden file always contains one extra row at the very
			// front: the t=0 state produced as a side effect of
			// initialization (sgp4init's own internal "propagate to zero
			// epoch" step), printed in a shorter 7-field format. That row
			// is unrelated to the verification loop's own tsince=startmfe
			// unless startmfe==0, in which case the two coincide and only
			// one row is printed.
			offset := 0
			if tc.startmfe != 0 {
				offset = 1
			}

			seq := expectedTsinceSeq(tc.startmfe, tc.stopmfe, tc.dt)
			for i, tsince := range seq {
				if i+offset >= len(rows) {
					t.Logf("golden data (%d rows) shorter than requested range (%d steps) -- reference stopped early too", len(rows), len(seq))
					break
				}
				want := rows[i+offset]
				if math.Abs(tsince-want.tsince) > 1e-6 {
					t.Fatalf("row %d: tsince mismatch, computed %.6f vs golden %.6f", i, tsince, want.tsince)
				}

				r, v, err := rec.Propagate(tsince)
				if err != nil {
					if i+offset == len(rows)-1 {
						// golden data itself stops here -- the reference hit
						// the same error condition at (or after) this point.
						t.Logf("expected propagation error at t=%.2f: %v", tsince, err)
						return
					}
					t.Fatalf("Propagate(%.2f): unexpected error: %v", tsince, err)
				}

				checkVector(t, "r", tsince, r, want.r, relTol, rAbsFloor)
				checkVector(t, "v", tsince, v, want.v, relTol, vAbsFloor)
			}
		})
	}
}
