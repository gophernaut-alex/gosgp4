package sgp4

import (
	"errors"
	"fmt"
	"math"
	"strconv"
	"strings"
)

// ErrSGP4XPNotSupported is returned by ParseTLE for a TLE whose ephemeris
// type is 4 (SGP4-XP). The reference (twoline2rv) silently skips these
// ("sgp4fix note that the ephtype must be 0 for SGP4. SGP4-XP uses 4.");
// ParseTLE surfaces it as an explicit error instead.
var ErrSGP4XPNotSupported = errors.New("sgp4: SGP4-XP TLE (ephtype=4) not supported")

// alpha5Digit maps an uppercase letter A-Z (indexed by letter-'A') to the
// leading digit it represents in the NORAD "alpha-5" satellite numbering
// scheme for catalog numbers >= 100000. I and O are skipped (map to 0,
// unused) to avoid visual confusion with 1 and 0. Ports the alpha5[26]
// table in twoline2rv (SGP4.cpp:2254).
var alpha5Digit = [26]int{
	10, 11, 12, 13, 14, 15, 16, 17, 0, 18, 19, 20, 21, 22, 0,
	23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
}

// ParseTLE parses a two-line element set (line1, line2 -- each the full
// line, without a trailing newline) into a fully-initialized *ElsetRec,
// replacing the non-interactive path of twoline2rv in SGP4.cpp (lines
// 2181-2465).
//
// Unlike the C++, which mutates its char buffers in place to normalize
// missing decimals/blanks before sscanf, this never mutates the caller's
// strings -- Go strings are immutable, so all parsing works on fixed-column
// substrings of the input directly. The interactive manual-entry modes
// (typerun other than 'v'/'c') and the verification-format trailing
// startmfe/stopmfe/deltamin fields are both out of scope for this port (the
// latter is parsed by test code directly, from a copy of SGP4-VER.TLE, not
// by this function).
func ParseTLE(line1, line2 string, which GravConstType, opsMode OpsMode) (*ElsetRec, error) {
	if len(line1) < 68 || len(line2) < 68 {
		return nil, fmt.Errorf("sgp4: TLE line too short")
	}

	satnumStr := strings.TrimSpace(line1[2:7])
	classification := byte('U')
	if c := line1[7]; c != ' ' {
		classification = c
	}
	intldesg := strings.TrimRight(line1[9:17], " ")

	epochyr, err := strconv.Atoi(strings.TrimSpace(line1[18:20]))
	if err != nil {
		return nil, fmt.Errorf("sgp4: parsing epoch year: %w", err)
	}
	epochdays, err := strconv.ParseFloat(strings.TrimSpace(line1[20:32]), 64)
	if err != nil {
		return nil, fmt.Errorf("sgp4: parsing epoch day: %w", err)
	}

	ndot, err := strconv.ParseFloat(strings.TrimSpace(line1[33:43]), 64)
	if err != nil {
		return nil, fmt.Errorf("sgp4: parsing ndot: %w", err)
	}
	nddot, err := parseAssumedDecimal(line1[44:52])
	if err != nil {
		return nil, fmt.Errorf("sgp4: parsing nddot: %w", err)
	}
	bstar, err := parseAssumedDecimal(line1[53:61])
	if err != nil {
		return nil, fmt.Errorf("sgp4: parsing bstar: %w", err)
	}

	ephtype := 0
	if s := strings.TrimSpace(line1[62:63]); s != "" {
		ephtype, err = strconv.Atoi(s)
		if err != nil {
			return nil, fmt.Errorf("sgp4: parsing ephtype: %w", err)
		}
	}
	elnum, err := strconv.Atoi(strings.TrimSpace(line1[64:68]))
	if err != nil {
		return nil, fmt.Errorf("sgp4: parsing element set number: %w", err)
	}

	// sgp4fix note that the ephtype must be 0 for SGP4. SGP4-XP uses 4.
	if ephtype == 4 {
		return nil, ErrSGP4XPNotSupported
	}

	inclo, err := strconv.ParseFloat(strings.TrimSpace(line2[8:16]), 64)
	if err != nil {
		return nil, fmt.Errorf("sgp4: parsing inclination: %w", err)
	}
	nodeo, err := strconv.ParseFloat(strings.TrimSpace(line2[17:25]), 64)
	if err != nil {
		return nil, fmt.Errorf("sgp4: parsing raan: %w", err)
	}
	ecco, err := strconv.ParseFloat("0."+strings.TrimSpace(line2[26:33]), 64)
	if err != nil {
		return nil, fmt.Errorf("sgp4: parsing eccentricity: %w", err)
	}
	argpo, err := strconv.ParseFloat(strings.TrimSpace(line2[34:42]), 64)
	if err != nil {
		return nil, fmt.Errorf("sgp4: parsing argument of perigee: %w", err)
	}
	mo, err := strconv.ParseFloat(strings.TrimSpace(line2[43:51]), 64)
	if err != nil {
		return nil, fmt.Errorf("sgp4: parsing mean anomaly: %w", err)
	}
	noKozaiRevPerDay, err := strconv.ParseFloat(strings.TrimSpace(line2[52:63]), 64)
	if err != nil {
		return nil, fmt.Errorf("sgp4: parsing mean motion: %w", err)
	}
	revnum := 0
	if len(line2) >= 68 {
		if s := strings.TrimSpace(line2[63:68]); s != "" {
			revnum, err = strconv.Atoi(s)
			if err != nil {
				return nil, fmt.Errorf("sgp4: parsing revolution number: %w", err)
			}
		}
	}

	// ---- find no, ndot, nddot ----
	const deg2rad = pi / 180.0
	const xpdotp = 1440.0 / (2.0 * pi) // 229.1831180523293

	noKozai := noKozaiRevPerDay / xpdotp
	nddot /= xpdotp * 1440.0 * 1440.0
	ndot /= xpdotp * 1440.0

	// ---- find standard orbital elements ----
	inclo *= deg2rad
	nodeo *= deg2rad
	argpo *= deg2rad
	mo *= deg2rad

	// ---------------- temp fix for years from 1957-2056 --------------
	// correct fix will occur when year is 4-digit in tle
	var year int
	if epochyr < 57 {
		year = epochyr + 2000
	} else {
		year = epochyr + 1900
	}

	mon, day, hr, minute, sec := Days2MDHMS(year, epochdays)
	jdsatepoch, jdsatepochF := JDay(year, mon, day, hr, minute, sec)

	// sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
	epoch := (jdsatepoch + jdsatepochF) - 2433281.5

	s, err := NewElsetRec(which, opsMode, satnumStr, epoch, bstar, ndot, nddot,
		ecco, argpo, inclo, mo, noKozai, nodeo)
	if err != nil {
		return nil, err
	}

	s.Classification = classification
	s.Intldesg = intldesg
	s.Epochyr = epochyr
	s.Epochdays = epochdays
	s.Jdsatepoch = jdsatepoch
	s.JdsatepochF = jdsatepochF
	s.Ephtype = ephtype
	s.Elnum = int32(elnum)
	s.Revnum = int32(revnum)
	s.Satnum = decodeSatnum(satnumStr)

	return s, nil
}

// parseAssumedDecimal parses an 8-character TLE-style assumed-decimal
// field with a trailing signed exponent digit (e.g. " 12345-3" meaning
// 0.12345e-3, or "-12345-3" meaning -0.12345e-3), as used for the nddot and
// bstar fields.
func parseAssumedDecimal(field string) (float64, error) {
	field = strings.TrimSpace(field)
	if field == "" {
		return 0, nil
	}
	sign := 1.0
	switch field[0] {
	case '-':
		sign = -1.0
		field = field[1:]
	case '+':
		field = field[1:]
	}
	if len(field) < 2 {
		return 0, fmt.Errorf("field too short: %q", field)
	}
	mantissaStr := field[:len(field)-2]
	expStr := field[len(field)-2:]

	mantissa, err := strconv.ParseFloat("0."+mantissaStr, 64)
	if err != nil {
		return 0, err
	}
	exp, err := strconv.Atoi(expStr)
	if err != nil {
		return 0, err
	}
	return sign * mantissa * math.Pow(10, float64(exp)), nil
}

// decodeSatnum converts a 5-character catalog number field (numeric, or
// NORAD "alpha-5" for catalog numbers >= 100000) to an int, replacing the
// satnum-decoding logic in twoline2rv (SGP4.cpp:2250-2260). satnum is pure
// metadata -- it never feeds into propagation math -- so unlike the rest of
// this port, exact bug-for-bug fidelity with the reference isn't the goal
// here; this implements the alpha-5 scheme correctly (leading-letter digit
// * 10000 + the last 4 characters as a number) rather than replicating what
// looks like a redundant double substr() in the original C++.
func decodeSatnum(s string) int {
	if s == "" {
		return 0
	}
	if s[0] >= '0' && s[0] <= '9' {
		n, _ := strconv.Atoi(strings.TrimSpace(s))
		return n
	}
	letter := s[0]
	if letter >= 'a' && letter <= 'z' {
		letter -= 'a' - 'A'
	}
	if letter < 'A' || letter > 'Z' {
		return 0
	}
	digit := alpha5Digit[letter-'A']
	suffix := "0000"
	if len(s) > 1 {
		suffix = s[1:]
	}
	n, _ := strconv.Atoi(strings.TrimSpace(suffix))
	return digit*10000 + n
}
