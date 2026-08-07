package sgp4

import "math"

// JDay converts a calendar date/time to a Julian date split into a whole-day
// part jd and a fractional-day part jdFrac, replacing jday_SGP4. The split
// (rather than a single float64) preserves sub-second precision that would
// otherwise be lost to float64 rounding once jd reaches ~2.4 million days.
func JDay(year, mon, day, hr, minute int, sec float64) (jd, jdFrac float64) {
	jd = 367.0*float64(year) -
		math.Floor((7*(float64(year)+math.Floor((float64(mon)+9)/12.0)))*0.25) +
		math.Floor(275*float64(mon)/9.0) +
		float64(day) + 1721013.5 // use - 678987.0 to go to mjd directly
	jdFrac = (sec + float64(minute)*60.0 + float64(hr)*3600.0) / 86400.0

	// check that the day and fractional day are correct
	if math.Abs(jdFrac) > 1.0 {
		dtt := math.Floor(jdFrac)
		jd += dtt
		jdFrac -= dtt
	}
	return jd, jdFrac
}

// Days2MDHMS converts a year plus a fractional day-of-year into month, day,
// hour, minute, second, replacing days2mdhms_SGP4. Uses the reference's own
// simplified "divisible by 4" leap-year rule (not full Gregorian), which is
// exact for the TLE-relevant year range (this does NOT match Gregorian at
// century boundaries like 1900/2100, but no valid two-digit TLE epoch year
// falls on one).
func Days2MDHMS(year int, days float64) (mon, day, hr, minute int, sec float64) {
	lmonth := [13]int{0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}

	dayofyr := int(math.Floor(days))
	// ----------------- find month and day of month ----------------
	if year%4 == 0 {
		lmonth[2] = 29
	}

	i := 1
	inttemp := 0
	for dayofyr > inttemp+lmonth[i] && i < 12 {
		inttemp += lmonth[i]
		i++
	}
	mon = i
	day = dayofyr - inttemp

	// ----------------- find hours minutes and seconds ---------------
	temp := (days - float64(dayofyr)) * 24.0
	hr = int(math.Floor(temp))
	temp = (temp - float64(hr)) * 60.0
	minute = int(math.Floor(temp))
	sec = (temp - float64(minute)) * 60.0
	return mon, day, hr, minute, sec
}

// GSTime returns the Greenwich sidereal time (radians, wrapped to
// [0, 2*pi)) for the given UT1 Julian date, replacing gstime_SGP4.
func GSTime(jdut1 float64) float64 {
	twopi := 2.0 * pi
	const deg2rad = pi / 180.0

	tut1 := (jdut1 - 2451545.0) / 36525.0
	temp := -6.2e-6*tut1*tut1*tut1 + 0.093104*tut1*tut1 +
		(876600.0*3600+8640184.812866)*tut1 + 67310.54841 // sec
	temp = math.Mod(temp*deg2rad/240.0, twopi) // 360/86400 = 1/240, to deg, to rad

	// ------------------------ check quadrants ---------------------
	if temp < 0.0 {
		temp += twopi
	}
	return temp
}

// InvJDay converts a split Julian date (jd, jdfrac) back into a calendar
// date/time, replacing invjday_SGP4.
func InvJDay(jd, jdfrac float64) (year, mon, day, hr, minute int, sec float64) {
	// check jdfrac for multiple days
	if math.Abs(jdfrac) >= 1.0 {
		dtt := math.Floor(jdfrac)
		jd += dtt
		jdfrac -= dtt
	}

	// check for fraction of a day included in the jd
	dt := jd - math.Floor(jd) - 0.5
	if math.Abs(dt) > 0.00000001 {
		jd -= dt
		jdfrac += dt
	}

	// --------------- find year and days of the year ---------------
	temp := jd - 2415019.5
	tu := temp / 365.25
	year = 1900 + int(math.Floor(tu))
	leapyrs := int(math.Floor(float64(year-1901) * 0.25))

	days := math.Floor(temp - (float64(year-1900)*365.0 + float64(leapyrs)))

	// ------------ check for case of beginning of a year -----------
	if days+jdfrac < 1.0 {
		year--
		leapyrs = int(math.Floor(float64(year-1901) * 0.25))
		days = math.Floor(temp - (float64(year-1900)*365.0 + float64(leapyrs)))
	}

	// ----------------- find remaining data ------------------------
	mon, day2, hr, minute, sec := Days2MDHMS(year, days+jdfrac)
	return year, mon, day2, hr, minute, sec
}
