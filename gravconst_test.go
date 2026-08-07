package sgp4

import "testing"

// Reference values transcribed directly from getgravconst in SGP4.cpp
// (lines 2097-2126). Every other test in this repo hardcodes WGS72, so
// this is the only place WGS72Old/WGS84 are exercised at all -- without
// it, a transcription error in either literal-constant branch (as WGS72Old's
// Tumin once had: 1.0/60.0 instead of 1.0/xke) is invisible to `go test`.
func TestGetGravConst(t *testing.T) {
	cases := []struct {
		name  string
		model GravConstType
		want  GravConst
	}{
		{
			name:  "WGS72Old",
			model: WGS72Old,
			want: GravConst{
				Mus:           398600.79964,
				Radiusearthkm: 6378.135,
				Xke:           0.0743669161,
				Tumin:         1.0 / 0.0743669161,
				J2:            0.001082616,
				J3:            -0.00000253881,
				J4:            -0.00000165597,
				J3oj2:         -0.00000253881 / 0.001082616,
			},
		},
		{
			name:  "WGS72",
			model: WGS72,
			want: GravConst{
				Mus:           398600.8,
				Radiusearthkm: 6378.135,
				Xke:           0.07436691613317342, // 60/sqrt(6378.135^3/398600.8)
				Tumin:         13.446839696959309,  // 1/Xke
				J2:            0.001082616,
				J3:            -0.00000253881,
				J4:            -0.00000165597,
				J3oj2:         -0.00000253881 / 0.001082616,
			},
		},
		{
			name:  "WGS84",
			model: WGS84,
			want: GravConst{
				Mus:           398600.5,
				Radiusearthkm: 6378.137,
				Xke:           0.07436685316871385, // 60/sqrt(6378.137^3/398600.5)
				Tumin:         13.446851082044981,  // 1/Xke
				J2:            0.00108262998905,
				J3:            -0.00000253215306,
				J4:            -0.00000161098761,
				J3oj2:         -0.00000253215306 / 0.00108262998905,
			},
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got, err := GetGravConst(tc.model)
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			almostEqual(t, "Mus", got.Mus, tc.want.Mus)
			almostEqual(t, "Radiusearthkm", got.Radiusearthkm, tc.want.Radiusearthkm)
			almostEqual(t, "Xke", got.Xke, tc.want.Xke)
			almostEqual(t, "Tumin", got.Tumin, tc.want.Tumin)
			almostEqual(t, "J2", got.J2, tc.want.J2)
			almostEqual(t, "J3", got.J3, tc.want.J3)
			almostEqual(t, "J4", got.J4, tc.want.J4)
			almostEqual(t, "J3oj2", got.J3oj2, tc.want.J3oj2)
		})
	}

	t.Run("unknown model", func(t *testing.T) {
		_, err := GetGravConst(GravConstType(99))
		if err != ErrUnknownGravConst {
			t.Fatalf("got %v, want ErrUnknownGravConst", err)
		}
	})
}
