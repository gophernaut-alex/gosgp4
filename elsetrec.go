package sgp4

// ElsetRec mirrors the C++ elsetrec struct (SGP4.h), field order and
// grouping preserved for side-by-side comparison with the original.
type ElsetRec struct {
	SatnumStr     string
	Satnum        int
	Epochyr       int
	Epochtynumrev int
	Error         int
	Operationmode byte
	Init          byte
	Method        byte

	// Near Earth
	Isimp   int
	Aycof   float64
	Con41   float64
	Cc1     float64
	Cc4     float64
	Cc5     float64
	D2      float64
	D3      float64
	D4      float64
	Delmo   float64
	Eta     float64
	Argpdot float64
	Omgcof  float64
	Sinmao  float64
	T       float64
	T2cof   float64
	T3cof   float64
	T4cof   float64
	T5cof   float64
	X1mth2  float64
	X7thm1  float64
	Mdot    float64
	Nodedot float64
	Xlcof   float64
	Xmcof   float64
	Nodecf  float64

	// Deep Space
	Irez  int
	D2201 float64
	D2211 float64
	D3210 float64
	D3222 float64
	D4410 float64
	D4422 float64
	D5220 float64
	D5232 float64
	D5421 float64
	D5433 float64
	Dedt  float64
	Del1  float64
	Del2  float64
	Del3  float64
	Didt  float64
	Dmdt  float64
	Dnodt float64
	Domdt float64
	E3    float64
	Ee2   float64
	Peo   float64
	Pgho  float64
	Pho   float64
	Pinco float64
	Plo   float64
	Se2   float64
	Se3   float64
	Sgh2  float64
	Sgh3  float64
	Sgh4  float64
	Sh2   float64
	Sh3   float64
	Si2   float64
	Si3   float64
	Sl2   float64
	Sl3   float64
	Sl4   float64
	Gsto  float64
	Xfact float64
	Xgh2  float64
	Xgh3  float64
	Xgh4  float64
	Xh2   float64
	Xh3   float64
	Xi2   float64
	Xi3   float64
	Xl2   float64
	Xl3   float64
	Xl4   float64
	Xlamo float64
	Zmol  float64
	Zmos  float64
	Atime float64
	Xli   float64
	Xni   float64

	A           float64
	Altp        float64
	Alta        float64
	Epochdays   float64
	Jdsatepoch  float64
	JdsatepochF float64
	Nddot       float64
	Ndot        float64
	Bstar       float64
	Rcse        float64
	Inclo       float64
	Nodeo       float64
	Ecco        float64
	Argpo       float64
	Mo          float64
	NoKozai     float64

	// sgp4fix add new variables from tle
	Classification byte
	Intldesg       string
	Ephtype        int
	Elnum          int32
	Revnum         int32

	// sgp4fix add unkozai'd variable
	NoUnkozai float64

	// sgp4fix add singly averaged variables
	Am    float64
	Em    float64
	Im    float64
	NodeM float64 // C++ "Om" (RAAN, singly averaged)
	ArgpM float64 // C++ "om" (argument of perigee, singly averaged)
	Mm    float64
	Nm    float64

	// sgp4fix add constant parameters to eliminate multiple calls during execution
	Tumin         float64
	Mus           float64
	Radiusearthkm float64
	Xke           float64
	J2            float64
	J3            float64
	J4            float64
	J3oj2         float64

	// Additional elements to capture relevant TLE and object information
	DiaMm      int32   // RSO dia in mm
	PeriodSec  float64 // Period in seconds
	Active     bool    // "Active S/C" flag
	NotOrbital bool    // "Orbiting S/C" flag
	RcsM2      float64 // "RCS (m^2)" storage
}
