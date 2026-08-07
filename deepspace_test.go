package sgp4

import "testing"

// See initl_test.go for how these reference vectors were generated (same
// three satellites: 23333/irez0, 26900/irez1, 08195/irez2).

type deepspaceFixture struct {
	name                            string
	epoch, ecco, argpo, inclo       float64
	nodeo, mo, noKozai, noUnkozai   float64
	xke, j2, mdot, nodedot, argpdot float64

	wantDscom DscomResult

	wantDpperYep, wantDpperYinclp, wantDpperYnodep, wantDpperYargpp, wantDpperYmp float64

	wantDpperNep, wantDpperNinclp, wantDpperNnodep, wantDpperNargpp, wantDpperNmp float64

	wantIrez                                                             int
	wantD2201, wantD2211, wantD3210, wantD3222, wantD4410, wantD4422     float64
	wantD5220, wantD5232, wantD5421, wantD5433                           float64
	wantDedt, wantDidt, wantDmdt, wantDnodt, wantDomdt                   float64
	wantDel1, wantDel2, wantDel3, wantXfact, wantXlamo, wantXli, wantXni float64

	wantDspaceEm, wantDspaceArgpm, wantDspaceInclm, wantDspaceMm float64
	wantDspaceNodem, wantDspaceNm, wantDspaceDndt                float64
	wantDspaceAtime, wantDspaceXli, wantDspaceXni                float64
}

var deepspaceFixtures = []deepspaceFixture{
	{
		name: "irez0_23333", epoch: 16376.499999990221,
		ecco: 0.97282979999999997, argpo: 0.53120841113699413, inclo: 0.50176470665584982,
		nodeo: 0.041399209857305497, mo: 0.02356194490192345,
		noKozai: 0.00031893671148723211, noUnkozai: 0.00031891772263956645,
		xke: 0.074366916133173422, j2: 0.001082616,
		mdot: 0.00031893671733400098, nodedot: -1.10243921949836e-07, argpdot: 1.7876183534731737e-07,

		wantDscom: DscomResult{
			Snodm: 0.04138738522378331, Cnodm: 0.99914317509771244, Sinim: 0.48097346707760225, Cosim: 0.87673515041165695,
			Sinomm: 0.50657559751980841, Cosomm: 0.862195548583643, Day: 34637.999999990221,
			E3: 0.0066181304675952462, Ee2: -0.0075031988777229568, Em: 0.97282979999999997, Emsq: 0.94639781976803994,
			Gam: 73.184170183980982, Rtemsq: 0.23152144659180079, Nm: 0.00031891772263956645,
			Se2: 0.046298098090509623, Se3: 0.042716268514968528,
			Sgh2: -0.043914801636940455, Sgh3: 0.047764725328650261, Sgh4: -0.00065367135928805328,
			Sh2: 0.050870556699581897, Sh3: 0.021941234461940802,
			Si2: 0.087500685337536258, Si3: 0.035116477158559256,
			Sl2: 0.36922521742036968, Sl3: -0.4025541103841882, Sl4: 0.0092599037908010079,
			Xgh2: -0.0068782929297384471, Xgh3: -0.0077755164029081215, Xgh4: -0.0003441197687525333,
			Xh2: -0.011795764942147404, Xh3: 0.011690374638795664,
			Xi2: -0.019565898713799566, Xi3: 0.020395664959604404,
			Xl2: 0.058253284233195868, Xl3: 0.065726324968919767, Xl4: 0.004874798177836254,
			Zmol: 5.6078681962094947, Zmos: 5.1960588437711053,
			S1: -0.0050815113207585494, S2: -0.0032482715059784094, S3: 0.0015040890359740973,
			S4: 0.0003482288694115901, S5: -0.0023285288268560956, S6: 0.73828418398592754, S7: -0.65119706026821522,
			Ss1: -0.031637362075717664, Ss2: -0.020223656953210326, Ss3: 0.00936442062636717,
			Ss4: 0.0021680642099106244, Ss5: -0.002895833784881056, Ss6: -0.73169972230466684, Ss7: -0.67509213335700569,
			Sz1: 3.995477449554401, Sz2: -19.714258476426789, Sz3: 25.489285045655244,
			Sz11: 0.83102261134876876, Sz12: -2.1633249995284927, Sz13: -0.037180336313098687,
			Sz21: -0.52083791822072945, Sz22: 1.2576992582814417, Sz23: 0.021626644758762421,
			Sz31: -1.0025138942857206, Sz32: -10.127652455171239, Sz33: 10.013009795074776,
			Z1: 25.678336187276685, Z2: -19.364972032878736, Z3: 3.8291226776677783,
			Z11: 1.7176024130716347, Z12: 3.0117400404782568, Z13: -1.4218619023444872,
			Z21: -1.0187923322117183, Z22: -1.8156987370725359, Z23: 0.7806838841820859,
			Z31: 10.127190382863452, Z32: -9.8761095559952405, Z33: -1.0371861032770866,
		},

		wantDpperYep: 0.97282979999999997, wantDpperYinclp: 0.50176470665584982,
		wantDpperYnodep: 0.041399209857305497, wantDpperYargpp: 0.53120841113699413, wantDpperYmp: 0.02356194490192345,

		wantDpperNep: 0.99027630581636905, wantDpperNinclp: 0.5278669890199682,
		wantDpperNnodep: 0.072075752983924435, wantDpperNargpp: 0.50659223947279997, wantDpperNmp: 0.0029266364731296741,

		wantIrez: 0,
		wantDedt: 2.9681356032420809e-09, wantDidt: -3.4390448847643113e-07, wantDmdt: -3.4380511056252244e-06,
		wantDnodt: -5.0539252314737997e-07, wantDomdt: 6.914571527622915e-07,

		wantDspaceEm: 0.40000296813560327, wantDspaceArgpm: 0.3006914571527623, wantDspaceInclm: 0.50142080216737339,
		wantDspaceMm: 0.096561948894374777, wantDspaceNodem: 0.19949460747685263, wantDspaceNm: 0.00031891772263956645,
		wantDspaceDndt: 0, wantDspaceAtime: 0, wantDspaceXli: 0, wantDspaceXni: 0.00031891772263956645,
	},
	{
		name: "irez1_26900", epoch: 20560.745032470208,
		ecco: 0.00033189999999999999, argpo: 1.5041142773932012, inclo: 0.00028623399732707005,
		nodeo: 4.6519621910221423, mo: 3.1810196413923451,
		noKozai: 0.0043752719594775955, noUnkozai: 0.0043751093799833923,
		xke: 0.074366916133173422, j2: 0.001082616,
		mdot: 0.0043752719564557393, nodedot: -1.6259675948021914e-07, argpdot: 3.2519953954951574e-07,

		wantDscom: DscomResult{
			Snodm: -0.99817485702518638, Cnodm: -0.060390022377448352, Sinim: 0.00028623399341854957, Cosim: 0.99999995903504968,
			Sinomm: 0.99777757582687476, Cosomm: 0.066632643479342502, Day: 38822.245032470208,
			E3: -4.2973571095838277e-07, Ee2: -8.5758630143127902e-07, Em: 0.00033189999999999999, Emsq: 1.1015760999999999e-07,
			Gam: 81.319882329294032, Rtemsq: 0.99999994492119348, Nm: 0.0043751093799833923,
			Se2: 4.9777105032906995e-06, Se3: -3.75812288083242e-06,
			Sgh2: 0.011702388132002264, Sgh3: 0.015966174072024737, Sgh4: -0.00020580595860649297,
			Sh2: -0.0014694899495776008, Sh3: 0.00070314570647507925,
			Si2: 0.0006596584265507404, Si3: 0.0013421801429428528,
			Sl2: -0.012195666993036194, Sl3: -0.017165372384922488, Sl4: 0.00048021395253585413,
			Xgh2: 0.0014686108663714017, Xgh3: -0.0027627719033653901, Xgh4: -0.00010834480947841367,
			Xh2: 0.00024439772292797517, Xh3: 0.00019800681546381559,
			Xi2: 0.00017005953655943311, Xi3: -0.00021812405267991613,
			Xl2: -0.0016726101706792886, Xl3: 0.0029843084633158676, Xl4: 0.00025280458130881177,
			Zmol: 4.6850958465720396, Zmos: 1.7751219687009936,
			S1: -5.4583549280336207e-07, S2: -5.4819278189344382e-05, S3: 0.00010963855033992794,
			S4: 0.00010963854430116743, S5: 0.031348946797160819, S6: 0.78557213000824921, S7: 0.39364947555104812,
			Ss1: -3.3983580926116924e-06, Ss2: -0.00034130345153828853, Ss3: 0.00068260686547940352,
			Ss4: 0.00068260682788223205, Ss5: 0.020006300997783919, Ss6: -0.73236992212690122, Ss7: 0.55293214817515635,
			Sz1: 3.9565757890765361, Sz2: 8.9331558249645067, Sz3: 16.529971438171735,
			Sz11: 2.0761535540749438, Sz12: -0.96638112444745938, Sz13: 0.10989722573499631,
			Sz21: -0.5818957556801565, Sz22: -2.1527616303826753, Sz23: 0.44819301621173591,
			Sz31: -1.1293252850056297, Sz32: 8.5718364173916815, Sz33: 10.565674984446208,
			Z1: 16.928702780372589, Z2: 7.6278378612881061, Z3: 3.3189439224563126,
			Z11: 0.26236145042624565, Z12: -1.5510924457273207, Z13: 2.2518445290479785,
			Z21: -0.99372477883819654, Z22: 2.2291220442909858, Z23: 0.8122714145073735,
			Z31: 11.110248359788896, Z32: 6.6975116996138553, Z33: -1.4892070642830382,
		},

		wantDpperYep: 0.00033189999999999999, wantDpperYinclp: 0.00028623399732707005,
		wantDpperYnodep: 4.6519621910221423, wantDpperYargpp: 1.5041142773932012, wantDpperYmp: 3.1810196413923451,

		wantDpperNep: 0.00033238418392316654, wantDpperNinclp: 0.00064022510573909263,
		wantDpperNnodep: 4.4576848941677341, wantDpperNargpp: 1.7032378821738905, wantDpperNmp: 3.1759860787505869,

		wantIrez: 1,
		wantDedt: -3.5218071145820684e-12, wantDidt: -3.073811485140335e-08, wantDmdt: -1.6136212111376092e-07,
		wantDnodt: 0, wantDomdt: 9.0887925667102694e-08,
		wantDel1: -6.3981362169008843e-13, wantDel2: 1.4105229361725126e-11, wantDel3: 1.978672345715981e-12,
		wantXfact: -0.0043750143829543297, wantXlamo: 1.0826342773273936, wantXli: 1.0826342773273936, wantXni: 0.0043751093799833923,

		wantDspaceEm: 0.3999999964781929, wantDspaceArgpm: 0.30009088792566707, wantDspaceInclm: 0.00025549588247566671,
		wantDspaceMm: 0.64600241558380789, wantDspaceNodem: 0.20000000000000001, wantDspaceNm: 0.0043751168251923902,
		wantDspaceDndt: 7.4452089979248903e-09, wantDspaceAtime: 720, wantDspaceXli: 1.0827046053460336, wantDspaceXni: 0.0043751147408247213,
	},
	{
		name: "irez2_08195", epoch: 20630.332154440228,
		ecco: 0.68771459999999995, argpo: 4.6210227393720391, inclo: 1.1197788134700339,
		nodeo: 4.8707200141378602, mo: 0.35300505852061709,
		noKozai: 0.0087480868880674655, noUnkozai: 0.008748547019630239,
		xke: 0.074366916133173422, j2: 0.001082616,
		mdot: 0.0087480868866331284, nodedot: -1.2845672158012262e-06, argpdot: -7.4370092653335332e-08,

		wantDscom: DscomResult{
			Snodm: -0.98749180502655509, Cnodm: 0.1576703364789844, Sinim: 0.90000405307080833, Cosim: 0.43588152571096384,
			Sinomm: -0.99582900776036032, Cosomm: -0.091239176360904298, Day: 38891.832154440228,
			E3: -0.00044080461331893458, Ee2: -0.0001093428053937556, Em: 0.68771459999999995, Emsq: 0.47295137105315993,
			Gam: 81.455185302464628, Rtemsq: 0.7259811491676903, Nm: 0.008748547019630239,
			Se2: 0.0025019415226424269, Se3: 0.0010072419885301313,
			Sgh2: -0.0012619155097006729, Sgh3: 0.0024557462815410938, Sgh4: -7.4719902613734867e-05,
			Sh2: 0.00094580372708250797, Sh3: -0.0039702236906196025,
			Si2: -0.0011105699174389215, Si3: 2.7357459043652346e-05,
			Sl2: 0.0040834926664333263, Sl3: -0.0021191158730878212, Sl4: 0.00028883025414443435,
			Xgh2: 0.00045492024438304162, Xgh3: -2.4648151078821564e-05, Xgh4: -3.9335661939747741e-05,
			Xh2: -0.000643834081488151, Xh3: 0.00030794292793486264,
			Xi2: 0.00011424082985422945, Xi3: 0.00012000896281451965,
			Xl2: -0.00078626437719944763, Xl3: -0.00045594235604769131, Xl4: 0.00015205224896677633,
			Zmol: 1.703291771990898, Zmos: 2.9721580403255459,
			S1: -0.0004106209626082095, S2: -3.7762516647879401e-05, S3: 5.482975046298305e-05,
			S4: 3.9805365249694133e-05, S5: 0.31998284015161677, S6: 0.13314323348134094, S7: 0.53675366513073586,
			Ss1: -0.0025565158178131155, Ss2: -0.00023510848184058264, Ss3: 0.00034136865165139445,
			Ss4: 0.00024782720601570433, Ss5: 0.26454981282358392, Ss6: -0.48932643115477137, Ss7: -0.19699506287266827,
			Sz1: 2.8238981585147602, Sz2: -5.9810598405552886, Sz3: 5.9277506400423574,
			Sz11: 0.55800941905082369, Sz12: 2.3618244410934377, Sz13: 0.49982891694923998,
			Sz21: 0.96685403313051566, Sz22: 2.0114198341083638, Sz23: -7.476532737814753,
			Sz31: -1.9673875696895442, Sz32: -2.5459583917124653, Sz33: 2.9871658891939767,
			Z1: 2.9000521897765901, Z2: 7.1700524857419747, Z3: 7.0578529475559773,
			Z11: 1.3641748695339637, Z12: -1.5126220389318885, Z13: -0.22482095843952615,
			Z21: -5.1187451627829477, Z22: -8.5247771949583004, Z23: -1.0413828036172537,
			Z31: 1.0839759284760986, Z32: 5.7143081281804999, Z33: 0.77436752613879101,
		},

		wantDpperYep: 0.68771459999999995, wantDpperYinclp: 1.1197788134700339,
		wantDpperYnodep: 4.8707200141378602, wantDpperYargpp: 4.6210227393720391, wantDpperYmp: 0.35300505852061709,

		wantDpperNep: 0.68712086232186131, wantDpperNinclp: 1.120084249358912,
		wantDpperNnodep: 4.8700123851223474, wantDpperNargpp: 4.6218691626583723, wantDpperNmp: 0.35183452479212357,

		wantIrez:  2,
		wantD2201: -1.1973595516231104e-11, wantD2211: 6.4532138341213982e-11,
		wantD3210: -3.8937227381311457e-12, wantD3222: -7.3648575380239341e-12,
		wantD4410: 2.5769601409463596e-12, wantD4422: 4.3614555927145512e-12,
		wantD5220: -2.5287894659528222e-12, wantD5232: 6.767712568551256e-13,
		wantD5421: -2.2806980465619702e-12, wantD5433: -1.659570821491405e-12,
		wantDedt: -2.8885476234610706e-08, wantDidt: -9.7841082692141341e-09, wantDmdt: 9.2707867425269653e-08,
		wantDnodt: -6.1243234896952638e-08, wantDomdt: -1.4155210733836969e-08,
		wantXfact: -0.0087535972220536828, wantXlamo: 2.6628995258087294, wantXli: 2.6628995258087294, wantXni: 0.008748547019630239,

		wantDspaceEm: 0.39997111452376544, wantDspaceArgpm: 0.29998584478926615, wantDspaceInclm: 1.1197690293617648,
		wantDspaceMm: 12.156870540492807, wantDspaceNodem: 0.19993875676510306, wantDspaceNm: 0.0087485477419244392,
		wantDspaceDndt: 7.2229420014502388e-10, wantDspaceAtime: 720, wantDspaceXli: 2.659263514896685, wantDspaceXni: 0.0087485474986700885,
	},
}

func checkDscomResult(t *testing.T, got, want DscomResult) {
	t.Helper()
	almostEqual(t, "Snodm", got.Snodm, want.Snodm)
	almostEqual(t, "Cnodm", got.Cnodm, want.Cnodm)
	almostEqual(t, "Sinim", got.Sinim, want.Sinim)
	almostEqual(t, "Cosim", got.Cosim, want.Cosim)
	almostEqual(t, "Sinomm", got.Sinomm, want.Sinomm)
	almostEqual(t, "Cosomm", got.Cosomm, want.Cosomm)
	almostEqual(t, "Day", got.Day, want.Day)
	almostEqual(t, "E3", got.E3, want.E3)
	almostEqual(t, "Ee2", got.Ee2, want.Ee2)
	almostEqual(t, "Em", got.Em, want.Em)
	almostEqual(t, "Emsq", got.Emsq, want.Emsq)
	almostEqual(t, "Gam", got.Gam, want.Gam)
	almostEqual(t, "Peo", got.Peo, want.Peo)
	almostEqual(t, "Pgho", got.Pgho, want.Pgho)
	almostEqual(t, "Pho", got.Pho, want.Pho)
	almostEqual(t, "Pinco", got.Pinco, want.Pinco)
	almostEqual(t, "Plo", got.Plo, want.Plo)
	almostEqual(t, "Rtemsq", got.Rtemsq, want.Rtemsq)
	almostEqual(t, "Se2", got.Se2, want.Se2)
	almostEqual(t, "Se3", got.Se3, want.Se3)
	almostEqual(t, "Sgh2", got.Sgh2, want.Sgh2)
	almostEqual(t, "Sgh3", got.Sgh3, want.Sgh3)
	almostEqual(t, "Sgh4", got.Sgh4, want.Sgh4)
	almostEqual(t, "Sh2", got.Sh2, want.Sh2)
	almostEqual(t, "Sh3", got.Sh3, want.Sh3)
	almostEqual(t, "Si2", got.Si2, want.Si2)
	almostEqual(t, "Si3", got.Si3, want.Si3)
	almostEqual(t, "Sl2", got.Sl2, want.Sl2)
	almostEqual(t, "Sl3", got.Sl3, want.Sl3)
	almostEqual(t, "Sl4", got.Sl4, want.Sl4)
	almostEqual(t, "S1", got.S1, want.S1)
	almostEqual(t, "S2", got.S2, want.S2)
	almostEqual(t, "S3", got.S3, want.S3)
	almostEqual(t, "S4", got.S4, want.S4)
	almostEqual(t, "S5", got.S5, want.S5)
	almostEqual(t, "S6", got.S6, want.S6)
	almostEqual(t, "S7", got.S7, want.S7)
	almostEqual(t, "Ss1", got.Ss1, want.Ss1)
	almostEqual(t, "Ss2", got.Ss2, want.Ss2)
	almostEqual(t, "Ss3", got.Ss3, want.Ss3)
	almostEqual(t, "Ss4", got.Ss4, want.Ss4)
	almostEqual(t, "Ss5", got.Ss5, want.Ss5)
	almostEqual(t, "Ss6", got.Ss6, want.Ss6)
	almostEqual(t, "Ss7", got.Ss7, want.Ss7)
	almostEqual(t, "Sz1", got.Sz1, want.Sz1)
	almostEqual(t, "Sz2", got.Sz2, want.Sz2)
	almostEqual(t, "Sz3", got.Sz3, want.Sz3)
	almostEqual(t, "Sz11", got.Sz11, want.Sz11)
	almostEqual(t, "Sz12", got.Sz12, want.Sz12)
	almostEqual(t, "Sz13", got.Sz13, want.Sz13)
	almostEqual(t, "Sz21", got.Sz21, want.Sz21)
	almostEqual(t, "Sz22", got.Sz22, want.Sz22)
	almostEqual(t, "Sz23", got.Sz23, want.Sz23)
	almostEqual(t, "Sz31", got.Sz31, want.Sz31)
	almostEqual(t, "Sz32", got.Sz32, want.Sz32)
	almostEqual(t, "Sz33", got.Sz33, want.Sz33)
	almostEqual(t, "Xgh2", got.Xgh2, want.Xgh2)
	almostEqual(t, "Xgh3", got.Xgh3, want.Xgh3)
	almostEqual(t, "Xgh4", got.Xgh4, want.Xgh4)
	almostEqual(t, "Xh2", got.Xh2, want.Xh2)
	almostEqual(t, "Xh3", got.Xh3, want.Xh3)
	almostEqual(t, "Xi2", got.Xi2, want.Xi2)
	almostEqual(t, "Xi3", got.Xi3, want.Xi3)
	almostEqual(t, "Xl2", got.Xl2, want.Xl2)
	almostEqual(t, "Xl3", got.Xl3, want.Xl3)
	almostEqual(t, "Xl4", got.Xl4, want.Xl4)
	almostEqual(t, "Nm", got.Nm, want.Nm)
	almostEqual(t, "Z1", got.Z1, want.Z1)
	almostEqual(t, "Z2", got.Z2, want.Z2)
	almostEqual(t, "Z3", got.Z3, want.Z3)
	almostEqual(t, "Z11", got.Z11, want.Z11)
	almostEqual(t, "Z12", got.Z12, want.Z12)
	almostEqual(t, "Z13", got.Z13, want.Z13)
	almostEqual(t, "Z21", got.Z21, want.Z21)
	almostEqual(t, "Z22", got.Z22, want.Z22)
	almostEqual(t, "Z23", got.Z23, want.Z23)
	almostEqual(t, "Z31", got.Z31, want.Z31)
	almostEqual(t, "Z32", got.Z32, want.Z32)
	almostEqual(t, "Z33", got.Z33, want.Z33)
	almostEqual(t, "Zmol", got.Zmol, want.Zmol)
	almostEqual(t, "Zmos", got.Zmos, want.Zmos)
}

func TestDscom(t *testing.T) {
	for _, f := range deepspaceFixtures {
		t.Run(f.name, func(t *testing.T) {
			got := dscom(f.epoch, f.ecco, f.argpo, 0.0, f.inclo, f.nodeo, f.noUnkozai)
			checkDscomResult(t, got, f.wantDscom)
		})
	}
}

func TestDpper(t *testing.T) {
	for _, f := range deepspaceFixtures {
		t.Run(f.name+"/init=y", func(t *testing.T) {
			dsc := dscom(f.epoch, f.ecco, f.argpo, 0.0, f.inclo, f.nodeo, f.noUnkozai)
			ep, inclp, nodep, argpp, mp := dpper(
				dsc.E3, dsc.Ee2, dsc.Peo, dsc.Pgho, dsc.Pho, dsc.Pinco, dsc.Plo, dsc.Se2, dsc.Se3,
				dsc.Sgh2, dsc.Sgh3, dsc.Sgh4, dsc.Sh2, dsc.Sh3, dsc.Si2, dsc.Si3, dsc.Sl2, dsc.Sl3, dsc.Sl4,
				0.0, dsc.Xgh2, dsc.Xgh3, dsc.Xgh4, dsc.Xh2, dsc.Xh3, dsc.Xi2, dsc.Xi3, dsc.Xl2, dsc.Xl3, dsc.Xl4,
				dsc.Zmol, dsc.Zmos, f.inclo,
				'y',
				f.ecco, f.inclo, f.nodeo, f.argpo, f.mo,
				'a',
			)
			almostEqual(t, "ep", ep, f.wantDpperYep)
			almostEqual(t, "inclp", inclp, f.wantDpperYinclp)
			almostEqual(t, "nodep", nodep, f.wantDpperYnodep)
			almostEqual(t, "argpp", argpp, f.wantDpperYargpp)
			almostEqual(t, "mp", mp, f.wantDpperYmp)
		})

		t.Run(f.name+"/init=n,t=200", func(t *testing.T) {
			dsc := dscom(f.epoch, f.ecco, f.argpo, 0.0, f.inclo, f.nodeo, f.noUnkozai)
			ep, inclp, nodep, argpp, mp := dpper(
				dsc.E3, dsc.Ee2, dsc.Peo, dsc.Pgho, dsc.Pho, dsc.Pinco, dsc.Plo, dsc.Se2, dsc.Se3,
				dsc.Sgh2, dsc.Sgh3, dsc.Sgh4, dsc.Sh2, dsc.Sh3, dsc.Si2, dsc.Si3, dsc.Sl2, dsc.Sl3, dsc.Sl4,
				200.0, dsc.Xgh2, dsc.Xgh3, dsc.Xgh4, dsc.Xh2, dsc.Xh3, dsc.Xi2, dsc.Xi3, dsc.Xl2, dsc.Xl3, dsc.Xl4,
				dsc.Zmol, dsc.Zmos, f.inclo,
				'n',
				dsc.Em, f.inclo, f.nodeo, f.argpo, f.mo,
				'a',
			)
			almostEqual(t, "ep", ep, f.wantDpperNep)
			almostEqual(t, "inclp", inclp, f.wantDpperNinclp)
			almostEqual(t, "nodep", nodep, f.wantDpperNnodep)
			almostEqual(t, "argpp", argpp, f.wantDpperNargpp)
			almostEqual(t, "mp", mp, f.wantDpperNmp)
		})
	}
}

func TestDsinit(t *testing.T) {
	for _, f := range deepspaceFixtures {
		t.Run(f.name, func(t *testing.T) {
			dsc := dscom(f.epoch, f.ecco, f.argpo, 0.0, f.inclo, f.nodeo, f.noUnkozai)
			xpidot := f.argpdot + f.nodedot
			initRes := initl(f.xke, f.j2, f.ecco, f.epoch, f.inclo, f.noKozai, 'a')
			s := &ElsetRec{}
			_, _, _, _, _, _, _ = s.dsinit(
				f.xke, dsc.Cosim, dsc.Emsq, f.argpo, dsc.S1, dsc.S2, dsc.S3, dsc.S4, dsc.S5, dsc.Sinim,
				dsc.Ss1, dsc.Ss2, dsc.Ss3, dsc.Ss4, dsc.Ss5, dsc.Sz1, dsc.Sz3, dsc.Sz11, dsc.Sz13, dsc.Sz21, dsc.Sz23,
				dsc.Sz31, dsc.Sz33, 0.0, 0.0, initRes.Gsto, f.mo, f.mdot, f.noUnkozai, f.nodeo, f.nodedot,
				xpidot, dsc.Z1, dsc.Z3, dsc.Z11, dsc.Z13, dsc.Z21, dsc.Z23, dsc.Z31, dsc.Z33, f.ecco,
				f.ecco*f.ecco, dsc.Em, 0.0, f.inclo, 0.0, dsc.Nm, 0.0,
			)
			if s.Irez != f.wantIrez {
				t.Fatalf("Irez: got %d, want %d", s.Irez, f.wantIrez)
			}
			almostEqual(t, "D2201", s.D2201, f.wantD2201)
			almostEqual(t, "D2211", s.D2211, f.wantD2211)
			almostEqual(t, "D3210", s.D3210, f.wantD3210)
			almostEqual(t, "D3222", s.D3222, f.wantD3222)
			almostEqual(t, "D4410", s.D4410, f.wantD4410)
			almostEqual(t, "D4422", s.D4422, f.wantD4422)
			almostEqual(t, "D5220", s.D5220, f.wantD5220)
			almostEqual(t, "D5232", s.D5232, f.wantD5232)
			almostEqual(t, "D5421", s.D5421, f.wantD5421)
			almostEqual(t, "D5433", s.D5433, f.wantD5433)
			almostEqual(t, "Dedt", s.Dedt, f.wantDedt)
			almostEqual(t, "Didt", s.Didt, f.wantDidt)
			almostEqual(t, "Dmdt", s.Dmdt, f.wantDmdt)
			almostEqual(t, "Dnodt", s.Dnodt, f.wantDnodt)
			almostEqual(t, "Domdt", s.Domdt, f.wantDomdt)
			almostEqual(t, "Del1", s.Del1, f.wantDel1)
			almostEqual(t, "Del2", s.Del2, f.wantDel2)
			almostEqual(t, "Del3", s.Del3, f.wantDel3)
			almostEqual(t, "Xfact", s.Xfact, f.wantXfact)
			almostEqual(t, "Xlamo", s.Xlamo, f.wantXlamo)
			almostEqual(t, "Xli", s.Xli, f.wantXli)
			almostEqual(t, "Xni", s.Xni, f.wantXni)
		})
	}
}

func TestDspace(t *testing.T) {
	for _, f := range deepspaceFixtures {
		t.Run(f.name, func(t *testing.T) {
			dsc := dscom(f.epoch, f.ecco, f.argpo, 0.0, f.inclo, f.nodeo, f.noUnkozai)
			xpidot := f.argpdot + f.nodedot
			initRes := initl(f.xke, f.j2, f.ecco, f.epoch, f.inclo, f.noKozai, 'a')
			s := &ElsetRec{}
			_, _, _, _, _, _, _ = s.dsinit(
				f.xke, dsc.Cosim, dsc.Emsq, f.argpo, dsc.S1, dsc.S2, dsc.S3, dsc.S4, dsc.S5, dsc.Sinim,
				dsc.Ss1, dsc.Ss2, dsc.Ss3, dsc.Ss4, dsc.Ss5, dsc.Sz1, dsc.Sz3, dsc.Sz11, dsc.Sz13, dsc.Sz21, dsc.Sz23,
				dsc.Sz31, dsc.Sz33, 0.0, 0.0, initRes.Gsto, f.mo, f.mdot, f.noUnkozai, f.nodeo, f.nodedot,
				xpidot, dsc.Z1, dsc.Z3, dsc.Z11, dsc.Z13, dsc.Z21, dsc.Z23, dsc.Z31, dsc.Z33, f.ecco,
				f.ecco*f.ecco, dsc.Em, 0.0, f.inclo, 0.0, dsc.Nm, 0.0,
			)
			// dspace additionally needs Argpo/Argpdot/Gsto/T/NoUnkozai set
			// on the receiver (it reads them directly -- see dspace's doc
			// comment). NoUnkozai in particular matters here even though
			// it looks unrelated: dspace's own "atime==0" warm-start reset
			// sets Xni = s.NoUnkozai, so leaving it at the Go zero value
			// on this from-scratch test *ElsetRec (dsinit doesn't set it --
			// only initl/NewElsetRec do) would silently zero out Xni.
			s.Argpo = f.argpo
			s.Argpdot = f.argpdot
			s.Gsto = initRes.Gsto
			s.NoUnkozai = f.noUnkozai
			s.T = 1000.0
			// Fixed, arbitrary (not pipeline-derived) integrator warm-start
			// state, matching the reference vectors' generation exactly:
			// a fresh atime=0 with xli=Xlamo (whatever dsinit just set it
			// to -- 0 for irez0, the resonance xlamo otherwise) and
			// xni=NoUnkozai.
			s.Atime = 0.0
			s.Xli = s.Xlamo
			s.Xni = f.noUnkozai

			em, argpm, inclm, mm, nodem, nm, dndt := s.dspace(0.4, 0.3, f.inclo, 0.1, 0.2, f.noUnkozai)
			almostEqual(t, "em", em, f.wantDspaceEm)
			almostEqual(t, "argpm", argpm, f.wantDspaceArgpm)
			almostEqual(t, "inclm", inclm, f.wantDspaceInclm)
			almostEqual(t, "mm", mm, f.wantDspaceMm)
			almostEqual(t, "nodem", nodem, f.wantDspaceNodem)
			almostEqual(t, "nm", nm, f.wantDspaceNm)
			almostEqual(t, "dndt", dndt, f.wantDspaceDndt)
			almostEqual(t, "Atime", s.Atime, f.wantDspaceAtime)
			almostEqual(t, "Xli", s.Xli, f.wantDspaceXli)
			almostEqual(t, "Xni", s.Xni, f.wantDspaceXni)
		})
	}
}
