def pdg_id(particle):

    pdg_id_dict = {
	## Quarks
	"d": 1,
	"u": 2,
	"s": 3,
	"c": 4,
	"b": 5,
	"t": 6,

	## Leptons
	"e-"    : 11,
	"e+"    : -11,
	"mu-"   : 13,
	"mu+"   : -13,
	"tau-"  : 15,
	"tau+"  : -15,
	"nu_e"  : 12,
	"nu_mu" : 14,
	"nu_tau": 16,

	## Gauge and Higgs bosons
	"g"     : 21,
	"photon": 22,
	"Z"     : 23,
	"W+"    : 24,
	"W-"    : -24,
	"H"     : 25,

	## Light I=1 mesons
	"pi0" : 111,
	"pi+" : 211,
	"pi-" : -211,
	"rho0": 113,
	"rho+": 213,
	"rho-": -213,

	## Strange mesons
	"KL0": 130,
    }

    d = pdg_id_dict   # shorthand

    pdg_id_dict["quarks"] = [d["d"], d["u"], d["s"], d["c"], d["b"], d["t"]]

    pdg_id_dict["charged_leptons"] = [d["e-"], d["e+"], d["mu-"], d["mu+"], d["tau-"], d["tau+"]]
    pdg_id_dict["neutral_leptons"] = [d["nu_e"], d["nu_mu"], d["nu_tau"]]
    pdg_id_dict["leptons"] =  pdg_id_dict["charged_leptons"] + pdg_id_dict["neutral_leptons"]

    pdg_id_dict["light_mesons"] = [d["pi0"], d["pi+"], d["pi-"], d["rho0"], d["rho+"], d["rho-"]]
    pdg_id_dict["strange_mesons"] = [d["KL0"]]
    pdg_id_dict["mesons"] = d["light_mesons"] + d["strange_mesons"]

    pdg_id_dict["hadrons"] = d["mesons"]


    return pdg_id_dict[particle]
