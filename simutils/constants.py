
# from https://refractiveindex.info/
# At 1.55um
n_SiO2 = 1.4440
n_Si3N4 = 1.9963

# approximations
# from https://www.computational-photonics.eu/oms.html
# TE0 mode -> only Ex=0, Ey!=0 field (x being along thickness)
# chosen thick and to be in the CMi standards
h_wg = [200e-3, 500e-3]
n_eff = [1.584222801, 1.79523596]

# to stay in the single mode
w_wg_max = [1.1893, 0.7265]
