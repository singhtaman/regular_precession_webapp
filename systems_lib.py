#########################################################################################################################
########################## Regularly precessing & Non-precessing systems #########################
#########################################################################################################################
# This file contains the default parameters for the three systems used in the paper with q = 0.6.
# The parameters are set for the three systems - System 1, System 2 and System 3 - are 50th percentile values of the
# posterior distributions obtained from the Precession code for a population of 10000 binaries with maximally spinning binaries
# that are isotropically oriented.
import numpy as np
error_handler = np.seterr(invalid="raise")

# Constants
solar_mass = 4.92624076 * 1e-6          # [solar_mass] = sec
giga_parsec = 1.02927125 * 1e17         # [giga_parsec] = sec
year = 31557600                         # [year] = sec

# These systems are used in the paper - Add arXiv number and/or DOI identifier

# Precessing system 2 - Edge-on -- cos_i_JN = 0

default_precession_params_sys2_RP = {
    "theta_S": np.pi / 4,
    "phi_S": 0.0,
    "theta_J": np.pi / 2,
    "phi_J": np.pi / 2,
    "mcz": 10 * solar_mass,
    "dist": 1.5 * giga_parsec,
    "eta": 0.234375,
    "t_c": 0.0,
    "phi_c": 0.0,
    "theta_tilde": 5.0,             # dimensionless
    "Omega_tilde": 2.0,             # dimensionless
    "gamma_P": 0.0,
}

# Non-precessing system 2 - Edge-on -- cos_i_JN = 0
default_precession_params_sys2_NP = {
    "theta_S": np.pi / 4,
    "phi_S": 0.0,
    "theta_J": np.pi / 2,
    "phi_J": np.pi / 2,
    "mcz": 10 * solar_mass,
    "dist": 1.5 * giga_parsec,
    "eta": 0.234375,
    "t_c": 0.0,
    "phi_c": 0.0,
    "theta_tilde": 0.0,
    "Omega_tilde": 0.0,
    "gamma_P": 0.0,
}

# Precessing system 3 - Random

default_precession_params_sys3_RP = {
    'theta_S' : np.pi/4, 
    'phi_S' : 0, 
    'theta_J' : 8*np.pi/9, 
    'phi_J' : np.pi/4, 
    "mcz": 10 * solar_mass,
    "dist": 1.5 * giga_parsec,
    "eta": 0.234375,
    "t_c": 0.0,
    "phi_c": 0.0,
    "theta_tilde": 5.0,             # dimensionless
    "Omega_tilde": 2.0,             # dimensionless
    "gamma_P": 0.0,
}

# Non-precessing system 3 - Random
default_precession_params_sys3_NP = {
    'theta_S' : np.pi/4, 
    'phi_S' : 0, 
    'theta_J' : 8*np.pi/9, 
    'phi_J' : np.pi/4, 
    "mcz": 10 * solar_mass,
    "dist": 1.5 * giga_parsec,
    "eta": 0.234375,
    "t_c": 0.0,
    "phi_c": 0.0,
    "theta_tilde": 0.0,
    "Omega_tilde": 0.0,
    "gamma_P": 0.0,
}

# Precessing system 1 - Face-on -- cos_i_JN = 1

default_precession_params_sys1_RP = {
    "theta_S": np.pi / 4,
    "phi_S": 0.0,
    "theta_J": np.pi / 4,
    "phi_J": 0,
    "mcz": 10 * solar_mass,
    "dist": 1.5 * giga_parsec,
    "eta": 0.234375,
    "t_c": 0.0,
    "phi_c": 0.0,
    "theta_tilde": 5.0,             # dimensionless
    "Omega_tilde": 2.0,             # dimensionless
    "gamma_P": 0.0,
}

# Non-precessing system 1 - Face-on -- cos_i_JN = 1
default_precession_params_sys1_NP = {
    "theta_S": np.pi / 4,
    "phi_S": 0.0,
    "theta_J": np.pi / 4,
    "phi_J": 0,
    "mcz": 10 * solar_mass,
    "dist": 1.5 * giga_parsec,
    "eta": 0.234375,
    "t_c": 0.0,
    "phi_c": 0.0,
    "theta_tilde": 0.0,
    "Omega_tilde": 0.0,
    "gamma_P": 0.0,
}



#default precessing parameters - Same as RP_params_edge_on_1 - System 2 - Edge-on (cos_i_JN = 0)

rp_params = {
    "theta_S": np.pi / 4,
    "phi_S": 0.0,
    "theta_J": np.pi / 2,
    "phi_J": np.pi / 2,
    "mcz": 20 * solar_mass,
    "dist": 1.5 * giga_parsec,
    "eta": 0.234375,
    "t_c": 0.0,
    "phi_c": 0.0,
    "theta_tilde": 4.0,             # dimensionless
    "Omega_tilde": 2.0,             # dimensionless
    "gamma_P": 0.0,
}}

#default non-precessing parameters - Same as NP_params_edge_on_1 - System 2 - Edge-on (cos_i_JN = 0)
np_params = {
    "theta_S": np.pi / 4,
    "phi_S": 0.0,
    "theta_J": np.pi / 2,
    "phi_J": np.pi / 2,
    "mcz": 20 * solar_mass,
    "dist": 1.5 * giga_parsec,
    "eta": 0.234375,
    "t_c": 0.0,
    "phi_c": 0.0,
    "theta_tilde": 0.0,
    "Omega_tilde": 0.0,
    "gamma_P": 0.0,
}