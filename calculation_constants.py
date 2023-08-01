#----------------------------------------------------------------------------
# Name:	Calculation Constants
#----------------------------------------------------------------------------
# Fractional water contents
WC_PLASMA = 0.94 # plasma
WC_RBC = 0.65 # red blood cells

# Constants from Dash and Bassingthwaight 2010
N0 = 1.7 # unitless
K_1 = 7.43e-7 # M
K_2 = 2.95*10**-5 # unitless
K_3 = 2.51*10**-5 # unitless
K_2_PRIME_1 = 5.5*10**-4 # M
K_2_PRIME_2 = 1*10**-6 # M
K_2_PRIME_3 = 1*10**-6 # M
K_2_PRIME_5 = 2.63*10**-8 # M
K_2_PRIME_6 = 1.91*10**-8 # M
K_PRIME_1 = 1.35*10**-3 # unitless
K_PRIME_2 = K_2*K_2_PRIME_2**-1 # M**-1
K_PRIME_3 = K_3*K_2_PRIME_3**-1 # M**-1

# Constants from Wagner & Pruss 1993
TEMP_CRITICAL = 647.096 # K
PRES_CRITICAL = 22.064e3 # kPa
A1 = -7.85951783
A2 = 1.84408259
A3 = -11.7866497
A4 = 22.6807411
A5 = -15.9618719
A6 = 1.80122502

# Standard temperature and pressure
STD_R = 8.3145 # J.K**-1.mol**-1
STD_TEMP = 273.15 # K
STD_PRES = 101.325 # kPa
STD_R_TEMP_PRES_1 = STD_R*STD_TEMP*STD_PRES**-1
STD_R_TEMP_PRES_1_1E2 = STD_R*STD_TEMP*STD_PRES**-1*1e2

# Standard bicarbonate
STD_BICARB = 24.5e-3 # M

# Standard unit conversion values
STD_UNIT_CONV_DICT={
	'kg': 1,
	'lb': 0.45359237,
	'm': 1,
	'ft': 0.3048,
	'K': 1,
	'mmol/l': 1e-3,
	'mol/l': 1,
	'mEq/l': 1e-3,
	'Eq/l': 1,
	'g/dl': 1,
	'g/l': 1,
	'fraction': 1,
	'%': 1e-2,
	'unitless': 1,
	'ml': 1e-3,
	'l': 1,
	'l/min': 1,
	'ml/min': 10**-3,
	'bpm': 1,
	'mlO2/min/kPa': 1,
}