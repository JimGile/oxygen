﻿#!/usr/local/shared_bin/python/python-2.7.2/ES5-64/bin/python
# -*- coding: UTF-8 -*-
#----------------------------------------------------------------------------
# Name:		 	Oxygen delivery model
# Authors:  	Andrew Bretherick, Kenneth Baillie
# Copyright:	Bretherick, Baillie 2018
# Contact:      j.k.baillie@baillielab.net
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Name:		Imported functions
#----------------------------------------------------------------------------
from scipy import optimize, integrate, special
from math import exp, e, log10, log
import numpy as np
import sys, timeit

module = sys.modules[__name__]

MAX_RUNS=2
TOL_ERR_SO2 =	1e-7
TOL_ERR_LUNG =	0.00001
TOL_ERR_ORGAN = 0.00001
TOL_ERR_PH =	0.00001
VO2_CORRECTION_SPEED = 3 # normally 3; faster has higher chance of oscillation
GLOBAL_DIFF_TOL = 0.045

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

#========================
RESET_VAR_LIST=[
	'PAO2',
	'PcO2',
	'PaO2',
	'PtO2',
	'PvO2',
	'CcO2',
	'CaO2',
	'CtO2',
	'CvO2',
	'SaO2',
	'PACO2',
	'PcCO2',
	'PaCO2',
	'PtCO2',
	'PvCO2',
	'CcCO2',
	'CaCO2',
	'CtCO2',
	'CvCO2',
	'pH_c',
	'pH_a',
	'pH_t',
	'pH_v',
	'HCO3_c',
	'HCO3_a',
	'HCO3_t',
	'HCO3_v',
	'trueVO2',
	'P50_a',
	'P50_v',
	]
RESET_VAR_DICT = {x:1 for x in RESET_VAR_LIST}
OLD_VALUES_DICT={}

CcCO2 = 1
CcO2 = 1
CvO2 = 1
CvCO2 = 1
CaCO2 = 1
CaO2 = 1
CtCO2 = 1
CtO2 = 1


# Standard unit conversion values
UNIT_CONV_VALS={
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

#----------------------------------------------------------------------------
# Name:		 Unit Conversion
# Purpose:	Convert to SI units
#----------------------------------------------------------------------------
def get_standard_unit_value(value,unit):
	if unit == 'C':
		return value+273.15
	elif unit == 'F':
		return (5*(value-32)*9**-1)+273.15
	else:
		return value*UNIT_CONV_VALS[unit]

#----------------------------------------------------------------------------
# Name:	Standard User input values
#----------------------------------------------------------------------------
RR=get_standard_unit_value(float(12.0),'bpm')
VT=get_standard_unit_value(float(0.475),'l')
VD=get_standard_unit_value(float(0.11),'l')
fio2=get_standard_unit_value(float(0.21),'fraction')
alt=get_standard_unit_value(float(8200),'m')
CO=get_standard_unit_value(float(6.5),'l/min')
pulm_shunt=get_standard_unit_value(float(0.03),'fraction')
DmO2=get_standard_unit_value(float(300),'mlO2/min/kPa')
Vc=get_standard_unit_value(float(0.075),'l')
Hb=get_standard_unit_value(float(150),'g/l')
BE=get_standard_unit_value(float(0),'mEq/l')
DPG=get_standard_unit_value(float(4.65),'mmol/l')
MCHC=get_standard_unit_value(float(340),'g/l')
VO2=get_standard_unit_value(float(0.25),'l/min')
Temp=get_standard_unit_value(float(309.65),'K')
tissue_shunt=get_standard_unit_value(float(0.05),'fraction')
RQ=get_standard_unit_value(float(0.8),'fraction')

ecmosites=[]
Qecmo=get_standard_unit_value(0,'l/min')
hetindex=0


#----------------------------------------------------------------------------
# Name:		Calculated Constants
#----------------------------------------------------------------------------
def calculatedconstants():
	global Pres,PH2O,PIO2
	global Wbl
	global VA, VQ
	global alphaO2,alphaCO2
	global HbMol, Hct
	global BE, trueVO2
	global PAO2, PACO2
	#------------------------------------------------------------------------
	# Name:		 Atmospheric pressure (kPa) from altitude (m)
	# Source:	West 1996
	#------------------------------------------------------------------------
	Pres = exp(6.63268-0.1112*(alt*1e-3)-0.00149*(alt*1e-3)**2)*0.1333 # kPa
	#------------------------------------------------------------------------
	# Name:		 Saturated vapour pressure of water (kPa)
	# Source:	   Wagner & Pruss 1993
	#------------------------------------------------------------------------
	tau = 1 - Temp*TEMP_CRITICAL**-1 # fraction
	PH2O = PRES_CRITICAL*e**(TEMP_CRITICAL*Temp**-1*(A1*tau+A2*tau**1.5+A3*\
		tau**3+A4*tau**3.5+A5*tau**4+A6*tau**7)) # kPa
	#------------------------------------------------------------------------
	# Name:		 PIO2
	#------------------------------------------------------------------------
	PIO2 = fio2*(Pres-PH2O) # kPa
	#------------------------------------------------------------------------
	# Name:		 [Hb]
	#------------------------------------------------------------------------
	HbMol = Hb*64458**-1 # M
	Hct  = float(Hb)/MCHC
	#------------------------------------------------------------------------
	# Name:		 Fractional water space of blood
	#------------------------------------------------------------------------
	Wbl = (1-Hct)*WC_PLASMA+Hct*WC_RBC # fraction
	#------------------------------------------------------------------------
	# Name:		 Alveolar Ventilation
	#------------------------------------------------------------------------
	VA = RR*(VT-VD) # l/min
	VQ = VA*(CO*(1-pulm_shunt))**-1
	#------------------------------------------------------------------------
	# Name:		 alpha O2 / CO2
	# Source:	   Dash & Bassingthwaight 2010
	#------------------------------------------------------------------------
	alphaO2 = alpha_o2_func(Temp, WC_PLASMA) # M/kPa
	alphaCO2 = alpha_co2_func(Temp, WC_PLASMA) # M/kPa
	#------------------------------------------------------------------------
	#----------- AND FIX INPUT UNITS
	BE = BE/1000 # input is mEq/l
	trueVO2 = VO2 # user-defined VO2 is an aspiration...
	#------------------------------------------------------------------------
	PACO2 = trueVO2*RQ*STD_PRES*Temp*STD_TEMP**-1*VA**-1 # Alveolar ventilation equation. Simplifies to PACO2=0.863*VCO2/VA under normal conditions
	PAO2 = PIO2 - ((PACO2*(1-fio2*(1-RQ)))*RQ**-1)# Alveolar gas Equation
	#------------------------------------------------------------------------
	#------------------------------------------------------------------------

#-###########################################################################
# Name:		 FUNCTIONS
#-###########################################################################
#----------------------------------------------------------------------------
# Name:			alpha_o2_func
# Source:		Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def alpha_o2_func(temp, wpl):
	return 0.1333**-1*(1.37-0.0137*(temp-310.15)+0.00058*(temp-310.15)**2)*\
		(1e-6*wpl**-1) # M/kPa
#----------------------------------------------------------------------------
# Name:		   alpha_co2_func
# Source:	   Kelman 1967
#----------------------------------------------------------------------------
def alpha_co2_func(temp, wpl):
	return 0.1333**-1*(3.07-0.057*(temp-310.15)+0.002*(temp-310.15)**2)*\
		(1e-5*wpl**-1) # M/kPa
#----------------------------------------------------------------------------
# Name:		   p50_func
# Source:	   Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def p50_func(ph, pn_co2, dpg, temp):
	p50_ph = -2.84*(ph+log10(0.69)-7.24)+1.18*(ph+log10(0.69)-7.24)**2
	p50_pn_co2 = 4.82e-2*(pn_co2-5.332)+3.64e-5*(pn_co2-5.332)**2
	p50_dpg = 1.06e2*(dpg-4.65e-3)-2.62e3*(dpg-4.65e-3)**2
	p50_temp = 1.99e-1*(temp-310.15)+5.78e-3*(temp-310.15)**2+9.33e-05*(temp-310.15)**3
	return 3.57+p50_pn_co2+p50_ph+p50_dpg+p50_temp # kPa
#----------------------------------------------------------------------------
# Name:		 SnO2
# Source:	   Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def sn_o2_f1(pn_o2, ph, pn_co2, dpg, temp): # Equation B.3
	'''Calculate O2 SATURATION from PARTIAL PRESSURE''' # fraction
	try:
		pno2p50 = (pn_o2*p50_func(ph, pn_co2, dpg, temp)**-1)**(1+N0)
		return (pno2p50)/(1+(pno2p50))
	except:
		print(pn_o2, ph, pn_co2, dpg, temp)
		print("failed at sn_o2_f1")
		sys.exit()
def sn_o2_f2_null(sats, cn_o2, p50_sn_o2):
	'''returns 0 '''
	return Wbl*alphaO2*p50_sn_o2*(sats*(1-sats)**-1)**((1+N0)**-1) + \
		(4*HbMol)*sats - (cn_o2*(STD_R*STD_TEMP*STD_PRES**-1*1e2)**-1)
def sn_o2_f2(cn_o2, p50_sn_o2):
	'''Calculate O2 SATURATION from CONTENT''' # fraction
	cn_o2 = max(cn_o2,0.1) # prevent negative values being fed to acidbase2
	p50_sn_o2 = max(p50_sn_o2,0.1) # prevent negative values being fed to acidbase2
	return optimize.brentq(sn_o2_f2_null,1e-15,(1-1e-15),args=(cn_o2, p50_sn_o2),rtol=TOL_ERR_SO2)
#----------------------------------------------------------------------------
# Name:		 Blood O2 content
# Source:	   Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def cn_o2_f1(pn_o2, ph, pn_co2, dpg, temp):
	'''Calculate O2 CONTENT from PARTIAL PRESSURE'''
	return (Wbl*p50_func(ph,pn_co2,dpg,temp)*(sn_o2_f1(pn_o2,ph,pn_co2,dpg,temp)*(1\
	-sn_o2_f1(pn_o2,ph,pn_co2,dpg,temp))**-1)**((1+N0)**-1)*alphaO2\
	+4*HbMol*sn_o2_f1(pn_o2,ph,pn_co2,dpg,temp))*(STD_R*STD_TEMP*STD_PRES**-1*1e2)
	# ml of O2 per 100ml blood STP
#----------------------------------------------------------------------------
# Name:		 Blood O2 partial pressure
# Source:	   Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def pn_o2_f1(sn_o2, p50_sn_o2):
	'''Calculate O2 PARTIAL PRESSURE from SATURATION'''
	return p50_sn_o2*(sn_o2*(1-sn_o2)**-1)**((1+N0)**-1) # kPa
def pn_o2_f2(cn_o2, p50_sn_o2):
	'''Calculate O2 PARTIAL PRESSURE from CONTENT'''
	sn_o2 = sn_o2_f2(cn_o2, p50_sn_o2)
	return p50_sn_o2*(sn_o2*(1-sn_o2)**-1)**((1+N0)**-1) # kPa
#----------------------------------------------------------------------------
# Name:		 CnCO2
# Source:	   Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def k_ratio_func(pn_co2, ph, temp, wpl):
	hrbc = 10**-(ph+log(0.69,10)) # M
	CO2 = pn_co2*alpha_co2_func(temp, wpl) # M
	return (K_PRIME_2*CO2*(1+K_2_PRIME_2*hrbc**-1)
	 	+(1+hrbc*K_2_PRIME_5**-1))*(K_PRIME_3*CO2*(1+K_2_PRIME_3*hrbc**-1)
		+(1+hrbc*K_2_PRIME_6**-1))**-1 # unitless
def cn_co2_dissolved(pn_co2):
	'''Calculate CO2 dissolved in WHOLE BLOOD'''
	return Wbl*alphaCO2*pn_co2 # M
def cn_co2_bicarb(ph, pn_co2):
	'''Calculate CO2 as bicarbonate in WHOLE BLOOD'''
	return ((1-Hct)*WC_PLASMA+Hct*WC_RBC*0.69)*hh_f1(pn_co2, ph) # M
def cn_co2_hb_bound(pn_co2, pn_o2, ph):
	'''Calculate Hb-CO2 from PARTIAL PRESSURE'''
	k_prime_4 = (alphaO2*pn_o2)**N0*k_ratio_func(pn_co2,ph,Temp,WC_PLASMA)*(p50_func(ph, pn_co2, DPG, Temp)*alphaO2)**-(1+N0)
	s_hb_co2 = ( \
	K_PRIME_2*alphaCO2*pn_co2*(1+K_2_PRIME_2*(10**-(ph+log(0.69,10)))**-1) + \
	K_PRIME_3*alphaCO2*pn_co2*(1+K_2_PRIME_3*(10**-(ph+log(0.69,10)))**-1)*k_prime_4*alphaO2*pn_o2
	)*( \
	K_PRIME_2*alphaCO2*pn_co2*(1+K_2_PRIME_2*(10**-(ph+log(0.69,10)))**-1) + \
	(1+(10**-(ph+log(0.69,10)))*K_2_PRIME_5**-1) + \
	k_prime_4*alphaO2*pn_o2*(K_PRIME_3*alphaCO2*pn_co2*(1+K_2_PRIME_3*(10**-(ph+log(0.69,10)))**-1) + \
	(1+(10**-(ph+log(0.69,10)))*K_2_PRIME_6**-1)) \
	)**-1
	return 4*HbMol*s_hb_co2 # M
def cn_co2_f1(ph, pn_co2, pn_o2):
	'''Calculate CO2 CONTENT from PARTIAL PRESSURE'''
	return (cn_co2_hb_bound(pn_co2,pn_o2,ph)+cn_co2_bicarb(ph,pn_co2)\
		+cn_co2_dissolved(pn_co2))*(STD_R*STD_TEMP*STD_PRES**-1*1e2) # ml CO2 per 100ml blood STP
def pn_co2_f1_null(pn_co2, cn_co2, ph, cn_o2):
	'''returns 0'''
	p50_pn_co2_null = p50_func(ph, pn_co2, DPG, Temp)
	return cn_co2_f1(ph, pn_co2, pn_o2_f2(cn_o2, p50_pn_co2_null))-cn_co2 # null
def pn_co2_f1(cn_co2, ph, cn_o2):
	'''Calculate CO2 PARTIAL PRESSURE from CONTENTS'''
	args = cn_co2,ph,cn_o2
	pn_co2 = optimize.newton(pn_co2_f1_null,AGE['PcCO2'],args=args,tol=0.01)
	return pn_co2 # kPa
#----------------------------------------------------------------------------
# Name:	 Henderson-Hasselbalch Equation
#----------------------------------------------------------------------------
def hh_f1(pn_co2, ph):
	'''Calculate CO2 as bicarbonate in SOLUTION'''
	return (K_1*alphaCO2*pn_co2)*(10**-ph)**-1 # M
#----------------------------------------------------------------------------
# Name:	 van Slyke Equation
# Source:   Siggaard-Andersen 1977
#----------------------------------------------------------------------------
def van_slyke_f1(ph, base_excess, sn_o2):
	'''returns bicarbonate from PLASMA pH @ 37 deg C and BE'''
	zeta = 1-(0.0143*Hb*1e-1)
	beta = 9.5+1.63*Hb*1e-1
	return ((base_excess*1e3 - 0.2*Hb*1e-1*(1-sn_o2))*zeta**-1 - beta*(ph-7.4) + STD_BICARB*1e3)*1e-3 # M
#----------------------------------------------------------------------------
# Name:	 simultaneous solution Henderson-Hasselbalch and van Slyke
#----------------------------------------------------------------------------
def acidbase_f1_null(ph,base_excess, pn_co2, pn_o2):
	'''returns 0'''
	sn_o2 = sn_o2_f1(pn_o2, ph, pn_co2, DPG, Temp)
	zeta = 1-(0.0143*Hb*1e-1)
	beta = 9.5+1.63*Hb*1e-1
	return ((base_excess*1e3 - 0.2*Hb*1e-1*(1-sn_o2))*zeta**-1 - beta*(ph-7.4) + STD_BICARB*1e3)*1e-3 - hh_f1(pn_co2, ph) # null

def acidbase_f1(base_excess, pn_co2, pn_o2):
	'''returns PLASMA pH from PARTIAL PRESSURE '''
	return optimize.brentq(acidbase_f1_null,1,14,args=(base_excess,pn_co2,pn_o2),rtol=TOL_ERR_PH) # pH units

def acidbase_f2_null(ph, base_excess, cn_co2, cn_o2):
	'''returns 0'''
	pn_co2 = pn_co2_f1(cn_co2, ph, cn_o2)
	sn_o2 = sn_o2_f2(cn_o2, p50_func(ph, pn_co2, DPG, Temp))
	zeta = 1-(0.0143*Hb*1e-1)
	beta = 9.5+1.63*Hb*1e-1
	return ((base_excess*1e3 - 0.2*Hb*1e-1*(1-sn_o2))*zeta**-1 - beta*(ph-7.4) + STD_BICARB*1e3)*1e-3 - hh_f1(pn_co2, ph) # null

def acidbase_f2(base_excess, cn_co2, cn_o2):
	'''returns PLASMA pH from CONTENT '''
	# prevent negative values being fed to acidbase_f2
	cn_co2 = max(cn_co2,0.1) 
	cn_o2 = max(cn_o2,0.1)
	return optimize.brentq(acidbase_f2_null,1,14,args=(base_excess, cn_co2, cn_o2),rtol=TOL_ERR_PH) # pH units

#-###########################################################################
# Name:		 COMPARTMENT SPECIFIC FUNCTIONS
#-###########################################################################
#-------------------------------------------------------------------------------
#		   Alveolar Gas Equation
#-------------------------------------------------------------------------------
def populate_alv_gas_eqn():
	global AGE; AGE={}
	global PACO2, PAO2
	AGE['PACO2'] = VO2*RQ*STD_PRES*Temp*STD_TEMP**-1*VA**-1 # Alveolar ventilation equation
	AGE['PcCO2'] = AGE['PACO2'] # assumes complete equilibrium
	AGE['PAO2'] = PIO2 - ((AGE['PACO2']*(1-fio2*(1-RQ)))*RQ**-1)# Alveolar gas Equation
	AGE['PcO2'] = AGE['PAO2'] # assumes complete equilibrium
	AGE['pH'] = acidbase_f1(BE,AGE['PcCO2'],AGE['PcO2'])
	AGE['CcCO2'] = cn_co2_f1(AGE['pH'], AGE['PcCO2'], AGE['PcO2'])
	AGE['CcO2'] = cn_o2_f1(AGE['PcO2'], AGE['pH'], AGE['PcCO2'], DPG, Temp)
	AGE['CvO2'] = AGE['CcO2'] - 100*(VO2*(CO*(1-pulm_shunt))**-1) 
	AGE['CvCO2'] = AGE['CcCO2'] + 100*(VO2*RQ)*(CO*(1-pulm_shunt))**-1
	PAO2 = AGE['PAO2']
	PACO2 = AGE['PACO2']
	#print(AGE)
#-------------------------------------------------------------------------------
# Name:		 Transit Time, initial value problem
# Source:	   Wagner and West 1972
#-------------------------------------------------------------------------------
def dxdt(inputs, t, vq, dm_o2_ivp, vc_ivp, cv_o2, cv_co2):
	cc_o2 = inputs[0]
	cc_co2 = inputs[1]
	pa_co2 = 0
	pa_o2 = fio2*(Pres-PH2O)
	if (cv_co2 - cc_co2) != (cc_o2 - cv_o2):
		rq = (cv_co2-cc_co2)*(cc_o2-cv_o2)**-1 # ratio
		pa_co2 = vq**-1*(cc_o2-cv_o2)*1e-2*rq*STD_PRES*Temp*STD_TEMP**-1 # kPa
		pa_o2 =  fio2*(Pres-PH2O)-((pa_co2*(1-fio2*(1-rq)))*rq**-1) # kPa
	ph_c = acidbase_f2(BE, cc_co2, cc_o2) # pH units
	pc_co2 = pn_co2_f1(cc_co2, ph_c, cc_o2) # kPa
	p50_dxdt = p50_func(ph_c, pc_co2, DPG, Temp) # kPa
	pc_o2 = pn_o2_f2(cc_o2, p50_dxdt) # kPa
	sats = sn_o2_f1(pc_o2, ph_c, pc_co2, DPG, Temp) # fraction
	k_prime_c = 1.25284e5+3.6917e4*e**(3.8200*sats) # M**-1.sec**-1
	rate_o2 = k_prime_c*alphaO2*60*(1-sats)*4*HbMol*STD_R*STD_TEMP*STD_PRES**-1 # ml(O2).ml(bld)**-1.kPa**-1.min**-1
	dl_o2 = (dm_o2_ivp**-1 + (rate_o2*vc_ivp*1e3)**-1)**-1 # ml(O2).min**-1.kPa**-1
	do2_dt = 100*(vc_ivp*1e3)**-1*dl_o2*(pa_o2-pc_o2) # ml(O2).100ml(bld)**-1.min**-1
	dm_co2 = dm_o2_ivp*20 # ml(CO2).min**-1.kPa**-1
	dl_co2 = dm_co2 # assumes infinitely fast rate of reaction of CO2
	dco2_dt = 100*(vc_ivp*1e3)**-1*dl_co2*(pa_co2-pc_co2) # ml.CO2.100mlBlood^-1.min^-1
	return [do2_dt, dco2_dt]
#-------------------------------------------------------------------------------
# Name:		 Transit Time, initial value problem
# Source:	   Wagner and West 1972, altered 
#-------------------------------------------------------------------------------
def dxdt_organ(inputs, t, dm_o2_ivp, vol_b, c_in_o2, c_in_co2, q, compound, qd, c_di_o2, c_di_co2):
	# no diffusion, no change. 
	if dm_o2_ivp==0:
		return [0,0] 
	# these are the only values that are changed with each iteration
	c_out_o2 = inputs[0] 
	c_out_co2 = inputs[1]
	c_do_o2 = c_di_o2
	c_do_co2 = c_di_co2
	if (c_in_co2 - c_out_co2) != (c_out_o2 - c_in_o2): 
		c_do_o2 = (qd*c_di_o2 - q*(c_out_o2-c_in_o2))/qd 		# FICK within organ
		c_out_o2 = (q*c_in_o2 + qd*(c_di_o2-c_do_o2))/q 		# FICK in blood 
		c_do_co2 = (qd*c_di_co2 - q*(c_out_co2-c_in_co2))/qd 	# FICK within organ
		c_out_co2 = (q*c_in_co2 + qd*(c_di_co2-c_do_co2))/q 	# FICK in blood
	p_organ_o2 = (c_do_o2*10*1000**-1) * STD_PRES 
	k = STD_PRES*Temp/STD_TEMP
	p_organ_co2 = (c_do_co2*10*1000**-1) * k 
	# equilibrate haemoglobin
	ph_c = acidbase_f2(BE, c_out_co2, c_out_o2) # pH units
	p_out_co2 = pn_co2_f1(c_out_co2, ph_c, c_out_o2) # kPa
	p50_dxdt = p50_func(ph_c, p_out_co2, DPG, Temp) # kPa
	p_out_o2 = pn_o2_f2(c_out_o2, p50_dxdt) # kPa
	sats = sn_o2_f1(p_out_o2, ph_c, p_out_co2, DPG, Temp) # fraction
	# now do O2
	k_prime_c = 1.25284e5+3.6917e4*e**(3.8200*sats) # M**-1.sec**-1
	theta = k_prime_c*alphaO2*60*(1-sats)*4*HbMol*STD_R*STD_TEMP*STD_PRES**-1 # ml(O2).ml(bld)**-1.kPa**-1.min**-1
	diffusion_o2 = (dm_o2_ivp**-1 + (theta*vol_b*1e3)**-1)**-1 # ml(O2).min**-1.kPa**-1
	do2_dt = 100*(vol_b*1e3)**-1*diffusion_o2*(p_organ_o2-p_out_o2) # ml(O2).100ml(bld)**-1.min**-1
	# now do CO2
	dm_co2 = dm_o2_ivp*20 # ml(CO2).min**-1.kPa**-1
	if dm_co2==0:
		return [0,0] # no diffusion, no change. 
	dl_co2 = dm_co2 # assumes infinitely fast rate of reaction of CO2
	dco2_dt = 100*(vol_b*1e3)**-1*dl_co2*(p_organ_co2-p_out_co2) # ml.CO2.100mlBlood^-1.min^-1
	return [do2_dt, dco2_dt]

def run_organ(dtype, preorgan_cn_o2, preorgan_cn_co2):
	if dtype=="lung":
		# organ settings for organ==LUNG
		organgas = "air"
		organflow = VA # organ flow rate l/min
		v_organ_blood = Vc
		q_organ_blood = CO*(1-pulm_shunt)
		organ_o2_in_content = fio2*(Pres-PH2O) # Content organ input O2 mls_gas/volume
		organ_co2_in_content = 0 # Content organ input CO2 mls_gas/volume
		membrane_diffusion = DmO2
		organtime = v_organ_blood/q_organ_blood # in minutes
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		post_organ_cn_o2, post_organ_cn_co2 = integrate.odeint(dxdt_organ,\
				[preorgan_cn_o2,preorgan_cn_co2],[0,organtime],\
				(membrane_diffusion, v_organ_blood, preorgan_cn_o2, preorgan_cn_co2, q_organ_blood,\
				organgas, organflow, organ_o2_in_content, organ_co2_in_content ),\
				rtol=TOL_ERR_ORGAN, mxstep=1000)[1]
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		# and get the organ values too:
		c_do_o2 = (organflow*organ_o2_in_content - q_organ_blood*(post_organ_cn_o2-preorgan_cn_o2))/organflow 		# FICK within organ
		c_do_co2 = (organflow*organ_co2_in_content - q_organ_blood*(post_organ_cn_co2-preorgan_cn_co2))/organflow
		post_organ_o2 = (c_do_o2*10*1000**-1) * STD_PRES*Temp/STD_TEMP 
		post_organ_co2 = (c_do_co2*10*1000**-1) * STD_PRES*Temp/STD_TEMP 
		return post_organ_cn_o2, post_organ_cn_co2, post_organ_o2, post_organ_co2

#-------------------------------------------------------------------------------
#		   mass balance
#-------------------------------------------------------------------------------
def lung_null(cv_o2__cv_co2):
	VQ = VA*(CO*(1-pulm_shunt))**-1
	[content['CvO2'],content['CvCO2']] = cv_o2__cv_co2
	[content['CcO2'],content['CcCO2']] = \
	integrate.odeint(dxdt,[content['CvO2'],content['CvCO2']],[0,Vc*(CO*(1-pulm_shunt))**-1],\
	(VQ,DmO2,Vc,content['CvO2'],content['CvCO2']),rtol=TOL_ERR_LUNG)[1]
	out = [content['CcO2'] - content['CvO2'] - 100*((trueVO2)*(CO*(1-pulm_shunt))**-1)]
	out.append(content['CvCO2'] - content['CcCO2'] - 100*((RQ*trueVO2)*(CO*(1-pulm_shunt))**-1))
	return out
#-------------------------------------------------------------------------------
#		 compartment contents
#-------------------------------------------------------------------------------
def updatebloodgascontents():
	global content; content = {} # this line starts again at the beginning. 
	#-------------------------------------------------------------------------------
	# Pulmonary capillaries & Veins
	#-------------------------------------------------------------------------------
	optimize.fsolve(lung_null, [AGE['CvO2'],AGE['CvCO2']], xtol=0.001)
	#-------------------------------------------------------------------------------
	# Arteries
	#-------------------------------------------------------------------------------
	content['CaCO2'] = pulm_shunt*content['CvCO2'] + (1-pulm_shunt)*content['CcCO2'] # mlO2/100mlblood
	content['CaO2'] = pulm_shunt*content['CvO2'] + (1-pulm_shunt)*content['CcO2'] # mlO2/100mlblood
	#-------------------------------------------------------------------------------
	# Tissues
	#-------------------------------------------------------------------------------
	QO2 = CO*(1-tissue_shunt)*content['CaO2']*100**-1 # lO2/min
	content['CtO2'] = 100*(QO2-VO2)/((1-tissue_shunt)*CO) # mlO2/100mlblood
	content['CtCO2'] = content['CaCO2'] + 100*VO2*RQ*(CO*(1-tissue_shunt))**-1 # ml.CO2/100ml.Blood
	#-------------------------------------------------------------------------------
	# copy to globals
	#-------------------------------------------------------------------------------
	# print("content", content)
	for name, value in content.items():
		setattr(module, name, value)

#-------------------------------------------------------------------------------
#		   compartment partial pressures, acid base and O2 saturation
#-------------------------------------------------------------------------------
def updatepartialpressures(compartments = ['c','a','t','v']):
	for i in compartments:
		setattr(module,"pH_%s"%i, eval("acidbase_f2(BE, C%sCO2, C%sO2)"%(i,i)))
		setattr(module,"P%sCO2"%i, eval("pn_co2_f1(C%sCO2, pH_%s, C%sO2)"%(i,i,i)))
		setattr(module,"P50_%s"%i, eval("p50_func(pH_%s, P%sCO2, DPG, Temp)"%(i,i)))
		setattr(module,"S%sO2"%i, eval("sn_o2_f2(C%sO2, P50_%s)"%(i,i)))
		setattr(module,"P%sO2"%i, eval("pn_o2_f1(S%sO2, P50_%s)"%(i,i)))
		setattr(module,"HCO3_%s"%i, eval("van_slyke_f1(pH_%s, BE, S%sO2)"%(i,i)))

#-------------------------------------------------------------------------------
# Name:		 Timecheck
#-------------------------------------------------------------------------------
def fail(message=''):
	print(0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",message)
	sys.exit()

def check_completion():
	old_vals_len = len(OLD_VALUES_DICT.values())
	global_dict = globals()
	diff_dict={}
	for key in RESET_VAR_LIST:
		if old_vals_len > 0:
			diff = abs(float(global_dict[key]) - float(OLD_VALUES_DICT[key]))
			if diff > GLOBAL_DIFF_TOL:
				diff_dict[key] = diff
		OLD_VALUES_DICT[key] = global_dict[key]
	if old_vals_len < 1:
		return 100
	return sum(diff_dict.values())

def formatoutput():
	if CtO2 <= 0:
		fail("Ct02<0")
	else: # online
		outputdata = [
				['PatmosO2',round(fio2*Pres,2)],
				['PIO2',round(fio2*(Pres - PH2O),2)],
				['PAO2',PAO2],
				['PcO2',PcO2],
				['PaO2',PaO2],
				['PtO2',PtO2],
				['PvO2',PvO2],
				['PaO2',round(PaO2,1)],
				['PaCO2',round(PaCO2,1)],
				['pH',round(pH_a,2)],
				['H+',round(10**9*10**-pH_a,2)],
				['HCO3',round(HCO3_a*1000,1)],
				['SaO2',round(SaO2*100,3)],
				]
		print('||'.join([str(x[1]) for x in outputdata]))


def format_output_json():
	if CtO2 <= 0:
		fail("Ct02<0")
	else: # online
		output_data = {
				'PatmosO2':round(fio2*Pres,4),
				'PIO2':round(fio2*(Pres - PH2O),4),
				'PAO2':round(PAO2,4),
				'PcO2':round(PcO2,4),
				'PtO2':round(PtO2,4),
				'PvO2':round(PvO2,4),
				'PaO2':round(PaO2,4),
				'PaCO2':round(PaCO2,4),
				'pH':round(pH_a,2),
				'H+':round(10**9*10**-pH_a,2),
				'HCO3':round(HCO3_a*1000,2),
				'SaO2':round(SaO2*100,2),
				}
		print(output_data)

#-------------------------------------------------------------------------------
# Name:		 Run model
#-------------------------------------------------------------------------------
def circulate_once(iterationnumber=0):
	global PAO2,PcO2,PaO2,PtO2,PvO2,CcO2,CaO2,CtO2,CvO2,SaO2,PACO2,PcCO2,PaCO2,PtCO2,PvCO2,CcCO2,CaCO2,CtCO2,CvCO2
	global pH_c,pH_a,pH_t,pH_v,HCO3_c,HCO3_a,HCO3_t,HCO3_v, trueVO2
	# ====== alveoli and pulmonary capillaries =======
	CcO2, CcCO2, PAO2, PACO2 = run_organ('lung', CvO2, CvCO2)
	updatepartialpressures(compartments = ['c'])
	# ======== arteries ========
	CaO2 = pulm_shunt*CvO2 + (1-pulm_shunt)*CcO2
	CaCO2 = pulm_shunt*CvCO2 + (1-pulm_shunt)*CcCO2 # mlO2/100mlblood
	updatepartialpressures(compartments = ['a'])
	# ========= tissues ==========
	q_tissue = CO*(1-tissue_shunt)
	CtO2 = (q_tissue*CaO2*10 - trueVO2*1000)/ (q_tissue*10)
	CtCO2 = (q_tissue*CaCO2*10 + trueVO2*RQ*1000)/ (q_tissue*10)
	updatepartialpressures(compartments = ['t'])
	# ========= veins ==========
	CvO2 = (CtO2*q_tissue + CaO2*CO*tissue_shunt)/CO
	CvCO2 = (CtCO2*q_tissue + CaCO2*CO*tissue_shunt)/CO
	updatepartialpressures(compartments = ['v'])
	# ========= set VO2 ==========
	critical_oer = 0.94 # unrealistic maximum
	critical_do2 = (CO*CaO2*10)*critical_oer/1000 # mlsO2/min
	if critical_do2 < VO2: #then set trueVO2 = criticalDO2, but prevent oscillation:
		trueVO2 = trueVO2 - float(trueVO2-critical_do2)/(iterationnumber/VO2_CORRECTION_SPEED+1) # only go a hundredth of the way to avoid oscillation
		#print("CO %s CaO2 %s criticalDO2 %s < VO2 %s so correcting downwards to %s"%(CO, CaO2, criticalDO2, VO2, trueVO2)
	elif trueVO2<VO2: # then it needs to come back up
		trueVO2 = VO2

#========================
#========================

# ------  ------  ------  ------  ------  ------  ------  ------ 
# getinputs()
calculatedconstants()
populate_alv_gas_eqn()
if __name__ == "__main__":
	try:
		print("pre-updatebloodgascontents")
		updatebloodgascontents()
		print("pre-updatepartialpressures")
		updatepartialpressures()
	except BaseException as e:
		print(content)
		raise e

	for i in range(MAX_RUNS):
		circulate_once(i)
		globaldiff = check_completion()
		print(i, globaldiff)
		if globaldiff < GLOBAL_DIFF_TOL:
			break
	format_output_json()

# ------  ------  ------  ------  ------  ------  ------  ------ 
