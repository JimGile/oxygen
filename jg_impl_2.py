#!/usr/local/shared_bin/python/python-2.7.2/ES5-64/bin/python
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
debugging = False

toleranceoferror_so2 =	1e-7
toleranceoferror_lung =	0.00001
toleranceoferror_organ = 0.00001
toleranceoferror_ph =	0.00001
VO2correctionspeed =	3 # normally 3; faster has higher chance of oscillation
globaldifftolerance =	0.01
num_lung_compartments = 20

#------------
Pbar = 101.325 # the barometric pressure at which these solubility coefficients apply
#------------

#----------------------------------------------------------------------------
# Name:	Calculation Constants
#----------------------------------------------------------------------------
# Fractional water contents
Wpl = 0.94 # plasma
Wrbc = 0.65 # red blood cells

# Constants from Dash and Bassingthwaight 2010
n0 = 1.7 # unitless
K_1 = 7.43e-7 # M
K_prime_1 = 1.35*10**-3 # unitless
K_2prime_1 = 5.5*10**-4 # M
K_2 = 2.95*10**-5 # unitless
K_2prime_2 = 1*10**-6 # M
K_prime_2 = K_2*K_2prime_2**-1 # M**-1
K_3 = 2.51*10**-5 # unitless
K_2prime_3 = 1*10**-6 # M
K_prime_3 = K_3*K_2prime_3**-1 # M**-1
K_2prime_5 = 2.63*10**-8 # M
K_2prime_6 = 1.91*10**-8 # M

# Constants from Wagner & Pruss 1993
Temp_critical = 647.096 # K
Pres_critical = 22.064e3 # kPa
a1 = -7.85951783
a2 = 1.84408259
a3 = -11.7866497
a4 = 22.6807411
a5 = -15.9618719
a6 = 1.80122502

# Standard temperature and pressure
R = 8.3145 # J.K**-1.mol**-1
STP_T = 273.15 # K
STP_P = 101.325 # kPa

# Standard bicarbonate
StdBicarb = 24.5e-3 # M

# Standard unit conversion values
unit_conv={
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
		return value*unit_conv[unit]

#----------------------------------------------------------------------------
# Name:	Standard User input values
#----------------------------------------------------------------------------
RR=get_standard_unit_value(float(12.0),'bpm')
VT=get_standard_unit_value(float(0.475),'l')
VD=get_standard_unit_value(float(0.11),'l')
fio2=get_standard_unit_value(float(0.21),'fraction')
alt=get_standard_unit_value(float(0),'m')
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
maxruns=2

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
	tau = 1 - Temp*Temp_critical**-1 # fraction
	PH2O = Pres_critical*e**(Temp_critical*Temp**-1*(a1*tau+a2*tau**1.5+a3*\
		tau**3+a4*tau**3.5+a5*tau**4+a6*tau**7)) # kPa
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
	Wbl = (1-Hct)*Wpl+Hct*Wrbc # fraction
	#------------------------------------------------------------------------
	# Name:		 Alveolar Ventilation
	#------------------------------------------------------------------------
	VA = RR*(VT-VD) # l/min
	VQ = VA*(CO*(1-pulm_shunt))**-1
	#------------------------------------------------------------------------
	# Name:		 alpha O2 / CO2
	# Source:	   Dash & Bassingthwaight 2010
	#------------------------------------------------------------------------
	alphaO2 = alpha_o2_func(Temp, Wpl) # M/kPa
	alphaCO2 = alpha_co2_func(Temp, Wpl) # M/kPa
	#------------------------------------------------------------------------
	#----------- AND FIX INPUT UNITS
	BE = BE/1000 # input is mEq/l
	trueVO2 = VO2 # user-defined VO2 is an aspiration...
	#------------------------------------------------------------------------
	PACO2 = trueVO2*RQ*STP_P*Temp*STP_T**-1*VA**-1 # Alveolar ventilation equation. Simplifies to PACO2=0.863*VCO2/VA under normal conditions
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
		pno2p50 = (pn_o2*p50_func(ph, pn_co2, dpg, temp)**-1)**(1+n0)
		return (pno2p50)/(1+(pno2p50))
	except:
		print(pn_o2, ph, pn_co2, dpg, temp)
		print("failed at sn_o2_f1")
		sys.exit()
def sn_o2_f2_null(sats, cn_o2, p50_sn_o2):
	'''returns 0 '''
	return Wbl*alphaO2*p50_sn_o2*(sats*(1-sats)**-1)**((1+n0)**-1) + \
		(4*HbMol)*sats - (cn_o2*(R*STP_T*STP_P**-1*1e2)**-1)
def sn_o2_f2(cn_o2, p50_sn_o2):
	'''Calculate O2 SATURATION from CONTENT''' # fraction
	cn_o2 = max(cn_o2,0.1) # prevent negative values being fed to acidbase2
	p50_sn_o2 = max(p50_sn_o2,0.1) # prevent negative values being fed to acidbase2
	return optimize.brentq(sn_o2_f2_null,1e-15,(1-1e-15),args=(cn_o2, p50_sn_o2),rtol=toleranceoferror_so2)
#----------------------------------------------------------------------------
# Name:		 Blood O2 content
# Source:	   Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def cn_o2_f1(pn_o2, ph, pn_co2, dpg, temp):
	'''Calculate O2 CONTENT from PARTIAL PRESSURE'''
	return (Wbl*p50_func(ph,pn_co2,dpg,temp)*(sn_o2_f1(pn_o2,ph,pn_co2,dpg,temp)*(1\
	-sn_o2_f1(pn_o2,ph,pn_co2,dpg,temp))**-1)**((1+n0)**-1)*alphaO2\
	+4*HbMol*sn_o2_f1(pn_o2,ph,pn_co2,dpg,temp))*(R*STP_T*STP_P**-1*1e2)
	# ml of O2 per 100ml blood STP
#----------------------------------------------------------------------------
# Name:		 Blood O2 partial pressure
# Source:	   Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def pn_o2_f1(sn_o2, p50_sn_o2):
	'''Calculate O2 PARTIAL PRESSURE from SATURATION'''
	return p50_sn_o2*(sn_o2*(1-sn_o2)**-1)**((1+n0)**-1) # kPa
def pn_o2_f2(cn_o2, p50_sn_o2):
	'''Calculate O2 PARTIAL PRESSURE from CONTENT'''
	sn_o2 = sn_o2_f2(cn_o2, p50_sn_o2)
	return p50_sn_o2*(sn_o2*(1-sn_o2)**-1)**(((1+n0))**-1) # kPa
#----------------------------------------------------------------------------
# Name:		 CnCO2
# Source:	   Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def k_ratio_func(pn_co2, ph, temp, wpl):
	hrbc = 10**-(ph+log(0.69,10)) # M
	CO2 = pn_co2*alpha_co2_func(temp, wpl) # M
	return (K_prime_2*CO2*(1+K_2prime_2*hrbc**-1)+(1+hrbc*K_2prime_5**-1))*(K_prime_3*CO2*(1+K_2prime_3*hrbc**-1)+(1+hrbc*K_2prime_6**-1))**-1 # unitless
def CnCO2_Dissolved(PnCO2):
	'''Calculate CO2 dissolved in WHOLE BLOOD'''
	return Wbl*alphaCO2*PnCO2 # M
def CnCO2_Bicarb(pH, PnCO2):
	'''Calculate CO2 as bicarbonate in WHOLE BLOOD'''
	return ((1-Hct)*Wpl+Hct*Wrbc*0.69)*HH_1(PnCO2, pH) # M
def CnCO2_HbBound(PnCO2, PnO2, pH):
	'''Calculate Hb-CO2 from PARTIAL PRESSURE'''
	K_prime_4 = (alphaO2*PnO2)**n0*k_ratio_func(PnCO2,pH,Temp,Wpl)*(p50_func(pH, PnCO2, DPG, Temp)*alphaO2)**-(1+n0)
	SHbCO2 = ( \
	K_prime_2*alphaCO2*PnCO2*(1+K_2prime_2*(10**-(pH+log(0.69,10)))**-1) + \
	K_prime_3*alphaCO2*PnCO2*(1+K_2prime_3*(10**-(pH+log(0.69,10)))**-1)*K_prime_4*alphaO2*PnO2
	)*( \
	K_prime_2*alphaCO2*PnCO2*(1+K_2prime_2*(10**-(pH+log(0.69,10)))**-1) + \
	(1+(10**-(pH+log(0.69,10)))*K_2prime_5**-1) + \
	K_prime_4*alphaO2*PnO2*(K_prime_3*alphaCO2*PnCO2*(1+K_2prime_3*(10**-(pH+log(0.69,10)))**-1) + \
	(1+(10**-(pH+log(0.69,10)))*K_2prime_6**-1)) \
	)**-1
	return 4*HbMol*SHbCO2 # M
def CnCO2_1(pH, PnCO2, PnO2):
	'''Calculate CO2 CONTENT from PARTIAL PRESSURE'''
	return (CnCO2_HbBound(PnCO2,PnO2,pH)+CnCO2_Bicarb(pH,PnCO2)\
		+CnCO2_Dissolved(PnCO2))*(R*STP_T*STP_P**-1*1e2) # ml CO2 per 100ml blood STP
def PnCO2_null(PnCO2, CnCO2, pH, CnO2):
	'''returns 0'''
	P50_PnCO2_null = p50_func(pH, PnCO2, DPG, Temp)
	return CnCO2_1(pH, PnCO2, pn_o2_f2(CnO2, P50_PnCO2_null))-CnCO2 # null
def PnCO2_1(CnCO2, pH, CnO2):
	'''Calculate CO2 PARTIAL PRESSURE from CONTENTS'''
	args = CnCO2,pH,CnO2
	PnCO2 = optimize.newton(PnCO2_null,AGE['PcCO2'],args=args,tol=0.01)
	return PnCO2 # kPa
#----------------------------------------------------------------------------
# Name:	 Henderson-Hasselbalch Equation
#----------------------------------------------------------------------------
def HH_1(PnCO2, pH):
	'''Calculate CO2 as bicarbonate in SOLUTION'''
	return (K_1*alphaCO2*PnCO2)*(10**-pH)**-1 # M
#----------------------------------------------------------------------------
# Name:	 van Slyke Equation
# Source:   Siggaard-Andersen 1977
#----------------------------------------------------------------------------
def vanSlyke_1(pH, BE, SnO2):
	'''returns bicarbonate from PLASMA pH @ 37 deg C and BE'''
	zeta = 1-(0.0143*Hb*1e-1)
	beta = 9.5+1.63*Hb*1e-1
	return ((BE*1e3 - 0.2*Hb*1e-1*(1-SnO2))*zeta**-1 - beta*(pH-7.4) + StdBicarb*1e3)*1e-3 # M
#----------------------------------------------------------------------------
# Name:	 simultaneous solution Henderson-Hasselbalch and van Slyke
#----------------------------------------------------------------------------
def acidbase_1_null(pH,BE,PnCO2,PnO2):
	'''returns 0'''
	SnO2 = sn_o2_f1(PnO2, pH, PnCO2, DPG, Temp)
	zeta = 1-(0.0143*Hb*1e-1)
	beta = 9.5+1.63*Hb*1e-1
	return ((BE*1e3 - 0.2*Hb*1e-1*(1-SnO2))*zeta**-1 - beta*(pH-7.4) + StdBicarb*1e3)*1e-3 - HH_1(PnCO2, pH) # null
def acidbase_1(BE,PnCO2,PnO2):
	'''returns PLASMA pH from PARTIAL PRESSURE '''
	try:
		return optimize.brentq(acidbase_1_null,1,14,args=(BE,PnCO2,PnO2),rtol=toleranceoferror_ph) # pH units
	except:
		print(BE, PnCO2, PnO2)
		print(acidbase_1_null(1, BE, PnCO2, PnO2))
		print(acidbase_1_null(14, BE, PnCO2, PnO2))
		print("pH", optimize.brentq(acidbase_1_null,1,14,args=(BE,PnCO2,PnO2),rtol=toleranceoferror_ph)) # pH units
		return optimize.brentq(acidbase_1_null,1,14,args=(BE,PnCO2,PnO2),rtol=toleranceoferror_ph) # pH units
def acidbase_2_null(pH,BE,CnCO2,CnO2):
	'''returns 0'''
	PnCO2 = PnCO2_1(CnCO2,pH,CnO2)
	SnO2 = sn_o2_f2(CnO2, p50_func(pH,PnCO2,DPG,Temp))
	zeta = 1-(0.0143*Hb*1e-1)
	beta = 9.5+1.63*Hb*1e-1
	return ((BE*1e3 - 0.2*Hb*1e-1*(1-SnO2))*zeta**-1 - beta*(pH-7.4) + StdBicarb*1e3)*1e-3 - HH_1(PnCO2, pH) # null
def acidbase_2(BE,CnCO2,CnO2):
	'''returns PLASMA pH from CONTENT '''
	try:
		return optimize.brentq(acidbase_2_null,1,14,args=(BE,CnCO2,CnO2),rtol=toleranceoferror_ph) # pH units
	except:
		print(BE,CnCO2,CnO2)
		print(acidbase_2_null(1,BE,CnCO2,CnO2))
		print(acidbase_2_null(14,BE,CnCO2,CnO2))
		return optimize.brentq(acidbase_2_null,1,14,args=(BE,CnCO2,CnO2),rtol=toleranceoferror_ph) # pH units

#-###########################################################################
# Name:		 COMPARTMENT SPECIFIC FUNCTIONS
#-###########################################################################
#-------------------------------------------------------------------------------
#		   Alveolar Gas Equation
#-------------------------------------------------------------------------------
def populatealvgaseqn():
	global AGE; AGE={}
	global PACO2, PAO2
	AGE['PACO2'] = VO2*RQ*STP_P*Temp*STP_T**-1*VA**-1 # Alveolar ventilation equation
	AGE['PcCO2'] = AGE['PACO2'] # assumes complete equilibrium
	AGE['PAO2'] = PIO2 - ((AGE['PACO2']*(1-fio2*(1-RQ)))*RQ**-1)# Alveolar gas Equation
	AGE['PcO2'] = AGE['PAO2'] # assumes complete equilibrium
	AGE['pH'] = acidbase_1(BE,AGE['PcCO2'],AGE['PcO2'])
	AGE['CcCO2'] = CnCO2_1(AGE['pH'], AGE['PcCO2'], AGE['PcO2'])
	AGE['CcO2'] = cn_o2_f1(AGE['PcO2'], AGE['pH'], AGE['PcCO2'], DPG, Temp)
	AGE['CvO2'] = AGE['CcO2'] - 100*(VO2*(CO*(1-pulm_shunt))**-1) 
	AGE['CvCO2'] = AGE['CcCO2'] + 100*(((VO2*RQ))*(CO*(1-pulm_shunt))**-1)
	PAO2 = AGE['PAO2']
	PACO2 = AGE['PACO2']
	#print(AGE)
#-------------------------------------------------------------------------------
# Name:		 Transit Time, initial value problem
# Source:	   Wagner and West 1972
#-------------------------------------------------------------------------------
def dxdt(inputs, t, VQ, DmO2_ivp, Vc_ivp, CvO2_local, CvCO2_local):
	CcO2_local = inputs[0]
	CcCO2_local = inputs[1]
	if (CvCO2_local - CcCO2_local) != (CcO2_local - CvO2_local):
		RQ = (CvCO2_local-CcCO2_local)*(CcO2_local-CvO2_local)**-1 # ratio
		PACO2_local = VQ**-1*(CcO2_local-CvO2_local)*1e-2*RQ*STP_P*Temp*STP_T**-1 # kPa
		PAO2_local =  fio2*(Pres-PH2O)-((PACO2_local*(1-fio2*(1-RQ)))*RQ**-1) # kPa
	else:
		RQ = 0
		PACO2_local = 0
		PAO2_local = fio2*(Pres-PH2O)
	CcCO2_local = max(CcCO2_local,0.1) # prevent negative values being fed to acidbase2
	CcO2_local = max(CcO2_local,0.1) # prevent negative values being fed to acidbase2
	pH_c = acidbase_2(BE, CcCO2_local, CcO2_local) # pH units
	PcCO2_local = PnCO2_1(CcCO2_local, pH_c, CcO2_local) # kPa
	P50_dxdt = p50_func(pH_c, PcCO2_local, DPG, Temp) # kPa
	PcO2_local = pn_o2_f2(CcO2_local, P50_dxdt) # kPa
	Sats = sn_o2_f1(PcO2_local, pH_c, PcCO2_local, DPG, Temp) # fraction
	K_prime_c = 1.25284e5+3.6917e4*e**(3.8200*Sats) # M**-1.sec**-1
	rateO2 = K_prime_c*alphaO2*60*(1-Sats)*4*HbMol*R*STP_T*STP_P**-1 # ml(O2).ml(bld)**-1.kPa**-1.min**-1
	DLO2 = (DmO2_ivp**-1 + (rateO2*Vc_ivp*1e3)**-1)**-1 # ml(O2).min**-1.kPa**-1
	dO2dt = 100*(Vc_ivp*1e3)**-1*DLO2*(PAO2_local-PcO2_local) # ml(O2).100ml(bld)**-1.min**-1
	DmCO2 = DmO2_ivp*20 # ml(CO2).min**-1.kPa**-1
	DLCO2 = DmCO2 # assumes infinitely fast rate of reaction of CO2
	dCO2dt = 100*(Vc_ivp*1e3)**-1*DLCO2*(PACO2_local-PcCO2_local) # ml.CO2.100mlBlood^-1.min^-1
	return [dO2dt, dCO2dt]
#-------------------------------------------------------------------------------
# Name:		 Transit Time, initial value problem
# Source:	   Wagner and West 1972, altered 
#-------------------------------------------------------------------------------
def dxdt_organ(inputs, t, DmO2_ivp, vol_b, CinputO2, CinputCO2, Q, compound, Qd, CdiO2, CdiCO2):
	CoutputO2 = inputs[0] # these are the only values that are changed with each iteration
	CoutputCO2 = inputs[1]
	solubilitycoefficientO2 = 1
	solubilitycoefficientCO2 = 1
	if (CinputCO2 - CoutputCO2) != (CoutputO2 - CinputO2): 
		CdoO2 = (Qd*CdiO2 - Q*(CoutputO2-CinputO2))/Qd 		# FICK within organ
		CoutputO2 = (Q*CinputO2 + Qd*(CdiO2-CdoO2))/Q 		# FICK in blood 
		CdoCO2 = (Qd*CdiCO2 - Q*(CoutputCO2-CinputCO2))/Qd 	# FICK within organ
		CoutputCO2 = (Q*CinputCO2 + Qd*(CdiCO2-CdoCO2))/Q 	# FICK in blood
	else :# first iteration, catch error
		CdoO2 = CdiO2
		CdoCO2 = CdiCO2
	PorganO2 = ((CdoO2*10*1000**-1)/solubilitycoefficientO2) * STP_P 
	k = STP_P*Temp/STP_T
	PorganCO2 = ((CdoCO2*10*1000**-1)/solubilitycoefficientCO2) * k 
	# equilibrate haemoglobin
	pH_c = acidbase_2(BE, CoutputCO2, CoutputO2) # pH units
	PoutputCO2 = PnCO2_1(CoutputCO2, pH_c, CoutputO2) # kPa
	P50_dxdt = p50_func(pH_c, PoutputCO2, DPG, Temp) # kPa
	PoutputO2 = pn_o2_f2(CoutputO2, P50_dxdt) # kPa
	Sats = sn_o2_f1(PoutputO2, pH_c, PoutputCO2, DPG, Temp) # fraction
	# now do O2
	if DmO2_ivp==0:
		return [0,0] # no diffusion, no change. 
	K_prime_c = 1.25284e5+3.6917e4*e**(3.8200*Sats) # M**-1.sec**-1
	theta = K_prime_c*alphaO2*60*(1-Sats)*4*HbMol*R*STP_T*STP_P**-1 # ml(O2).ml(bld)**-1.kPa**-1.min**-1
	DiffusionO2 = (DmO2_ivp**-1 + (theta*vol_b*1e3)**-1)**-1 # ml(O2).min**-1.kPa**-1
	dO2dt = 100*(vol_b*1e3)**-1*DiffusionO2*(PorganO2-PoutputO2) # ml(O2).100ml(bld)**-1.min**-1
	# now do CO2
	DmCO2 = DmO2_ivp*20 # ml(CO2).min**-1.kPa**-1
	if DmCO2==0:
		return [0,0] # no diffusion, no change. 
	DLCO2 = DmCO2 # assumes infinitely fast rate of reaction of CO2
	dCO2dt = 100*(vol_b*1e3)**-1*DLCO2*(PorganCO2-PoutputCO2) # ml.CO2.100mlBlood^-1.min^-1
	return [dO2dt, dCO2dt]


def runorgan(dtype, preorganCnO2, preorganCnCO2, heter_stat=1e-200, ncomp=20):
	if dtype=="lung":
		# organ settings for organ==LUNG
		organgas = "air"
		organflow = VA # organ flow rate l/min
		Vorganblood = Vc
		Qorganblood = CO*(1-pulm_shunt)
		organO2inputcontent = fio2*(Pres-PH2O) # Content organ input O2 mls_gas/volume
		organCO2inputcontent = 0 # Content organ input CO2 mls_gas/volume
		membranediffusion = DmO2
		organtime = Vorganblood/Qorganblood # in minutes
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		postorganCnO2,postorganCnCO2 = integrate.odeint(dxdt_organ,\
				[preorganCnO2,preorganCnCO2],[0,organtime],\
				(membranediffusion, Vorganblood, preorganCnO2, preorganCnCO2, Qorganblood,\
				organgas, organflow, organO2inputcontent, organCO2inputcontent ),\
				rtol=toleranceoferror_organ, mxstep=100000)[1]
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
		# and get the organ values too:
		CdoO2 = (organflow*organO2inputcontent - Qorganblood*(postorganCnO2-preorganCnO2))/organflow 		# FICK within organ
		CdoCO2 = (organflow*organCO2inputcontent - Qorganblood*(postorganCnCO2-preorganCnCO2))/organflow
		PorganO2 = ((CdoO2*10*1000**-1)) * STP_P*Temp/STP_T 
		PorganCO2 = ((CdoCO2*10*1000**-1)) * STP_P*Temp/STP_T 
		return postorganCnO2, postorganCnCO2, PorganO2, PorganCO2


#-------------------------------------------------------------------------------
#		   mass balance
#-------------------------------------------------------------------------------
def lung_null(CvO2CvCO2):
	VQ = VA*(CO*(1-pulm_shunt))**-1
	[content['CvO2'],content['CvCO2']] = CvO2CvCO2
	[content['CcO2'],content['CcCO2']] = \
	integrate.odeint(dxdt,[content['CvO2'],content['CvCO2']],[0,Vc*(CO*(1-pulm_shunt))**-1],\
	(VQ,DmO2,Vc,content['CvO2'],content['CvCO2']),rtol=toleranceoferror_lung)[1]
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
		if getattr(module, "C%sO2"%i) < 0: setattr(module, "C%sO2"%i, 0.1) # set floor to prevent error in acidbase_2 
		if getattr(module, "C%sCO2"%i) < 0: setattr(module, "C%sCO2"%i, 0.1) 
		setattr(module,"pH_%s"%i, eval("acidbase_2(BE, C%sCO2, C%sO2)"%(i,i)))
		setattr(module,"P%sCO2"%i, eval("PnCO2_1(C%sCO2, pH_%s, C%sO2)"%(i,i,i)))
		setattr(module,"P50_%s"%i, eval("p50_func(pH_%s, P%sCO2, DPG, Temp)"%(i,i)))
		setattr(module,"S%sO2"%i, eval("sn_o2_f2(C%sO2, P50_%s)"%(i,i)))
		setattr(module,"P%sO2"%i, eval("PnO2_1(S%sO2, P50_%s)"%(i,i)))
		setattr(module,"HCO3_%s"%i, eval("vanSlyke_1(pH_%s, BE, S%sO2)"%(i,i)))

#-------------------------------------------------------------------------------
# Name:		 Timecheck
#-------------------------------------------------------------------------------
def fail(message=''):
	print(0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",message)
	sys.exit()

def timecheck(p,t):
	global timelimit
	tc = timeit.default_timer() - t
	if offline_mode or debugging:
		print("position:",p,"time:",tc,"ms")
	if tc > timelimit:
		fail("Killed at position: %s because time=%sms"%(p,tc))

def check_completion():
	global oldvalues
	gdic = globals()
	diff={}
	for key in gdic:
		try: 
			resetdic[key]
			diff[key] = abs(float(gdic[key]) - float(oldvalues[key]))
		except:
			pass
		oldvalues[key] = gdic[key]
	if len(diff.values()) < 2:
		return 100
	return sum(diff.values())

def clearglobals():
	global oldvalues, resetdic
	oldvalues={}
	gdic = globals()
	for key in globals():
		try: 
			resetdic[key]
			float(gdic[key])
		except: 
			continue
		setattr(module, key, 1)

def printglobals():
	print("===================")
	gdic = globals()
	for key in globals():
		try: float(gdic[key])
		except: continue
		print("%s: %s"%(key, gdic[key]))
	print("===================")

def formatoutput(method="full", label=''):
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


#-------------------------------------------------------------------------------
# Name:		 Run model
#-------------------------------------------------------------------------------
def run_all():  # obselete
	tic = timeit.default_timer()
	calculatedconstants()
	timecheck(1,tic) # check time and kill if too long
	populatealvgaseqn()
	timecheck(2,tic) # check time and kill if too long
	updatebloodgascontents()
	timecheck(3,tic) # check time and kill if too long
	updatepartialpressures()
	timecheck(4,tic) # check time and kill if too long
	formatoutput()

def circulate_once(iterationnumber=0):
	global PAO2,PcO2,PaO2,PtO2,PvO2,CcO2,CaO2,CtO2,CvO2,SaO2,PACO2,PcCO2,PaCO2,PtCO2,PvCO2,CcCO2,CaCO2,CtCO2,CvCO2
	global pH_c,pH_a,pH_t,pH_v,HCO3_c,HCO3_a,HCO3_t,HCO3_v,ecmodelivery, trueVO2
	# ====== alveoli and pulmonary capillaries =======
	CcO2, CcCO2, PAO2, PACO2 = runorgan('lung', CvO2, CvCO2)
	updatepartialpressures(compartments = ['c'])
	# ======== arteries ========
	CaO2 = pulm_shunt*CvO2 + (1-pulm_shunt)*CcO2
	CaCO2 = pulm_shunt*CvCO2 + (1-pulm_shunt)*CcCO2 # mlO2/100mlblood
	updatepartialpressures(compartments = ['a'])
	# ========= tissues ==========
	Qtissue = CO*(1-tissue_shunt)
	CtO2 = (Qtissue*CaO2*10 - trueVO2*1000)/ (Qtissue*10)
	CtCO2 = (Qtissue*CaCO2*10 + trueVO2*RQ*1000)/ (Qtissue*10)
	updatepartialpressures(compartments = ['t'])
	# ========= veins ==========
	CvO2 = (CtO2*Qtissue + CaO2*CO*tissue_shunt)/CO
	CvCO2 = (CtCO2*Qtissue + CaCO2*CO*tissue_shunt)/CO
	updatepartialpressures(compartments = ['v'])
	# ========= set VO2 ==========
	criticalOER = 0.94 # unrealistic maximum
	criticalDO2 = (CO*CaO2*10)*criticalOER/1000 # mlsO2/min
	if criticalDO2 < VO2: #then set trueVO2 = criticalDO2, but prevent oscillation:
		trueVO2 = trueVO2 - float(trueVO2-criticalDO2)/(iterationnumber/VO2correctionspeed+1) # only go a hundredth of the way to avoid oscillation
		#print("CO %s CaO2 %s criticalDO2 %s < VO2 %s so correcting downwards to %s"%(CO, CaO2, criticalDO2, VO2, trueVO2)
	elif trueVO2<VO2: # then it needs to come back up
		trueVO2 = VO2

#========================
#========================
#========================

resetlist=[
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
	'ecmodelivery',
	'trueVO2',
	'P50_a',
	'P50_v',
	]
resetdic = {x:1 for x in resetlist}
oldvalues={}
# ------  ------  ------  ------  ------  ------  ------  ------ 
CcCO2 = 1
CcO2 = 1
CvO2 = 1
CvCO2 = 1
CaCO2 = 1
CaO2 = 1
CtCO2 = 1
CtO2 = 1
ecmodelivery = 0
# ------  ------  ------  ------  ------  ------  ------  ------ 
# getinputs()
calculatedconstants()
populatealvgaseqn()
if __name__ == "__main__":
	try:
		print("pre-updatebloodgascontents")
		updatebloodgascontents()
		print("pre-updatepartialpressures")
		updatepartialpressures()
	except:
		print(content)
		pass

	for i in range(maxruns):
		circulate_once(i)
		globaldiff = check_completion()
		#print(i, globaldiff)
		if globaldiff < globaldifftolerance:
			break
	formatoutput()

# ------  ------  ------  ------  ------  ------  ------  ------ 
