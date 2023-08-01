#----------------------------------------------------------------------------
# Imported modules
#----------------------------------------------------------------------------
import timeit
from scipy import optimize, integrate, special
from math import exp, e, log10, log

# ###########################################################################
# Global Constants
# ###########################################################################
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
STD_UNIT_CONV_VALS = {
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

# ###########################################################################
# Global functions
# ###########################################################################
def get_std_input_parms_dict():
    return {
        'RR':[12.0,'bpm'],
        'VT':[0.475,'l'],
        'VD':[0.11,'l'],
        'fio2':[0.21,'fraction'],
        'alt':[0,'m'],
        'CO':[6.5,'l/min'],
        'Vc':[0.075,'l'],
        'VO2':[0.25,'l/min'],
        'Temp':[309.65,'K'],
    }

#----------------------------------------------------------------------------
# Name:		 Unit Conversion
# Purpose:	Convert to SI units
#----------------------------------------------------------------------------
def get_standard_unit_value(value, unit):
    if unit == 'C':
        return value+273.15
    elif unit == 'F':
        return (5*(value-32)*9**-1)+273.15
    else:
        return value*STD_UNIT_CONV_VALS[unit]
    
#----------------------------------------------------------------------------
# Name:			alpha_o2_func
# Source:		Dash & Bassingthwaight 2010
#----------------------------------------------------------------------------
def alpha_o2_func(temp, wpl):
    return 0.1333**-1*(1.37-0.0137*(temp-310.15)+0.00058*(temp-310.15)**2)*(1e-6*wpl**-1) # M/kPa

#----------------------------------------------------------------------------
# Name:		   alpha_co2_func
# Source:	   Kelman 1967
#----------------------------------------------------------------------------
def alpha_co2_func(temp, wpl):
    return 0.1333**-1*(3.07-0.057*(temp-310.15)+0.002*(temp-310.15)**2)*(1e-5*wpl**-1) # M/kPa



# ###########################################################################
# Class Def
# ###########################################################################
class OxygenCascade:

    #----------------------------------------------------------------------------
    # Class level constants
    #----------------------------------------------------------------------------
    MAX_RUNS=2
    TOL_ERR_SO2 =	1e-7
    TOL_ERR_LUNG =	0.00001
    TOL_ERR_ORGAN = 0.00001
    TOL_ERR_PH =	0.00001
    VO2_CORRECTION_SPEED = 3 # normally 3; faster has higher chance of oscillation
    GLOBAL_DIFF_TOL = 0.045

    def __init__(self, input_parms_dict) -> None:
        self.AGE={}
        self.content = {}

        # Init instance vars with standard values
        self.pulm_shunt=get_standard_unit_value(float(0.03),'fraction')
        self.DmO2=get_standard_unit_value(float(300),'mlO2/min/kPa')
        self.Hb=get_standard_unit_value(float(150),'g/l')
        self.BE=get_standard_unit_value(float(0),'mEq/l')
        self.DPG=get_standard_unit_value(float(4.65),'mmol/l')
        self.MCHC=get_standard_unit_value(float(340),'g/l')
        self.tissue_shunt=get_standard_unit_value(float(0.05),'fraction')
        self.RQ=get_standard_unit_value(float(0.8),'fraction')
        self.ecmosites=[]
        self.Qecmo=get_standard_unit_value(0,'l/min')
        self.hetindex=0

        # Init user input instance vars with standard values
        self.RR=get_standard_unit_value(float(12.0),'bpm')
        self.VT=get_standard_unit_value(float(0.475),'l')
        self.VD=get_standard_unit_value(float(0.11),'l')
        self.fio2=get_standard_unit_value(float(0.21),'fraction')
        self.alt=get_standard_unit_value(float(0),'m')
        self.CO=get_standard_unit_value(float(6.5),'l/min')
        self.Vc=get_standard_unit_value(float(0.075),'l')
        self.VO2=get_standard_unit_value(float(0.25),'l/min')
        self.Temp=get_standard_unit_value(float(309.65),'K')

        # Override user input instance vars with user input values
        self.input_parms = input_parms_dict
        for key in input_parms_dict:
            setattr(self, key, get_standard_unit_value(float(input_parms_dict[key][0]),input_parms_dict[key][1]))

        #------------------------------------------------------------------------
        # Calculated instance variables based on user input data
        #------------------------------------------------------------------------
        #------------------------------------------------------------------------
        # Name:		 Atmospheric pressure (kPa) from altitude (m)
        # Source:	West 1996
        #------------------------------------------------------------------------
        self.Pres = exp(6.63268-0.1112*(self.alt*1e-3)-0.00149*(self.alt*1e-3)**2)*0.1333 # kPa
        #------------------------------------------------------------------------
        # Name:		 Saturated vapour pressure of water (kPa)
        # Source:	   Wagner & Pruss 1993
        #------------------------------------------------------------------------
        tau = 1 - self.Temp*TEMP_CRITICAL**-1 # fraction
        self.PH2O = PRES_CRITICAL*e**(TEMP_CRITICAL*self.Temp**-1*(A1*tau+A2*tau**1.5+A3*tau**3+A4*tau**3.5+A5*tau**4+A6*tau**7)) # kPa
        #------------------------------------------------------------------------
        # Name:		 PIO2
        #------------------------------------------------------------------------
        self.PIO2 = self.fio2*(self.Pres-self.PH2O) # kPa
        #------------------------------------------------------------------------
        # Name:		 [Hb]
        #------------------------------------------------------------------------
        self.HbMol = self.Hb*64458**-1 # M
        self.Hct  = float(self.Hb)/self.MCHC
        #------------------------------------------------------------------------
        # Name:		 Fractional water space of blood
        #------------------------------------------------------------------------
        self.Wbl = (1-self.Hct)*WC_PLASMA+self.Hct*WC_RBC # fraction
        #------------------------------------------------------------------------
        # Name:		 Alveolar Ventilation
        #------------------------------------------------------------------------
        self.VA = self.RR*(self.VT-self.VD) # l/min
        self.VQ = self.VA*(self.CO*(1-self.pulm_shunt))**-1
        #------------------------------------------------------------------------
        # Name:		 alpha O2 / CO2
        # Source:	   Dash & Bassingthwaight 2010
        #------------------------------------------------------------------------
        self.alphaO2 = alpha_o2_func(self.Temp, WC_PLASMA) # M/kPa
        self.alphaCO2 = alpha_co2_func(self.Temp, WC_PLASMA) # M/kPa
        #------------------------------------------------------------------------
        #----------- AND FIX INPUT UNITS
        self.BE = self.BE/1000 # input is mEq/l
        self.trueVO2 = self.VO2 # user-defined VO2 is an aspiration...
        #------------------------------------------------------------------------
        # Alveolar ventilation equation. Simplifies to PACO2=0.863*VCO2/VA under normal conditions
        self.PACO2 = self.trueVO2*self.RQ*STD_PRES*self.Temp*STD_TEMP**-1*self.VA**-1 
        self.PAO2 = self.PIO2 - ((self.PACO2*(1-self.fio2*(1-self.RQ)))*self.RQ**-1)# Alveolar gas Equation




    #-###########################################################################
    # Class level methods
    #-###########################################################################
    #----------------------------------------------------------------------------
    # Name:		   p50_func
    # Source:	   Dash & Bassingthwaight 2010
    #----------------------------------------------------------------------------
    def p50_func(self, ph, pn_co2, dpg, temp):
        p50_ph = -2.84*(ph+log10(0.69)-7.24)+1.18*(ph+log10(0.69)-7.24)**2
        p50_pn_co2 = 4.82e-2*(pn_co2-5.332)+3.64e-5*(pn_co2-5.332)**2
        p50_dpg = 1.06e2*(dpg-4.65e-3)-2.62e3*(dpg-4.65e-3)**2
        p50_temp = 1.99e-1*(temp-310.15)+5.78e-3*(temp-310.15)**2+9.33e-05*(temp-310.15)**3
        return 3.57+p50_pn_co2+p50_ph+p50_dpg+p50_temp # kPa
    
    #----------------------------------------------------------------------------
    # Name:		 SnO2
    # Source:	   Dash & Bassingthwaight 2010
    #----------------------------------------------------------------------------
    def sn_o2_f1(self, pn_o2, ph, pn_co2, dpg, temp): # Equation B.3
        '''Calculate O2 SATURATION from PARTIAL PRESSURE''' # fraction
        try:
            pno2p50 = (pn_o2*self.p50_func(ph, pn_co2, dpg, temp)**-1)**(1+N0)
            return (pno2p50)/(1+(pno2p50))
        except BaseException as e:
            print(pn_o2, ph, pn_co2, dpg, temp)
            print("failed at sn_o2_f1")
            raise e
    def sn_o2_f2_null(self, sats, cn_o2, p50_sn_o2):
        '''returns 0 '''
        return self.Wbl*self.alphaO2*p50_sn_o2*(sats*(1-sats)**-1)**((1+N0)**-1) + \
            (4*self.HbMol)*sats - (cn_o2*(STD_R*STD_TEMP*STD_PRES**-1*1e2)**-1)
    def sn_o2_f2(self, cn_o2, p50_sn_o2):
        '''Calculate O2 SATURATION from CONTENT''' # fraction
        cn_o2 = max(cn_o2,0.1) # prevent negative values being fed to acidbase2
        p50_sn_o2 = max(p50_sn_o2,0.1) # prevent negative values being fed to acidbase2
        return optimize.brentq(self.sn_o2_f2_null,1e-15,(1-1e-15),args=(cn_o2, p50_sn_o2),rtol=self.TOL_ERR_SO2)
    #----------------------------------------------------------------------------
    # Name:		 Blood O2 content
    # Source:	   Dash & Bassingthwaight 2010
    #----------------------------------------------------------------------------
    def cn_o2_f1(self, pn_o2, ph, pn_co2, dpg, temp):
        '''Calculate O2 CONTENT from PARTIAL PRESSURE'''
        return (self.Wbl*self.p50_func(ph,pn_co2,dpg,temp)*(self.sn_o2_f1(pn_o2,ph,pn_co2,dpg,temp)*(1\
        -self.sn_o2_f1(pn_o2,ph,pn_co2,dpg,temp))**-1)**((1+N0)**-1)*self.alphaO2\
        +4*self.HbMol*self.sn_o2_f1(pn_o2,ph,pn_co2,dpg,temp))*(STD_R*STD_TEMP*STD_PRES**-1*1e2)
        # ml of O2 per 100ml blood STP
    #----------------------------------------------------------------------------
    # Name:		 Blood O2 partial pressure
    # Source:	   Dash & Bassingthwaight 2010
    #----------------------------------------------------------------------------
    def pn_o2_f1(self, sn_o2, p50_sn_o2):
        '''Calculate O2 PARTIAL PRESSURE from SATURATION'''
        return p50_sn_o2*(sn_o2*(1-sn_o2)**-1)**((1+N0)**-1) # kPa
    def pn_o2_f2(self, cn_o2, p50_sn_o2):
        '''Calculate O2 PARTIAL PRESSURE from CONTENT'''
        sn_o2 = self.sn_o2_f2(cn_o2, p50_sn_o2)
        return p50_sn_o2*(sn_o2*(1-sn_o2)**-1)**((1+N0)**-1) # kPa
    #----------------------------------------------------------------------------
    # Name:		 CnCO2
    # Source:	   Dash & Bassingthwaight 2010
    #----------------------------------------------------------------------------
    def k_ratio_func(self, pn_co2, ph, temp, wpl):
        hrbc = 10**-(ph+log(0.69,10)) # M
        CO2 = pn_co2*alpha_co2_func(temp, wpl) # M
        return (K_PRIME_2*CO2*(1+K_2_PRIME_2*hrbc**-1)
            +(1+hrbc*K_2_PRIME_5**-1))*(K_PRIME_3*CO2*(1+K_2_PRIME_3*hrbc**-1)
            +(1+hrbc*K_2_PRIME_6**-1))**-1 # unitless
    def cn_co2_dissolved(self, pn_co2):
        '''Calculate CO2 dissolved in WHOLE BLOOD'''
        return self.Wbl*self.alphaCO2*pn_co2 # M
    def cn_co2_bicarb(self, ph, pn_co2):
        '''Calculate CO2 as bicarbonate in WHOLE BLOOD'''
        return ((1-self.Hct)*WC_PLASMA+self.Hct*WC_RBC*0.69)*self.hh_f1(pn_co2, ph) # M
    def cn_co2_hb_bound(self, pn_co2, pn_o2, ph):
        '''Calculate Hb-CO2 from PARTIAL PRESSURE'''
        k_prime_4 = (self.alphaO2*pn_o2)**N0*self.k_ratio_func(pn_co2,ph,self.Temp,WC_PLASMA)*(self.p50_func(ph, pn_co2, self.DPG, self.Temp)*self.alphaO2)**-(1+N0)
        s_hb_co2 = ( \
        K_PRIME_2*self.alphaCO2*pn_co2*(1+K_2_PRIME_2*(10**-(ph+log(0.69,10)))**-1) + \
        K_PRIME_3*self.alphaCO2*pn_co2*(1+K_2_PRIME_3*(10**-(ph+log(0.69,10)))**-1)*k_prime_4*self.alphaO2*pn_o2
        )*( \
        K_PRIME_2*self.alphaCO2*pn_co2*(1+K_2_PRIME_2*(10**-(ph+log(0.69,10)))**-1) + \
        (1+(10**-(ph+log(0.69,10)))*K_2_PRIME_5**-1) + \
        k_prime_4*self.alphaO2*pn_o2*(K_PRIME_3*self.alphaCO2*pn_co2*(1+K_2_PRIME_3*(10**-(ph+log(0.69,10)))**-1) + \
        (1+(10**-(ph+log(0.69,10)))*K_2_PRIME_6**-1)) \
        )**-1
        return 4*self.HbMol*s_hb_co2 # M
    def cn_co2_f1(self, ph, pn_co2, pn_o2):
        '''Calculate CO2 CONTENT from PARTIAL PRESSURE'''
        return (self.cn_co2_hb_bound(pn_co2,pn_o2,ph)+self.cn_co2_bicarb(ph,pn_co2)\
            +self.cn_co2_dissolved(pn_co2))*(STD_R*STD_TEMP*STD_PRES**-1*1e2) # ml CO2 per 100ml blood STP
    def pn_co2_f1_null(self, pn_co2, cn_co2, ph, cn_o2):
        '''returns 0'''
        p50_pn_co2_null = self.p50_func(ph, pn_co2, self.DPG, self.Temp)
        return self.cn_co2_f1(ph, pn_co2, self.pn_o2_f2(cn_o2, p50_pn_co2_null))-cn_co2 # null
    def pn_co2_f1(self, cn_co2, ph, cn_o2):
        '''Calculate CO2 PARTIAL PRESSURE from CONTENTS'''
        args = cn_co2,ph,cn_o2
        pn_co2 = optimize.newton(self.pn_co2_f1_null,self.AGE['PcCO2'],args=args,tol=0.01)
        return pn_co2 # kPa
    #----------------------------------------------------------------------------
    # Name:	 Henderson-Hasselbalch Equation
    #----------------------------------------------------------------------------
    def hh_f1(self, pn_co2, ph):
        '''Calculate CO2 as bicarbonate in SOLUTION'''
        return (K_1*self.alphaCO2*pn_co2)*(10**-ph)**-1 # M
    #----------------------------------------------------------------------------
    # Name:	 van Slyke Equation
    # Source:   Siggaard-Andersen 1977
    #----------------------------------------------------------------------------
    def van_slyke_f1(self, ph, base_excess, sn_o2):
        '''returns bicarbonate from PLASMA pH @ 37 deg C and BE'''
        zeta = 1-(0.0143*self.Hb*1e-1)
        beta = 9.5+1.63*self.Hb*1e-1
        return ((base_excess*1e3 - 0.2*self.Hb*1e-1*(1-sn_o2))*zeta**-1 - beta*(ph-7.4) + STD_BICARB*1e3)*1e-3 # M
    #----------------------------------------------------------------------------
    # Name:	 simultaneous solution Henderson-Hasselbalch and van Slyke
    #----------------------------------------------------------------------------
    def acidbase_f1_null(self, ph,base_excess, pn_co2, pn_o2):
        '''returns 0'''
        sn_o2 = self.sn_o2_f1(pn_o2, ph, pn_co2, self.DPG, self.Temp)
        zeta = 1-(0.0143*self.Hb*1e-1)
        beta = 9.5+1.63*self.Hb*1e-1
        return ((base_excess*1e3 - 0.2*self.Hb*1e-1*(1-sn_o2))*zeta**-1 - beta*(ph-7.4) + STD_BICARB*1e3)*1e-3 - self.hh_f1(pn_co2, ph) # null

    def acidbase_f1(self, base_excess, pn_co2, pn_o2):
        '''returns PLASMA pH from PARTIAL PRESSURE '''
        return optimize.brentq(self.acidbase_f1_null,1,14,args=(base_excess,pn_co2,pn_o2),rtol=self.TOL_ERR_PH) # pH units

    def acidbase_f2_null(self, ph, base_excess, cn_co2, cn_o2):
        '''returns 0'''
        pn_co2 = self.pn_co2_f1(cn_co2, ph, cn_o2)
        sn_o2 = self.sn_o2_f2(cn_o2, self.p50_func(ph, pn_co2, self.DPG, self.Temp))
        zeta = 1-(0.0143*self.Hb*1e-1)
        beta = 9.5+1.63*self.Hb*1e-1
        return ((base_excess*1e3 - 0.2*self.Hb*1e-1*(1-sn_o2))*zeta**-1 - beta*(ph-7.4) + STD_BICARB*1e3)*1e-3 - self.hh_f1(pn_co2, ph) # null

    def acidbase_f2(self, base_excess, cn_co2, cn_o2):
        '''returns PLASMA pH from CONTENT '''
        # prevent negative values being fed to acidbase_f2
        cn_co2 = max(cn_co2,0.1) 
        cn_o2 = max(cn_o2,0.1)
        return optimize.brentq(self.acidbase_f2_null,1,14,args=(base_excess, cn_co2, cn_o2),rtol=self.TOL_ERR_PH) # pH units

    #-###########################################################################
    # Name:		 COMPARTMENT SPECIFIC FUNCTIONS
    #-###########################################################################
    #-------------------------------------------------------------------------------
    #		   Alveolar Gas Equation
    #-------------------------------------------------------------------------------
    def populate_alv_gas_eqn(self):
        self.AGE['PACO2'] = self.VO2*self.RQ*STD_PRES*self.Temp*STD_TEMP**-1*self.VA**-1 # Alveolar ventilation equation
        self.AGE['PcCO2'] = self.AGE['PACO2'] # assumes complete equilibrium
        self.AGE['PAO2'] = self.PIO2 - ((self.AGE['PACO2']*(1-self.fio2*(1-self.RQ)))*self.RQ**-1)# Alveolar gas Equation
        self.AGE['PcO2'] = self.AGE['PAO2'] # assumes complete equilibrium
        self.AGE['pH'] = self.acidbase_f1(self.BE,self.AGE['PcCO2'],self.AGE['PcO2'])
        self.AGE['CcCO2'] = self.cn_co2_f1(self.AGE['pH'], self.AGE['PcCO2'], self.AGE['PcO2'])
        self.AGE['CcO2'] = self.cn_o2_f1(self.AGE['PcO2'], self.AGE['pH'], self.AGE['PcCO2'], self.DPG, self.Temp)
        self.AGE['CvO2'] = self.AGE['CcO2'] - 100*(self.VO2*(self.CO*(1-self.pulm_shunt))**-1) 
        self.AGE['CvCO2'] = self.AGE['CcCO2'] + 100*(self.VO2*self.RQ)*(self.CO*(1-self.pulm_shunt))**-1
        self.PAO2 = self.AGE['PAO2']
        self.PACO2 = self.AGE['PACO2']
        #print(AGE)
    #-------------------------------------------------------------------------------
    # Name:		 Transit Time, initial value problem
    # Source:	   Wagner and West 1972
    #-------------------------------------------------------------------------------
    def dxdt(self, inputs, t, vq, dm_o2_ivp, vc_ivp, cv_o2, cv_co2):
        cc_o2 = inputs[0]
        cc_co2 = inputs[1]
        pa_co2 = 0
        pa_o2 = self.fio2*(self.Pres-self.PH2O)
        if (cv_co2 - cc_co2) != (cc_o2 - cv_o2):
            rq = (cv_co2-cc_co2)*(cc_o2-cv_o2)**-1 # ratio
            pa_co2 = vq**-1*(cc_o2-cv_o2)*1e-2*rq*STD_PRES*self.Temp*STD_TEMP**-1 # kPa
            pa_o2 =  self.fio2*(self.Pres-self.PH2O)-((pa_co2*(1-self.fio2*(1-rq)))*rq**-1) # kPa
        ph_c = self.acidbase_f2(self.BE, cc_co2, cc_o2) # pH units
        pc_co2 = self.pn_co2_f1(cc_co2, ph_c, cc_o2) # kPa
        p50_dxdt = self.p50_func(ph_c, pc_co2, self.DPG, self.Temp) # kPa
        pc_o2 = self.pn_o2_f2(cc_o2, p50_dxdt) # kPa
        sats = self.sn_o2_f1(pc_o2, ph_c, pc_co2, self.DPG, self.Temp) # fraction
        k_prime_c = 1.25284e5+3.6917e4*e**(3.8200*sats) # M**-1.sec**-1
        rate_o2 = k_prime_c*self.alphaO2*60*(1-sats)*4*self.HbMol*STD_R*STD_TEMP*STD_PRES**-1 # ml(O2).ml(bld)**-1.kPa**-1.min**-1
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
    def dxdt_organ(self, inputs, t, dm_o2_ivp, vol_b, c_in_o2, c_in_co2, q, compound, qd, c_di_o2, c_di_co2):
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
        k = STD_PRES*self.Temp/STD_TEMP
        p_organ_co2 = (c_do_co2*10*1000**-1) * k 
        # equilibrate haemoglobin
        ph_c = self.acidbase_f2(self.BE, c_out_co2, c_out_o2) # pH units
        p_out_co2 = self.pn_co2_f1(c_out_co2, ph_c, c_out_o2) # kPa
        p50_dxdt = self.p50_func(ph_c, p_out_co2, self.DPG, self.Temp) # kPa
        p_out_o2 = self.pn_o2_f2(c_out_o2, p50_dxdt) # kPa
        sats = self.sn_o2_f1(p_out_o2, ph_c, p_out_co2, self.DPG, self.Temp) # fraction
        # now do O2
        k_prime_c = 1.25284e5+3.6917e4*e**(3.8200*sats) # M**-1.sec**-1
        theta = k_prime_c*self.alphaO2*60*(1-sats)*4*self.HbMol*STD_R*STD_TEMP*STD_PRES**-1 # ml(O2).ml(bld)**-1.kPa**-1.min**-1
        diffusion_o2 = (dm_o2_ivp**-1 + (theta*vol_b*1e3)**-1)**-1 # ml(O2).min**-1.kPa**-1
        do2_dt = 100*(vol_b*1e3)**-1*diffusion_o2*(p_organ_o2-p_out_o2) # ml(O2).100ml(bld)**-1.min**-1
        # now do CO2
        dm_co2 = dm_o2_ivp*20 # ml(CO2).min**-1.kPa**-1
        if dm_co2==0:
            return [0,0] # no diffusion, no change. 
        dl_co2 = dm_co2 # assumes infinitely fast rate of reaction of CO2
        dco2_dt = 100*(vol_b*1e3)**-1*dl_co2*(p_organ_co2-p_out_co2) # ml.CO2.100mlBlood^-1.min^-1
        return [do2_dt, dco2_dt]

    def run_organ(self, dtype, preorgan_cn_o2, preorgan_cn_co2):
        if dtype=="lung":
            # organ settings for organ==LUNG
            organgas = "air"
            organflow = self.VA # organ flow rate l/min
            v_organ_blood = self.Vc
            q_organ_blood = self.CO*(1-self.pulm_shunt)
            organ_o2_in_content = self.fio2*(self.Pres-self.PH2O) # Content organ input O2 mls_gas/volume
            organ_co2_in_content = 0 # Content organ input CO2 mls_gas/volume
            membrane_diffusion = self.DmO2
            organtime = v_organ_blood/q_organ_blood # in minutes
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            post_organ_cn_o2, post_organ_cn_co2 = integrate.odeint(self.dxdt_organ,\
                    [preorgan_cn_o2,preorgan_cn_co2],[0,organtime],\
                    (membrane_diffusion, v_organ_blood, preorgan_cn_o2, preorgan_cn_co2, q_organ_blood,\
                    organgas, organflow, organ_o2_in_content, organ_co2_in_content ),\
                    rtol=self.TOL_ERR_ORGAN, mxstep=1000)[1]
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            # and get the organ values too:
            c_do_o2 = (organflow*organ_o2_in_content - q_organ_blood*(post_organ_cn_o2-preorgan_cn_o2))/organflow 		# FICK within organ
            c_do_co2 = (organflow*organ_co2_in_content - q_organ_blood*(post_organ_cn_co2-preorgan_cn_co2))/organflow
            post_organ_o2 = (c_do_o2*10*1000**-1) * STD_PRES*self.Temp/STD_TEMP 
            post_organ_co2 = (c_do_co2*10*1000**-1) * STD_PRES*self.Temp/STD_TEMP 
            return post_organ_cn_o2, post_organ_cn_co2, post_organ_o2, post_organ_co2

    #-------------------------------------------------------------------------------
    #		   mass balance
    #-------------------------------------------------------------------------------
    def lung_null(self, cv_o2__cv_co2):
        VQ = self.VA*(self.CO*(1-self.pulm_shunt))**-1
        [self.content['CvO2'], self.content['CvCO2']] = cv_o2__cv_co2
        [self.content['CcO2'], self.content['CcCO2']] = \
        integrate.odeint(self.dxdt,[self.content['CvO2'], self.content['CvCO2']], [0,self.Vc*(self.CO*(1-self.pulm_shunt))**-1],\
        (VQ, self.DmO2, self.Vc, self.content['CvO2'], self.content['CvCO2']), rtol=self.TOL_ERR_LUNG)[1]
        out = [self.content['CcO2'] - self.content['CvO2'] - 100*((self.trueVO2)*(self.CO*(1-self.pulm_shunt))**-1)]
        out.append(self.content['CvCO2'] - self.content['CcCO2'] - 100*((self.RQ*self.trueVO2)*(self.CO*(1-self.pulm_shunt))**-1))
        return out
    #-------------------------------------------------------------------------------
    #		 compartment contents
    #-------------------------------------------------------------------------------
    def updatebloodgascontents(self):
        #-------------------------------------------------------------------------------
        # Pulmonary capillaries & Veins
        #-------------------------------------------------------------------------------
        optimize.fsolve(self.lung_null, [self.AGE['CvO2'], self.AGE['CvCO2']], xtol=0.001)

        #-------------------------------------------------------------------------------
        # Arteries
        #-------------------------------------------------------------------------------
        self.content['CaCO2'] = self.pulm_shunt*self.content['CvCO2'] + (1-self.pulm_shunt)*self.content['CcCO2'] # mlO2/100mlblood
        self.content['CaO2'] = self.pulm_shunt*self.content['CvO2'] + (1-self.pulm_shunt)*self.content['CcO2'] # mlO2/100mlblood
        #-------------------------------------------------------------------------------
        # Tissues
        #-------------------------------------------------------------------------------
        QO2 = self.CO*(1-self.tissue_shunt)*self.content['CaO2']*100**-1 # lO2/min
        self.content['CtO2'] = 100*(QO2-self.VO2)/((1-self.tissue_shunt)*self.CO) # mlO2/100mlblood
        self.content['CtCO2'] = self.content['CaCO2'] + 100*self.VO2*self.RQ*(self.CO*(1-self.tissue_shunt))**-1 # ml.CO2/100ml.Blood
        #-------------------------------------------------------------------------------
        # copy to globals
        #-------------------------------------------------------------------------------
        for name, value in self.content.items():
            setattr(self, name, value)

    #-------------------------------------------------------------------------------
    #		   compartment partial pressures, acid base and O2 saturation
    #-------------------------------------------------------------------------------
    def updatepartialpressures(self, compartments = ['c','a','t','v']):
        for i in compartments:
            setattr(self,"pH_%s"%i, eval("self.acidbase_f2(self.BE, self.C%sCO2, self.C%sO2)"%(i,i)))
            setattr(self,"P%sCO2"%i, eval("self.pn_co2_f1(self.C%sCO2, self.pH_%s, self.C%sO2)"%(i,i,i)))
            setattr(self,"P50_%s"%i, eval("self.p50_func(self.pH_%s, self.P%sCO2, self.DPG, self.Temp)"%(i,i)))
            setattr(self,"S%sO2"%i, eval("self.sn_o2_f2(self.C%sO2, self.P50_%s)"%(i,i)))
            setattr(self,"P%sO2"%i, eval("self.pn_o2_f1(self.S%sO2, self.P50_%s)"%(i,i)))
            setattr(self,"HCO3_%s"%i, eval("self.van_slyke_f1(self.pH_%s, self.BE, self.S%sO2)"%(i,i)))

    #-------------------------------------------------------------------------------
    # Name:		 Timecheck
    #-------------------------------------------------------------------------------
    def fail(self, message=''):
        print(0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",0,"||",message)

    def format_output_json(self):
        output_data = {
                'PatmosO2':round(self.fio2*self.Pres,4),
                'PIO2':round(self.fio2*(self.Pres - self.PH2O),4),
                'PAO2':round(self.PAO2,4),
                'PcO2':round(self.PcO2,4),
                'PtO2':round(self.PtO2,4),
                'PvO2':round(self.PvO2,4),
                'PaO2':round(self.PaO2,4),
                'PaCO2':round(self.PaCO2,4),
                'pH':round(self.pH_a,2),
                'H+':round(10**9*10**-self.pH_a,2),
                'HCO3':round(self.HCO3_a*1000,2),
                'SaO2':round(self.SaO2*100,2),
                }
        return output_data

    #-------------------------------------------------------------------------------
    # Name:	Circulate in the lungs one time
    #-------------------------------------------------------------------------------
    def circulate_once(self, iterationnumber=0):
        # ====== alveoli and pulmonary capillaries =======
        self.CcO2, self.CcCO2, self.PAO2, self.PACO2 = self.run_organ('lung', self.CvO2, self.CvCO2)
        self.updatepartialpressures(compartments = ['c'])
        # ======== arteries ========
        self.CaO2 = self.pulm_shunt*self.CvO2 + (1-self.pulm_shunt)*self.CcO2
        self.CaCO2 = self.pulm_shunt*self.CvCO2 + (1-self.pulm_shunt)*self.CcCO2 # mlO2/100mlblood
        self.updatepartialpressures(compartments = ['a'])
        # ========= tissues ==========
        q_tissue = self.CO*(1-self.tissue_shunt)
        self.CtO2 = (q_tissue*self.CaO2*10 - self.trueVO2*1000)/ (q_tissue*10)
        self.CtCO2 = (q_tissue*self.CaCO2*10 + self.trueVO2*self.RQ*1000)/ (q_tissue*10)
        self.updatepartialpressures(compartments = ['t'])
        # ========= veins ==========
        self.CvO2 = (self.CtO2*q_tissue + self.CaO2*self.CO*self.tissue_shunt)/self.CO
        self.CvCO2 = (self.CtCO2*q_tissue + self.CaCO2*self.CO*self.tissue_shunt)/self.CO
        self.updatepartialpressures(compartments = ['v'])
        # ========= set VO2 ==========
        critical_oer = 0.94 # unrealistic maximum
        critical_do2 = (self.CO*self.CaO2*10)*critical_oer/1000 # mlsO2/min
        if critical_do2 < self.VO2: #then set trueVO2 = criticalDO2, but prevent oscillation:
            self.trueVO2 = self.trueVO2 - float(self.trueVO2-critical_do2)/(iterationnumber/self.VO2_CORRECTION_SPEED+1) # only go a hundredth of the way to avoid oscillation
        elif self.trueVO2 < self.VO2: # then it needs to come back up
            self.trueVO2 = self.VO2

    #-------------------------------------------------------------------------------
    # Name:	Run full model
    #-------------------------------------------------------------------------------
    def run_model(self):
        start_time = timeit.default_timer()
        print(start_time)
        try:
            #print("pre-populate_alv_gas_eqn")
            self.populate_alv_gas_eqn()
            #print("pre-updatebloodgascontents")
            self.updatebloodgascontents()
            #print("pre-updatepartialpressures")
            self.updatepartialpressures()
            for i in range(self.MAX_RUNS):
                self.circulate_once(i)

        except BaseException as e:
            print(self.content)
            print(e)
            raise e
            
        print(timeit.default_timer() - start_time)
        return self.format_output_json()



    #========================
