import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from math import *

# Stanard global variables
std_tp = 37
std_ph = 7.4
std_paco2 = 40
std_pao2 = 84

def adjust_pa02(pao2, tp, ph, paco2):
    return pao2 * 10**(0.024 * (std_tp - tp) + 0.40*(ph - std_ph) + 0.06*(log10(std_paco2) - log10(paco2)))

def hb_o2_sat_pct(pao2):
    a1 = -8.5322289 * 10**3
    a2 = 2.1214010 * 10**3
    a3 = -6.7073989 * 10**1
    a4 = 9.3596087 * 10**5
    a5 = -3.1346258 * 10**4
    a6 = 2.3961674 * 10**3
    a7 = -6.7104406 * 10**1    
    return 100*(a1*pao2 + a2*pao2**2 + a3*pao2**3 + pao2**4) / (a4 + a5*pao2 + a6*pao2**2 + a7*pao2**3 + pao2**4)

def hb_o2_diss_curve(pao2, tp, ph, paco2):
    return hb_o2_sat_pct(adjust_pa02(pao2, tp, ph, paco2))

def hb_o2_sat_pct_p50(pao2, tp, ph, paco2):
    return hb_o2_diss_curve(pao2, tp, ph, paco2) - 50

def find_p50_pao2(tp, ph, paco2):
    return optimize.fsolve(hb_o2_sat_pct_p50, 26, args=(tp, ph, paco2))

std_spo2 = hb_o2_diss_curve(std_pao2, std_tp, std_ph, std_paco2)
std_p50_pao2 = find_p50_pao2(std_tp, std_ph, std_paco2)
std_p50_spo2 = hb_o2_diss_curve(std_p50_pao2, std_tp, std_ph, std_paco2)


max_pao2 = 110
x_pao2 = np.linspace(0, max_pao2, max_pao2*2)
y_spo2_std = hb_o2_diss_curve(x_pao2, std_tp, std_ph, std_paco2)
y_spo2_tp = hb_o2_diss_curve(x_pao2, std_tp-2, std_ph, std_paco2)
plt.plot(std_pao2, std_spo2, 'o', markersize=4, label='Standard PaO2')
plt.plot(std_p50_pao2, std_p50_spo2, 'o', markersize=4, label='Standard P50')
plt.plot(x_pao2, y_spo2_std, label='Standard Curve')
plt.plot(x_pao2, y_spo2_tp, label='Decreased Temperature Curve')
plt.grid(True)
plt.xlabel("PaO2 mmHg")
plt.xlim(left=0, right=110)
plt.xticks(np.arange(0,110,10))
plt.ylabel("SaO2 %")
plt.ylim(bottom=-10, top=100)
plt.yticks(np.arange(0,110,10))
plt.legend(loc='lower right')
plt.show()


