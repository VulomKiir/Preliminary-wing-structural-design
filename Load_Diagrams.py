import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, pi
from scipy.integrate import cumulative_trapezoid
from scipy import interpolate


# general values
span = 47.85
t_c = 0.1194
density = 77.4337752982149
root_area = 4.345851384
S = 228.9
sweep = 35* pi /180
# chord values
root_chord = 7.4
taper_ratio = 0.3
tip_chord = root_chord * taper_ratio
quarter_chord = 0.25
x_c_shear_center = 0.35
x_c_centroid = 0.42
x_c_flap = 1.23 * quarter_chord


# load case
choice  = int(input('1: Take off (n = 2.5) \n2: Cruise (n = 2.5)\n3: Ground (n = 0)\n4: Cruise (n = -1) \n5: Take off (n = -1)\nChoice: '))
if choice == 1:
    title = 'Take off (n = 2.5)'
    rho = 1.225
    v = 77
    alpha = 14
    dCl = 0.499
    n = 2.5
elif choice == 2:
    title = 'Cruise (n = 2.5)'
    rho = 0.31641
    v = 0.514444 * 470.326
    alpha = 5.8
    dCl = 0
    n = 2.5
elif choice == 3:
    title = 'Ground (n = 0)'
    rho = 1.225
    v = 0
    alpha = 0
    dCl = 0
    n = 0
elif choice == 4:
    title = 'Cruise (n = -1)'
    rho = 0.31641
    v = 0.514444 * 470.326
    alpha = 5.8
    dCl = 0
    n = -1
elif choice == 5:
    title = 'Take off (n = -1)'
    rho = 1.225
    v = 77
    alpha = 14
    dCl = 0.499
    n = -1

#pitch moment change with flaps
def dCm(y):
    dCm = np.zeros_like(y)
    y_positions = np.where((y > 3.5) & (y < 15.4))[0]
    dCm[y_positions] = dCl * (x_c_flap - quarter_chord)
    return dCm

# chord function
def c(y):
    c = root_chord + ((2 * (tip_chord - root_chord) / span) * y)
    return c

# engine values
eng_thrust = 375000
eng_weight = 7277 * 9.80665
y_engine =  (1/3) *(span/2)
x_engine = 3 + x_c_centroid*c(y_engine)
z_engine = 2


# mass distribution
def mass_distribution(y):
    mass_distribution = 9.81 * density * root_area * (
                (root_chord - y * ((root_chord - tip_chord) / (span / 2))) / root_chord) ** 2
    return mass_distribution


# second file outputs AoA = 0
C_L0 = 0.215252
C_M0 = -0.350077

file = "MainWing_a=0.00_v=10.00ms.csv"
data = np.genfromtxt(fname=file, dtype=None, delimiter=",", skip_header=39, skip_footer=903, encoding='unicode_escape')

ylist = []
Cllst = []
Cdlst = []
Cmlst = []

for i in range(18):
    index = float(data[i][0])
    ylist.append(index)

for i in range(18):
    index = float(data[i][3])
    Cllst.append(index)

for i in range(18):
    index = float(data[i][5])
    Cdlst.append(index)

for i in range(18):
    index = float(data[i][7])
    Cmlst.append(index)

# converting newly generated lists into array
ylist = np.array(ylist)
Cllst = np.array(Cllst)
Cdlst = np.array(Cdlst)
Cmlst = np.array(Cmlst)

lift_coef_dist = interpolate.interp1d(ylist, Cllst, kind='cubic', fill_value='extrapolate')
drag_coef_dist = interpolate.interp1d(ylist, Cdlst, kind='cubic', fill_value='extrapolate')
moment_coef_dist = interpolate.interp1d(ylist, Cmlst, kind='cubic', fill_value='extrapolate')

# second file outputs AoA = 20

C_L20 = 1.048435
C_M20 = -1.51776

file2 = "MainWing_a=20.00_v=10.00ms.csv"
data2 = np.genfromtxt(fname=file2, dtype=None, delimiter=",", skip_header=39, skip_footer=903,
                      encoding='unicode_escape')

ylist2 = []
Cllst2 = []
Cdlst2 = []
Cmlst2 = []

for i in range(18):
    index = float(data2[i][0])
    ylist2.append(index)

for i in range(18):
    index = float(data2[i][3])
    Cllst2.append(index)

for i in range(18):
    index = float(data2[i][5])
    Cdlst2.append(index)

for i in range(18):
    index = float(data2[i][7])
    Cmlst2.append(index)

# converting newly generated lists into array
ylist2 = np.array(ylist2)
Cllst2 = np.array(Cllst2)
Cdlst2 = np.array(Cdlst2)
Cmlst2 = np.array(Cmlst2)

lift_coef_dist2 = interpolate.interp1d(ylist2, Cllst2, kind='cubic', fill_value='extrapolate')
drag_coef_dist2 = interpolate.interp1d(ylist2, Cdlst2, kind='cubic', fill_value='extrapolate')
moment_coef_dist2 = interpolate.interp1d(ylist2, Cmlst2, kind='cubic', fill_value='extrapolate')

# spanwise position Y
step = 100000
y = np.linspace(0, ylist.max(), step)


def ylist_limit():
    return ylist.max()

# dynamic pressure
def dynamic_pres(rho, V):
    q = 0.5 * rho * V ** 2
    return q


# lift critical values
def lift_dist(C_L, q, c):
    L = C_L * q * c
    return L



L1 = np.array(lift_dist(lift_coef_dist(y), dynamic_pres(rho, v), c(y)))
L2 = np.array(lift_dist(lift_coef_dist2(y), dynamic_pres(rho, v), c(y)))



# torque critical values
def moment_dist(C_M, q, S, c, span):
    M = C_M * q * S * c / span
    return M


M1 = np.array(moment_dist(moment_coef_dist(y), dynamic_pres(rho, v), S, c(y), span))
M2 = np.array(moment_dist(moment_coef_dist2(y), dynamic_pres(rho, v), S, c(y), span))


# lift distribution function at any angle of attack between 0 and 10 degrees
def lift_distribution_any(y, alpha, rho, v):
    C_Ld = C_L0 + (alpha * (C_L20 - C_L0) / 20)
    cl_distribution = dCl + lift_coef_dist(y) + ((C_Ld - C_L0) / (C_L20 - C_L0)) * (lift_coef_dist2(y) - lift_coef_dist(y))
    lift_distribution = n * np.array(lift_dist(cl_distribution, dynamic_pres(rho, v), c(y)))
    return lift_distribution

#pitch moment change with flaps
def dCm(y):
    dCm = np.zeros_like(y)
    y_positions = np.where((y > 3.5) & (y < 15.4))[0]
    dCm[y_positions] = dCl * (x_c_flap - quarter_chord)
    return dCm

# torque distribution function at any angle of attack between 0 and 10 degrees
def torque_distribution_any(y, alpha, rho, v):
    C_Md = C_M0 + (alpha * (C_M20 - C_M0) / 20)
    cm_distribution = -dCm(y) + moment_coef_dist(y) + ((C_Md - C_M0) / (C_M20 - C_M0)) * (
                moment_coef_dist2(y) - moment_coef_dist(y))
    torque_distribution = n * np.array(moment_dist(cm_distribution, dynamic_pres(rho, v), S, c(y), span))
    return torque_distribution


# engine step functions
def engine_weight(y):
    engine_weight = np.zeros_like(y)
    y_positions = y > y_engine
    engine_weight[y_positions] = eng_weight
    return engine_weight


def engine_thrust(y):
    engine_thrust = np.zeros_like(y)
    y_positions = y > y_engine
    if choice != 3:
        engine_thrust[y_positions] = eng_thrust
    return engine_thrust


# shear function
def shear(y, alpha, rho, v):
    L = lift_distribution_any(y, alpha, rho, v)
    integrated_lift = -cumulative_trapezoid(L, y, initial=0)
    integrated_mass = cumulative_trapezoid(mass_distribution(y), y, initial=0)
    shear = integrated_lift + integrated_mass + engine_weight(y)
    shear = shear - shear[-1]
    return shear


# bending function
def bending(y, alpha, rho, v):
    integrated_shear = cumulative_trapezoid(shear(y, alpha, rho, v), y, initial=0) - engine_thrust(y) * sin(sweep) * z_engine
    moment = integrated_shear - integrated_shear[-1]
    return moment

def bending_inter(y_values, y):
    return np.interp(y, y_values, bending(y_values, alpha, rho, v))

# distance from centroid to shear center
def d_weight(y):
    chord = c(y)
    d_weight = chord * (-x_c_shear_center + x_c_centroid) * cos(3.14 * alpha / 180)
    return d_weight


# torsion function
def torsion(y, alpha, rho, v):
    distance = d_weight(y)
    T = -torque_distribution_any(y, alpha, rho, v)
    integrated_factor = T - mass_distribution(y) * distance
    integrated_torsion = -cumulative_trapezoid(integrated_factor, y, initial=0) - engine_weight(
        y) * x_engine * cos(sweep) + engine_thrust(y) * cos(sweep) * z_engine
    torsion = integrated_torsion - integrated_torsion[-1]
    return torsion

def torsion_inter(y_values, y):
    return np.interp(y, y_values, torsion(y_values, alpha, rho, v))


# subplot 1 - shear
plt.subplot(3, 1, 1)
plt.suptitle(f"Internal Load Diagrams - {title}")
plt.fill_between(y, shear(y, alpha, rho, v) / 1000, color='lightskyblue', alpha=0.7)
plt.gca().spines['right'].set_visible(False)
if shear(y, alpha, rho, v)[0] > 0:
    plt.gca().spines['top'].set_visible(False)
    plt.gca().xaxis.tick_bottom()
else:
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().xaxis.tick_top()
plt.gca().xaxis.set_label_position("top")
plt.ylabel("Shear (kN)")

# subplot 2 - bending
plt.subplot(3, 1, 2)
plt.fill_between(y, bending(y, alpha, rho, v) / 1000, color='lightgreen', alpha=0.7)
plt.gca().spines['right'].set_visible(False)
if bending(y, alpha, rho, v)[1] > 0:
    plt.gca().spines['top'].set_visible(False)
    plt.gca().xaxis.tick_bottom()
else:
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().xaxis.tick_top()
plt.gca().xaxis.set_label_position("top")
plt.ylabel("Bending (kNm)")

# subplot 3 - torsion
plt.subplot(3, 1, 3)
plt.fill_between(y, torsion(y, alpha, rho, v) / 1000, color='salmon', alpha=0.7)
plt.gca().spines['right'].set_visible(False)
if torsion(y, alpha, rho, v)[0] > 0:
    plt.gca().spines['top'].set_visible(False)
    plt.gca().xaxis.tick_bottom()
else:
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().xaxis.tick_top()
plt.xlabel("Spanwise Position (m)")
plt.ylabel("Torsion (kNm)")

# plot generation
plt.tight_layout()
plt.show()