import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import Load_Diagrams as ld
from Functions import *

design_philosophy_nr = int(input("Enter design philosophy nr:"))

# Constants from previous workpages
b = 47.85 #span [m]
hb = b / 2 #half span
c_r = 7.4 #root chord [m]
c_t = 2.22 #tip chord [m]
y_end_add_sp = hb/3 + 1 # the y choordinate of the end of additional spar

# Material constants
E = 72.4e9 #Pa
G = 28e9 #Pa

# Wingbox coords in fraction of c
x_start = 0
x_end = 0.4
x_additional_spar = 0.2
z_lower_left = 0
z_lower_right = 0.0033
z_upper_left = 0.0974
z_upper_right = 0.1050

# Interactive parameters
num_str_top = 0 # Number of equally distributed stringers on the top side
num_str_bot = 0 # Number of equally distributed stringers on the bottom side
num_tot = num_str_top + num_str_bot # Total number of stringers
str_area = 0.003 #area of each individual stringer [m^2]
thickness_m = [0,0,0,0,0] # skin and spar thickness [m] as a function of root chord (skin: lower, upper; spar: left, right, additional)

# Plotting parameters
num_i = 1000
y_values = np.linspace(0, ld.ylist_limit(), num_i+1)

max_deflection = 0 #max wingtip delfection [m]
max_twist = 0 #max wing twist [deg]
safety_factor = 1.5

# Wingbox corners
wingbox_corners = np.array([
    [x_start, z_lower_left],
    [x_start, z_upper_left],
    [x_end, z_upper_right],
    [x_end, z_lower_right],
    [x_start, z_lower_left]
])

if design_philosophy_nr == 1:
    num_str_top = 18
    num_str_bot = 18
    num_tot = num_str_top + num_str_bot + 4
    thickness_m = [0.0015, 0.0015, 0.0027, 0.0027, 0.0027] #thickness in meters
    thickness = [thickness_m/c_r for thickness_m in thickness_m] #thickness as a function of root chord
    thickness_a = thickness
    str_area = 0.002
elif design_philosophy_nr == 2:
    num_str_top = 16
    num_str_bot = 16
    num_tot = num_str_top + num_str_bot + 4
    thickness_m = [0.008, 0.008, 0.015, 0.015, 0.015] #thickness in meters
    thickness = [thickness_m/c_r for thickness_m in thickness_m] #thickness as a function of root chord
    thickness_a = thickness
    str_area = 0.002
elif design_philosophy_nr == 3:
    num_str_top = 17
    num_str_bot = 17
    num_tot = num_str_top + num_str_bot + 4
    thickness_m = [0.0025, 0.0025, 0.005, 0.005, 0.005] #thickness in meters
    thickness = [thickness_m/c_r for thickness_m in thickness_m] #thickness as a function of root chord
    thickness_a = thickness
    str_area = 0.002
else:
    raise ValueError("The case number is wrong")


# Stringer and centroid positions
str_pos_top = position_of_stringers(num_str_top, x_start, z_upper_left, x_end, z_upper_right)
str_pos_bot = position_of_stringers(num_str_bot, x_start, z_lower_left, x_end, z_lower_right)
centroid_x, centroid_z = centroid(x_end, z_upper_left, z_upper_right, z_lower_right, thickness_a[1], thickness_a[0],
                                  thickness_a[2], thickness_a[3], x_additional_spar, thickness_a[4], num_tot,
                                  num_str_top, num_str_bot, str_area, c_t)

def wingbox_plot():
    fig, ax = plt.subplots(figsize=(10, 6))  # Plotting set-up
    ax.plot(wingbox_corners[:, 0], wingbox_corners[:, 1], 'k-', label="Wingbox")  # Plot wingbox

    # Plot stringers
    ax.scatter(str_pos_top[:, 0], str_pos_top[:, 1], color='blue', label="Top Stringers")
    ax.scatter(str_pos_bot[:, 0], str_pos_bot[:, 1], color='red', label="Bottom Stringers")

    # Plot centroid
    ax.scatter(centroid_x, centroid_z, color='green', s=100, label="Centroid")

    # Labels and legend
    ax.set_xlabel("x/c")
    ax.set_ylabel("z/c")
    ax.set_title(f"Wingbox with stringers and centroid locations (design {design_philosophy_nr})")
    ax.legend()
    ax.grid(True)
    ax.axis('equal')

    plt.savefig("wingbox.svg", format="svg")
    plt.show()

wingbox_plot()

def moment_diagram_plot():
    ixx = np.zeros(num_i + 1)
    izz = np.zeros(num_i + 1)

    for i in range(num_i + 1):
      y = i * hb / num_i
      c = c_r - (c_r - c_t) * y / hb
      thickness2 = thickness
      thickness2 = np.array(thickness2)
      thickness2 *= c
      centroid_x, centroid_z = centroid(x_end, z_upper_left, z_upper_right, z_lower_right, thickness_a[1], thickness_a[0],
                                          thickness_a[2], thickness_a[3], x_additional_spar, thickness_a[4], num_tot,
                                          num_str_top, num_str_bot, str_area, c)
      if y <= y_end_add_sp:
          add_z_pos = 0.0033 * (0.2 / 0.4)
          sides_pos = [[[0.0, 0.0], [0.0, 0.1050]], [[0.0, 0.1050], [0.4, 0.1050]], [[0.4, 0.1050], [0.4, 0.0033]],
                       [[0.4, 0.0033], [0.0, 0.0]], [[0.2, add_z_pos], [0.2, 0.1050]]]  # For each side [[x_s, z_s], [x_e, z_e]] - positions given as chord perecentage
      else:
          sides_pos = [[[0.0, 0.0], [0.0, 0.1050]], [[0.0, 0.1050], [0.4, 0.1050]], [[0.4, 0.1050], [0.4, 0.0033]],
                       [[0.4, 0.0033], [0.0, 0.0]]]  # For each side [[x_s, z_s], [x_e, z_e]] - positions given as chord perecentage
      ixx[i], izz[i] = moment_main(c, centroid_x, centroid_z, num_str_bot, num_str_top, str_area, thickness, sides_pos)

    # Plotting results
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    ax1.plot(y_values, ixx)
    ax1.set_title(f"Plot of $I_{{xx}}$ vs y (design {design_philosophy_nr})")
    ax1.set_xlabel("y")
    ax1.set_ylabel("Ixx")

    ax2.plot(y_values, izz)
    ax2.set_title(f"Plot of $I_{{zz}}$ vs y (design {design_philosophy_nr})")
    ax2.set_xlabel("y")
    ax2.set_ylabel("I_{{zz}}")

    plt.tight_layout()
    plt.savefig("moment_of_inertia.svg", format="svg")
    plt.show()

    return ixx, izz

ixx, izz = moment_diagram_plot()


#---------------TORSIONAL STIFFNESS CALCULATIONS----------------#

def wingbox_gemoetry():
    # Wingbox geometry with each wall named
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(wingbox_corners[:, 0], wingbox_corners[:, 1], 'k-', label="Wingbox")
    ax.plot([0.2, 0.2], [0.0033, 0.102], color='orange', label="Spar")
    ax.text(0.1, -0.009, 'a', fontsize=12, color='black')  # Bottom wall left
    ax.text(0.205, 0.05, 's', fontsize=12, color='black')  # Spar
    ax.text(0.3, -0.009, 'b', fontsize=12, color='black')  # Bottom wall right
    ax.text(0.3, 0.110, 'd', fontsize=12, color='black')  # Top wall right
    ax.text(0.1, 0.102, 'e', fontsize=12, color='black')  # Top wall left
    ax.text(-0.01, 0.05, 'f', fontsize=12, color='black')  # Left wall
    ax.text(0.405, 0.05, 'g', fontsize=12, color='black') #Right wall

    ax.set_xlabel("x/c")
    ax.set_ylabel("z/c")
    ax.set_title(f"Wingbox geometry with the spar")
    ax.legend()
    ax.grid(True)
    ax.axis('equal')
    plt.savefig("wingbox_geometry.svg", format="svg")
    plt.show()

wingbox_gemoetry()

def torsional_constant():
    q_1values = []
    q_2values = []
    twist_rate_values = []
    J_yv = []

    for i in range(num_i + 1):
        y = i * hb / num_i
        c = c_r - ((c_r - c_t) * y / hb)

        # The wingbox
        a = 0.2 * c  # left cell, bottom wall
        b = 0.2 * c  # right cell, bottom wall
        s = 0.0996 * c  # middle spar
        d = 0.2 * c  # right cell, top wall
        e = 0.2 * c  # left cell, top wall
        f = 0.0974 * c  # left wall
        g = 0.1017 * c  # right wall
        t_1 = thickness_m[0]  # top, bottom
        t_2 = thickness_m[2]  # left, right
        t_3 = thickness_m[4]  # middle spar
        t = 1  # unit torque

        A_1 = 0.5*(f+s)*0.2*c   # area of left cell
        A_2 = 0.5*(s+g)*0.2*c  # area of right cell

        if y < y_end_add_sp:
            #
            # MULTI-CELL
            #
            q_1, q_2, twist_rate = torsion_const_solve_system(A_1, A_2, a, b, s, d, e, t_1, t_2, t_3, G, t, f, g)

            q_1values.append(q_1)
            q_2values.append(q_2)
            twist_rate_values.append(twist_rate)

            J_y = t / (G * twist_rate)
            J_yv.append(J_y)
        else:
            #
            # SINGLE-CELL
            #
            A = A_1 + A_2
            J_y = 4 * A**2 / ((a + b + e + d) / t_1 + (g + f) / t_2)
            J_yv.append(J_y)

    plt.plot(y_values, J_yv, color='red', label=r'$J(y)$')
    plt.title(f'Torsional constant along the half wingspan (design {design_philosophy_nr})')
    plt.xlabel('y')
    plt.ylabel('J(y)')
    plt.grid()
    plt.tight_layout()
    plt.savefig("torsion_constant.svg", format="svg")
    plt.show()

    return J_yv

J_yv = torsional_constant()

#---------------------BENDING AND TWIST DISTRIBUTION------------------------#

#------------------------LOAD CASES------------------------#

#Deflection
def calculate_deflection_and_twist_general(y_values, ixx, J_yv, E, G, ld, case_number):

    I_xx_func_case1 = sp.interpolate.interp1d(y_values, ixx, kind="cubic", fill_value="extrapolate")

    def get_moment_values1(y_values, y):
        return ld.bending_inter(y_values, y)

    moment_values1 = [get_moment_values1(y_values, y) for y in y_values]
    moment_func1 = sp.interpolate.interp1d(y_values, moment_values1, kind="cubic", fill_value="extrapolate")

    # 1. Define dv/dy as an integral
    def calc_dv_dy_case1(y):
        return - sp.integrate.quad(lambda y: moment_func1(y) / (E * I_xx_func_case1(y)), 0, y)[0]

    # Create an array of y-values for dv/dy

    dv_dy_values_case1 = [calc_dv_dy_case1(y) for y in y_values]

    # Interpolate dv/dy to make it callable as a continuous function
    dv_dy_func_case1 = sp.interpolate.interp1d(y_values, dv_dy_values_case1, kind="cubic", fill_value="extrapolate")

    # 2. Define v(y) as the integral of dv/dy
    def calc_deflection_case1(y):
        return sp.integrate.quad(lambda x: dv_dy_func_case1(x), 0, y, epsabs=1e-4, epsrel=1e-4)[0]

    # Create an array of v(y) values by evaluating the integral
    v_y_values_case1 = [calc_deflection_case1(y) for y in y_values]
    # Interpolate v(y) to make it callable as a continuous function
    v_y_func_case1 = sp.interpolate.interp1d(y_values, v_y_values_case1, kind="cubic", fill_value="extrapolate")

    # Twist

    J_func_case1 = sp.interpolate.interp1d(y_values, J_yv, kind="cubic", fill_value="extrapolate")

    def get_torsion_values1(y_values, y):
        return ld.torsion_inter(y_values, y)

    torsion_values1 = [get_torsion_values1(y_values, y) for y in y_values]
    torsion_func1 = sp.interpolate.interp1d(y_values, torsion_values1, kind="cubic", fill_value="extrapolate")

    def calc_theta_case1(y):
        return sp.integrate.quad(lambda y: torsion_func1(y) / (G * J_func_case1(y)), 0, y)[0]

    # Create an array of theta(y) values by evaluating the integral
    theta_values_case1 = [math.degrees(calc_theta_case1(y)) for y in y_values]
    # Interpolate v(y) to make it callable as a continuous function
    theta_func_case1 = sp.interpolate.interp1d(y_values, theta_values_case1, kind="cubic", fill_value="extrapolate")
    return v_y_values_case1, theta_values_case1


def plotting_cases():
    global max_deflection, max_twist
    # Get values
    v_y_func_case1, theta_func_case1 = calculate_deflection_and_twist_general(y_values, ixx, J_yv, E, G, ld, case_number=1)
    max_deflection = max(v_y_func_case1, key=abs)
    max_twist = max(theta_func_case1, key=abs)

    # Plot everything
    plt.figure()

    plt.plot(y_values, v_y_func_case1)
    plt.title(f'v(y) (design {design_philosophy_nr}) - {ld.title}')
    plt.xlabel('y [m]')  # Label for x-axis
    plt.ylabel('v(y) [m]')  # Label for y-axis
    plt.yticks(np.arange(round(min(v_y_func_case1), 1), round(max(v_y_func_case1), 1), 0.5))  # Set y scale
    plt.xticks(np.arange(0, hb, 5))  # Set x scale
    plt.grid()
    plt.savefig("bending.svg", format="svg")
    plt.show()

    plt.plot(y_values, theta_func_case1)
    plt.title(f' θ(y) (design {design_philosophy_nr}) - {ld.title}')
    plt.xlabel('y [m]')  # Label for x-axis
    plt.ylabel('θ(y) [deg]')  # Label for y-axis
    if design_philosophy_nr == 1:
        plt.yticks(np.arange(np.floor(min(theta_func_case1)/1) * 1, np.ceil(max(theta_func_case1)/1) * 1, 1))  # Set y scale
    elif design_philosophy_nr == 2:
        plt.yticks(np.arange(np.floor(min(theta_func_case1) / 0.2) * 0.2, np.ceil(max(theta_func_case1) / 0.2) * 0.2, 0.2))  # Set y scale
    else:
        plt.yticks(np.arange(np.floor(min(theta_func_case1) / 0.5) * 0.5, np.ceil(max(theta_func_case1) / 0.5) * 0.5, 0.5))  # Set y scale
    plt.xticks(np.arange(0, hb, 5))  # Set x scale
    plt.grid()
    plt.savefig("twist.svg", format="svg")
    plt.show()

plotting_cases()

def Area_wingbox(x_e, z_u_s, z_u_e, z_l_e, t_u, t_l, t_L, t_R, x_a, t_a, number, chord):  # [ _s = start, _e = end, _u = upper, _l = lower, t = thickness, _L = left, _R = Right] positions of all corners
    """
    Calculates the total area of the wingbox actually occupied by material at a certain chord position
    This function can later be used to help decide which design option is best
    """

    angle_u = abs(math.atan((z_u_e - z_u_s) / x_e))
    angle_l = math.atan(z_l_e / x_e)
    # additional spar
    xc_a = x_a
    length_a = z_u_s + x_a * math.tan(angle_u) + x_a * math.tan(abs(angle_l))
    A_a = length_a * t_a * chord**2
    # Area of each wing box spar component
    A_l = z_u_s * t_L * chord**2 # left
    A_r = (z_u_e - z_l_e) * t_R * chord**2  # right
    A_up = x_e / math.cos(angle_u) * t_u * chord**2  # upper
    A_low = x_e / math.cos(abs(angle_l)) * t_l * chord**2 # lower
    A_wingbox = A_l + A_r + A_up + A_low + A_a + str_area * number
    return A_wingbox

Area_wingbox_v = Area_wingbox(x_end, z_upper_left, z_upper_right, z_lower_right, thickness_a[1], thickness_a[0], thickness_a[2], thickness_a[3], x_additional_spar, thickness_a[4], num_tot, c_t)
print('')
print(f"The cross-sectional area of the wingbox at the wing tip (design philosophy {design_philosophy_nr}) is {round(Area_wingbox_v, 3)} m^2")
print(f'YOU CHOSE {ld.title}')
print(f"The maximum deflection (design philosophy {design_philosophy_nr}) is {round(max_deflection, 2)} m")
print(f"The maximum twist (design philosophy {design_philosophy_nr}) is {round(max_twist, 2)} deg")
print('')
print(f'Considering a safety factor of {safety_factor}:')
if max_deflection * safety_factor <= (b * 0.15):
    print(f'Design philosophy {design_philosophy_nr} WORKS for bending')
else:
    print(f'Design philosophy {design_philosophy_nr} DOES NOT WORK for bending')
if max_twist * safety_factor <= 10:
    print(f'Design philosophy {design_philosophy_nr} WORKS for torsion')
else:
    print(f'Design philosophy {design_philosophy_nr} DOES NOT WORK for torsion')
