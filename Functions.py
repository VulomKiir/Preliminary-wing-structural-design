import numpy as np
import scipy as sp
import math
import Load_Diagrams as ld

def points_translation(x, z, c_x, c_z):
    x_new = x - c_x
    z_new = z - c_z

    return x_new, z_new

def position_of_stringers(number, x_s, z_s, x_e, z_e):
    pos = np.zeros((number, 2))  # Initialize array for stringer positions.

    for i in range(number):
        # Linearly interpolate x and z positions between start and end points.
        pos[i][0] = x_s + i * (x_e - x_s) / (number - 1)
        pos[i][1] = z_s + i * (z_e - z_s) / (number - 1)

    return pos

def centroid(x_e, z_u_s, z_u_e, z_l_e, t_u, t_l, t_L, t_R, x_a, t_a, number, num_top, num_bot, str_area, chord):  # [ _s = start, _e = end, _u = upper, _l = lower, t = thickness, _L = left, _R = Right] positions of all corners
    angle_u = abs(math.atan((z_u_e - z_u_s) / x_e))
    angle_l = math.atan(z_l_e / x_e)
    # additional spar
    xc_a = x_a
    length_a = z_u_s + x_a * math.tan(angle_u) + x_a * math.tan(abs(angle_l))
    zc_a = x_a * math.tan(angle_l) + length_a / 2
    A_a = length_a * t_a

    # Area of each wing box spar component
    A_l = z_u_s * t_L  # left
    A_r = (z_u_e - z_l_e) * t_R  # right
    A_up = x_e / math.cos(angle_u) * t_u  # upper
    A_low = x_e / math.cos(abs(angle_l)) * t_l  # lower
    A_wingbox = A_l + A_r + A_up + A_low + A_a

    # Centroid location z
    zc_l = z_u_s / 2  # centroid loc of left spar
    zc_r = (z_u_e + z_l_e) / 2
    zc_up = z_u_s + t_u / 2 + 0.5 * x_e * math.tan(angle_u)
    zc_low = t_l / 2 + 0.5 * x_e * math.tan(angle_l)
    lwz = zc_l * A_l + zc_r * A_r + zc_up * A_up + zc_low * A_low + zc_a * A_a  # todo SUM of weight 1 --Z--

    # Centroid location x
    xc_l = 0
    xc_r = x_e
    xc_up = x_e / 2
    xc_low = x_e / 2
    lwx = xc_l * A_l + xc_r * A_r + xc_up * A_up + xc_low * A_low + xc_a * A_a  # todo SUM of weight 2 --X--

    # WEIGHT Area centroid product for corner spars of area A
    A_i_spar = (chord / math.sqrt(str_area)) ** 2
    x_c_low_l = 0 * A_i_spar
    z_c_low_l = 0 * A_i_spar
    x_c_up_l = 0 * A_i_spar
    z_c_up_l = z_u_s * A_i_spar
    x_c_low_r = x_e * A_i_spar
    z_c_low_r = x_e * math.tan(angle_l) * A_i_spar
    x_c_up_r = x_e * A_i_spar
    z_c_up_r = (z_u_s + x_e * math.tan(angle_u)) * A_i_spar
    lswx = x_c_low_l + x_c_low_r + x_c_up_l + x_c_up_r  # todo SUM of weight 3 --X--
    lswz = z_c_low_l + z_c_low_r + z_c_up_l + z_c_up_r  # todo SUM of weight 4 --Z--
    # print("lswz: ", z_c_low_l, z_c_low_r, z_c_up_l, z_c_up_r)

    # WEIGHT Area centroid product for all other spars of total area A
    A_r_bot = num_bot * A_i_spar
    A_r_top = num_top * A_i_spar
    x_c_bot_rest = x_e / 2 * A_r_bot
    x_c_top_rest = x_e / 2 * A_r_top
    z_c_bot_rest = x_e / 2 * math.tan(angle_l) * A_r_bot
    z_c_top_rest = (z_u_s + x_e / 2 * math.tan(angle_u)) * A_r_top
    lsswx = x_c_bot_rest + x_c_top_rest  # todo SUM of weight 5 --X--
    lsswz = z_c_bot_rest + z_c_top_rest  # todo SUM of weight 6 --Z--
    # print("lsswz: ", z_c_bot_rest, z_c_top_rest)

    # TOTAL area
    A_tot = A_i_spar * number + A_wingbox

    # print(lwz, lswz, lsswz)
    # CENTROID
    centroid_x = (lwx + lswx + lsswx) / A_tot
    centroid_z = (lwz + lswz + lsswz) / A_tot

    return centroid_x, centroid_z

def moment_of_inertia(sides, str_pos_bot, str_pos_top, str_area, a, thickness):  # Calculating moment of inertia for a specified position y
  # sides: numpy Array, dim: 3?
  # a defines axis, a = 1 for Ixx // a = 0 for Izz

  ixx = 0  # ixx/izz

  # Stringers - Steiners term
  for i in range(len(str_pos_bot)):
    ixx += str_area * str_pos_bot[i][a] ** 2
  for j in range(len(str_pos_top)):
    ixx += str_area * str_pos_top[j][a] ** 2

  # Calculating moment of inertia of sides of the box - thin walled assumptions
  for k in range(len(sides)):  # [[x_s, z_s], [x_e, z_e]]
    area = thickness[k] * np.linalg.norm(sides[k][0] - sides[k][1])
    z_avg = (sides[k][1][a] - sides[k][1][a]) / 2

    # Main inertia calc
    lens = np.linalg.norm(sides[k][0] - sides[k][1])
    #print(lens, sides[k][0], sides[k][1])
    Bet = (sides[k][1][a] + sides[k][0][a]) / lens  # (0.5*(y or x))/(0.5 * dist)  z - Ixx, x - Izz
    ixx += (lens ** 3) * thickness[k] * (Bet ** 2) / 12

    # Steiners term
    ixx += area * z_avg ** 2

  return ixx

def moment_main(c, x_c, z_c, num_str_bot, num_str_top, str_area, thickness, sides_pos):

    # x, z positions of the sides
    top_side = 1 # Index in slides_pos
    bottom_side = 3  # Index in slides_pos
    sides_pos = np.array(sides_pos)


    for s in range(len(sides_pos)):
        sides_pos[s][1] *= c
        sides_pos[s][0] *= c
    
    x_c *= c
    z_c *= c

    # Translating to centroid coord. sys.
    for i in range(len(sides_pos)):
        sides_pos[i][0][0], sides_pos[i][0][1] = points_translation(sides_pos[i][0][0], sides_pos[i][0][1], x_c, z_c)
        sides_pos[i][1][0], sides_pos[i][1][1] = points_translation(sides_pos[i][1][0], sides_pos[i][1][1], x_c, z_c)

    # Getting the positions of stringers
    str_pos_bot = position_of_stringers(num_str_bot, sides_pos[bottom_side][0][0], sides_pos[bottom_side][0][1], sides_pos[bottom_side][1][0], sides_pos[bottom_side][1][1])
    str_pos_top = position_of_stringers(num_str_top, sides_pos[top_side][0][0], sides_pos[top_side][0][1], sides_pos[top_side][1][0], sides_pos[top_side][1][1])

      # Calculating the moment of inertia
    ixx = moment_of_inertia(sides_pos, str_pos_bot, str_pos_top, str_area, 1, thickness)
    izz = moment_of_inertia(sides_pos, str_pos_bot, str_pos_top, str_area, 0, thickness)
    return ixx, izz

def torsion_const_solve_system(A_1, A_2, a, b, s, d, e, t_1, t_2, t_3, G, t, f, g):
    leftside = np.array([
        [(1 / (2 * A_1)) * ((s / (G * t_3)) + ((e + f + a) / (G * t_1))), (1 / (2 * A_1)) * (-s / (G * t_3)), -1],
        [(1 / (2 * A_2)) * (s / (G * t_3)), (1 / (2 * A_2)) * (((g + d + b) / (G * t_2)) + (-s / (G * t_3))), -1],
        [2 * A_1, 2 * A_2, 0]], )  # q1,q2,dtdy

    rightside = np.array([0, 0, t])
    solution = np.linalg.solve(leftside, rightside)
    q_1, q_2, twist_rate = solution
    return q_1, q_2, twist_rate
