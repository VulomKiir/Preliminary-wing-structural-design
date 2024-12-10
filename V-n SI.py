import math
import numpy as np
import matplotlib.pyplot as plt

''' FOR 2.0, I WILL NOT USE CN'''

Altitude = int(input('Type 1 for Sea Level \n Type 2 for Cruise Altitude \n'))

if Altitude == 1:
    rho = 1.225
elif Altitude == 2:
    rho = 0.316

#Parameters
Match_W_S_Ratio = 6086                                          #WS ratio obtained from matching diagram, N/m^2
S = 319.1

# Mach, FL, Temp = 0.82, 39000, 216.65                          #Flight conditions
CL_max = 1.63
CL_min = -1.67 * 0.8                                            #Using DATCOM method together with graphs from WP2
CL_max_HLD = 2.066                                              #Maximum C_Lmax attained by aircraft, with HLD devices extended

OEW = 65888                                                     #Operating Empty Mass
MTOW = 150566                                                   #Maximum TakeOff Weight
MPW = 49442                                                     #Maximum Payload Weight

def a(Temp):
    return math.sqrt(1.4*287*Temp)

#Speeds in TAS, m/s

if Altitude == 1:
    V_C = 0.82 * a(216.65) * math.sqrt(0.316/1.225)
    V_D = 1.25*V_C
elif Altitude == 2:
    V_C = 0.82 * a(216.65)
    V_D = V_C + 0.05*a(216.65)

#TAS to EAS
TAStoEAS = math.sqrt(rho/1.225)

#Finding n_max & n:_min, using FAR25
DesMTOW = Match_W_S_Ratio * S / 9.80665 * 2.2046226218487757

n_max = 2.1 + 2400/(DesMTOW+10000)
print(DesMTOW)
print(n_max)

if n_max < 2.5:
    n_max = 2.5
elif n_max > 3.8:
    n_max = 3.8

n_min = -1

for i in [OEW, MTOW, OEW + MPW]:                 #Draw manoeuver diagram for 3 given

    WS_ratio = i/S

    #Speeds
    V_S0 = math.sqrt(2*WS_ratio/(rho*CL_max_HLD))            #Stall speed with flaps extended
    V_S1 = math.sqrt(2*WS_ratio/(rho*CL_max))                #Stall speed with flaps retracted

    #Maneuver load diagram
    V1_array = []
    N1_array = []

    V2_array = []
    N2_array = []

    V3_array = []
    N3_array = []

    V = np.linspace(0,400, 10000)

    for v in V:
        n = (0.5 * rho * (v)**2 * CL_max)/WS_ratio      #Test some speeds until speed at which n=2.5 is reached
        V1_array = np.append(V1_array, v)               #Append values to array for plotting
        N1_array = np.append(N1_array, n)

        if n >= n_max:
            if v <= V_S1 * n_max**0.5:                  #Minimum value for V_A
                V_A = V_S1 * n_max**0.5
            else:
                V_A = v           #If v at n=2.5 is greater than minimum, then assign that value as V_A
            break
        if n <= 2:
            Vn2 = v               #Speed at which n = 2 for flaps down

    for v in V:
        n = (0.5 * rho * (v)**2 * CL_min)/WS_ratio      #Draw line for negative loadings
        V2_array = np.append(V2_array, v)
        N2_array = np.append(N2_array, n)
    
        if n <= n_min:
            V_Sm1 = v
            break

    for v in V:
        n = (0.5 * rho * (v)**2 * CL_max_HLD)/WS_ratio
        V3_array = np.append(V3_array, v)
        N3_array = np.append(N3_array, n)
    
        if n >= 2:
            break

    #Conversion
    mtokn = 1.94384

    V1_array = np.append(V1_array, [V_C, V_D, V_D]) * TAStoEAS  * mtokn       #Convert values to KEAS
    N1_array = np.append(N1_array, [n_max, n_max, 0])

    V2_array = np.append(V2_array, [V_C, V_D]) *TAStoEAS  * mtokn            #Convert values to KEAS
    N2_array = np.append(N2_array, [n_min, 0])

    V3_array = np.append(V3_array, [Vn2]) * TAStoEAS * mtokn                  #Convert values to KEAS
    N3_array = np.append(N3_array, [2])

    plt.figure(figsize=(8, 6))

    plt.plot(V1_array, N1_array, 'blue')
    plt.plot(V2_array, N2_array, 'blue')
    plt.plot(V3_array, N3_array, 'blue')
    plt.plot([0, V_S1 * TAStoEAS * 1.94384 ], [1, 1], 'gray', linestyle='dashed')
    plt.plot([V_S0 * TAStoEAS * 1.94384 , V_S0 * TAStoEAS * 1.94384 ], [1, 0], 'gray', linestyle='dashed')
    plt.plot([V_Sm1 * TAStoEAS * 1.94384 , V_Sm1 * TAStoEAS * 1.94384 ], [-1, 0], 'gray', linestyle='dashed')
    plt.plot([V_S1 * TAStoEAS * 1.94384 , V_S1 * TAStoEAS * 1.94384 ], [1, 0], 'gray', linestyle='dashed')
    plt.plot([V_A * TAStoEAS * 1.94384 , V_A * TAStoEAS * 1.94384 ], [2.5, 0], 'gray', linestyle='dashed')
    plt.plot([V_C * TAStoEAS * 1.94384 , V_C * TAStoEAS * 1.94384 ], [0, -1], 'gray', linestyle='dashed')
    plt.plot([0, V_D * TAStoEAS * 1.94384 ], [0, 0], 'black', linestyle='dashed')
    plt.xlim(0,320)
    plt.ylim(-1.2,3.2)
    plt.xlabel('Airspeed, KEAS')
    plt.ylabel('Load factor, n')

    #Text on graph
    plt.text(V_S0* TAStoEAS * 1.94384 - 10, -0.15, '$V_{S0}$', fontsize = 8)
    plt.text(V_Sm1* TAStoEAS * 1.94384 +2, -0.9, '$V_{S{-}1}$', fontsize = 8)
    plt.text(V_S1* TAStoEAS * 1.94384 + 2, 0.1, '$V_{S1}$', fontsize = 8)
    plt.text(V_A* TAStoEAS * 1.94384 + 2, 0.1, '$V_{A}$', fontsize = 8)
    plt.text(V_C* TAStoEAS * 1.94384 - 3, 0.1, '$V_{C}$', fontsize = 8)
    plt.text(V_D* TAStoEAS * 1.94384 + 2, 0.1, '$V_{D}$', fontsize = 8)

    print( f'V_A - V_D = {V_A * TAStoEAS * mtokn} - {V_D * TAStoEAS * mtokn}')
    print( f'V_S-1 - V_C = {V_Sm1 * TAStoEAS * mtokn} - {V_C * TAStoEAS * mtokn}')
    print( f'V_S0, V_S1 = {V_S0 * TAStoEAS * mtokn}, {V_S1 * TAStoEAS * mtokn}')

    plt.grid()
    plt.show()