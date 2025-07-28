import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#defining the properties for the egg

k = 0.5 #thermal conductivity
rho = 1000 #density
c = 3000 #specific heat
alpha = k /( rho*c)

#defining the parameters for boundary and initial conditions

R = 0.02 #chicken egg
T_i = 5 #intial temp in celcius
T_f = 80 #final target temp
T_surface = 100 #temperature for the surface water
dt = 0.2 #time step
dr= 0.0005 #steps in spherical coordinates

r = np.arange(0,R+dr,dr)
N = len(r)
T = np.ones(N) * T_i # need to define my somewhat empty ones matrix

def forward_step (T, alpha, r, dr, dt):
    T_update = T.copy()
    for i in range (1, len(r) - 1):
        d2Tdr2 = (T[i+1] - 2*T[i] +T[i-1]) / (dr**2)
        dTdr = (T[i+1] - T[i-1]) / (2*dr)
        T_update[i] = T[i] + dt*alpha*(d2Tdr2 + ((2/r[i])*dTdr)) #isolate for T_i (k+1)

    T_update[0] = T_update[1]
    T_update[-1] = T_surface                                                                            
    return T_update
   
#in the main code you will have a while loop for the temperture

time = 0.0
max_time = 2400
T_rec = [T.copy()]
times = [time]
over_duration = 0.0

cook_time = None

while time < max_time:

    T = forward_step(T, alpha, r, dr, dt)
    time += dt
    T_rec.append(T.copy())
    times.append(time)

    if np.min(T[:-1]) >= T_f:
        over_duration += dt
        if over_duration >= 10.0:
            cook_time = time
            break

    else:
        over_duration = 0.0

#plotting

T_rec = np.array(T_rec)

# Plot temperature at the center of the egg over time

center_index = 0  
center_temps = [temp[center_index] for temp in T_rec]

plt.figure(figsize=(8, 5))
plt.plot(times, center_temps, label="Center Temperature")
plt.axhline(T_f, color='r', linestyle='--', label="Target Temp (80°C)")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (°C)")
plt.title("Egg Center Temperature Over Time")
plt.legend()
plt.tight_layout()
plt.show()
plt.close()

                                   
print(cook_time)
