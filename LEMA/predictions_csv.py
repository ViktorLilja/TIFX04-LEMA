# Predicts acceleration for a specific set of given parameters
# and prints result to csv file for plotting

import matplotlib.pyplot as plt
import numpy as np
from stage import Stage
from projectile import Projectile
from experiment import Experiment

magnetType = "50mm"     # Type of magnet used, 50mm or 2x20mm

# Requested simulations:
#            SP1  SP2  SP3  SP4  SP5
n         = [668, 470, 358, 313, 274]
v0        = [  0, 3.6, 5.8, 7.4, 8.8] # (m/s)
t_start   = [  0,  34,  53.812,  69,  83.735] # (ms) 
t_end     = [ 33,  52,  68,  83, 100] # (ms)
dx        = [3.4, 3.8, 4.0, 4.2, 4.2] # (cm)

times    = []
speeds   = []
currents = []

for i in range(0,5):

    # Set up experiment
    proj = Projectile(type=magnetType, x0=-dx[i]/100, xdot0=v0[i])
    stages = [Stage(n=n[i],
                    gap=proj.gap(),
                    dx=-10e-2,
                    uC0=300
                    )]
    experiment =  Experiment(stages, proj, tspan=[t_start[i]/1000, t_end[i]/1000])


    # Get result
    experiment.simulate()
    t = experiment.getTime()*1000
    v = experiment.getSpeed()
    i = experiment.getCurrent(0)

    times.append(t)
    speeds.append(v)
    currents.append(i)


with open('LEMA/output/simulated_acceleration.csv', "w") as f:
    # Header
    for i in range(1,6):
        f.write("Simulated time %d (ms), Simulated speed %d (m/s), Simulated current %d (A)," % (i, i, i))
        if i == 5: f.write("\n")
    
    # Data
    row = 0
    written = True
    while written:
        written = False
        for i in range(0,5):
            if row < len(times[i]):
                written = True
                f.write("%e, %e, %e," % (times[i][row], speeds[i][row], currents[i][row]))
            else:
                f.write(",,,")
        f.write("\n")
        row = row + 1