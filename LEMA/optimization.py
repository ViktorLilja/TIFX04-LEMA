import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from stage import Stage
from projectile import Projectile
from experiment import Experiment

magnetType = "50mm"     # Type of magnet used, 50mm or 2x20mm
initialSpeed = 3.0      # Initial speed of projectile

# Helper method for setting up experiment
def setupExperiment(param):
    turns = 100*param[0]
    startPos = 1e-2*param[1]

    # Set up experiment
    proj = Projectile(type=magnetType, 
                      x0=startPos, 
                      xdot0=initialSpeed)
    stages = [Stage(n=turns,
                    gap=proj.gap(),
                    dx=-10e-2,
                    uC0=300
                    )]
    return Experiment(stages, proj)


def minus_efficiency(param):
    experiment = setupExperiment(param)
    
    # Run simulation
    experiment.simulate()

    # Get results
    effiency = experiment.getEfficiency()
    print("Turns: %.2f, Startpos %.3f cm, Efficiency: %.1f%%" % 
         (100*param[0], param[1], 100*effiency))
    return -effiency

# RUN OPTIMIZATION
res = minimize(minus_efficiency, np.array([2.0, -3]))
if res.success: print("Optimization successful!")
else: print(res.message)


# SHOW RESULT
param = res.x
experiment = setupExperiment(param)

experiment.simulate()
t = experiment.getTime()
v = experiment.getSpeed()

# Plot
plt.subplot(2,1,1)
plt.plot(1000*t, v)
plt.ylabel("Speed [m/s]")

plt.subplot(2,1,2)
for i, stage in enumerate(experiment.stages):
    I = experiment.getCurrent(i)
    plt.plot(1000*t, I)
plt.ylabel("Current [A]")
plt.xlabel("Time [ms]")

plt.show()