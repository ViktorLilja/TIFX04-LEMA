# Generate design guide plot by optimizing coil turns N
# and starting position dx for a range of initial
# velocities xdot0.

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from stage import Stage
from projectile import Projectile
from experiment import Experiment

magnetType = "50mm"     # Type of magnet used, 50mm or 2x20mm

def setupExperiment(param, const):
    turns    = 100*param[0]
    startPos = 1e-2*param[1]

    xdot0    = const[0]
    uC0      = const[1]

    # Set up experiment
    proj = Projectile(type=magnetType,
                      x0=startPos,
                      xdot0=xdot0)
    stages = [Stage(n=turns,
                    gap=proj.gap(),
                    uC0=uC0
                    )]
    return Experiment(stages, proj,tspan=[0, 0.05/(1+xdot0/10)])


def minus_efficiency(param, const):
    experiment = setupExperiment(param, const)
    
    # Run simulation
    experiment.simulate()

    # Get results
    effiency = experiment.getEfficiency()
    print("Turns: %.2f, Startpos %.3f cm, Efficiency: %.1f%%" % 
         (100*param[0], param[1], 100*effiency))
    return -effiency

# RUN OPTIMIZATION
n_points = 21
v0_range   = np.linspace(0,10,n_points)
n_range    = np.zeros_like(v0_range)
eta_range  = np.zeros_like(v0_range)
dx_range   = np.zeros_like(v0_range)
vfin_range = np.zeros_like(v0_range)


for i in range(0, n_points):

    # Constant values
    v0    = v0_range[i]  # Initial projectile velocity
    uC0   = 300          # Capacitor voltage
    const = np.array([v0, uC0])

    # Starting guess for parameters to optimize
    turns0 = 6/(1+v0/3)+0.4  # Number of coil turns in hundreds
    pos0   = -4.3   # Starting position relative coil in cm
    param0 = np.array([turns0, pos0])

    res = minimize(lambda param: minus_efficiency(param, const), param0)
    param = res.x

    if res.success: print("Optimization successful!")
    else: print(res.message)

    exp = setupExperiment(param, const)
    exp.simulate()
    vfin = exp.getFinalVelocity()

    n_range[i]     = param[0]
    dx_range[i]    = -param[1]
    eta_range[i]   = -100*res.fun
    vfin_range[i]  = vfin

plt.plot(v0_range, n_range, "k-")
plt.plot(v0_range, vfin_range-v0_range, "k.-")
plt.plot(v0_range, eta_range, "k--" )
plt.plot(v0_range, dx_range, "k-.")
#plt.xticks(range(0,11))
plt.grid()
plt.minorticks_on()
plt.grid(b=True, which='minor', color='lightgray')
plt.grid(b=True, which='major', color='gray')
plt.title("Optimal LEMA, ett spolpar")
plt.xlabel("Ingångsfart (m/s)")
plt.legend(["Spolvarv (hundratal)",
            "Fartökning (m/s)",
            "Verkningsgrad (%)",
            "Startposition (cm)"])
plt.show()


# Write csv for plotting
with open("LEMA/output/designplot.csv", "w") as f:
    
    # Header
    f.write("v0 (m/s),n (hundratal),deltav (m/s),eta (%), Dx (cm)\n")

    for i in range(0, len(v0_range)):
        f.write("%e," % v0_range[i])
        f.write("%e," % n_range[i])
        f.write("%e," % (vfin_range[i]-v0_range[i]))
        f.write("%e," % eta_range[i])
        f.write("%e\n" % dx_range[i])

# SHOW RESULT
experiment = setupExperiment(param, const)

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