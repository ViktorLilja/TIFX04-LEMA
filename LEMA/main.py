# Example code demonstrating how to simulate one accelerator stage

import matplotlib.pyplot as plt
import numpy as np
from stage import Stage
from projectile import Projectile
from experiment import Experiment

# Settings
magnetType = "50mm"     # Type of magnet used, 50mm or 2x20mm
initialSpeed = 0      # Initial speed of projectile

startPos = -3.5e-2      # Position of projectile relative to coil
turns = 670             # Turns of wire in coil

# Set up experiment
proj = Projectile(type=magnetType, x0=startPos, xdot0=initialSpeed)
stages = [Stage(n=turns,
                gap=proj.gap(),
                dx=-10e-2,
                uC0=300
                )]
experiment =  Experiment(stages, proj)

print("Inductance: %.2f mH" % (1000*stages[0].L))
print("Resistance: %.2f Ohm" % (stages[0].R))

# Get result
experiment.simulate()
t = experiment.getTime()
v = experiment.getSpeed()
print("Efficiency: %.2f%%" % (100*experiment.getEfficiency()))

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

plt.subplot(2,1,1)
plt.title("Energy accounting")
plt.plot(1000*t, experiment.getCapacitorEnergy())
plt.plot(1000*t, experiment.getInductorEnergy())
plt.plot(1000*t, experiment.getKineticEnergy())
plt.plot(1000*t, experiment.getGeneratedHeat())
plt.legend(["Capacitor", "Inductor", "Kinetic", "Heat"])
plt.ylabel("Energy [J]")

plt.subplot(2,1,2)
sum = np.zeros_like(t)
sum += experiment.getCapacitorEnergy()
sum += experiment.getInductorEnergy()
sum += experiment.getKineticEnergy()
sum += experiment.getGeneratedHeat()
plt.plot(1000*t, sum)
plt.ylabel("Total energy [J]")
plt.xlabel("Time [ms]")

plt.show()