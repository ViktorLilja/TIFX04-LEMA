# Experiment class, runs an experiment using given coil and projectile

import numpy as np
from scipy.integrate import solve_ivp
from stage import Stage

class Experiment:

    def __init__(self,
                 stages, 
                 projectile,
                 tspan = [0,0.05],
                 rtol = 1e-6
                 ):

        self.stages = stages
        self.proj = projectile
        self.tspan = tspan
        self.rtol = rtol

        # Latest simulation result
        self.result = None          # Object containing simulation result
        self.simFinished = False    # False if last simulation finished too early

    def simulate(self):
        # State variables:
        # x, xdot, vec_uC, vec_I
        # x:        projectile position
        # xdot:     projectile speed
        # vec_uC:   vector of capacitor voltages for all stages
        # vec_I:    vector of coil currents for all stages

        # Preparing empty arrays
        n = len(self.stages)
        x0 = np.zeros(shape=(1,))
        xdot0 = np.zeros(shape=(1,))
        vec_uC0 = np.zeros(shape=(n,))
        vec_I0  = np.zeros(shape=(n,))

        # Get initial values from components
        x0[0] = self.proj.x0
        xdot0[0] = self.proj.xdot0
        for i, stage in enumerate(self.stages):
            vec_uC0[i] = stage.uC0
            vec_I0[i] = stage.i0

        # Compile initial values
        y0 = np.concatenate((x0, xdot0, vec_uC0, vec_I0))
        
        # Solve ODE system
        self.result = solve_ivp(self.__ydot__, self.tspan, y0, rtol=self.rtol)

        # Check if capacitors are empty at end of simulation
        self.simFinished = True
        for i, stage in enumerate(self.stages):
            uC_final = self.getVoltage(i)[-1]
            if uC_final > 1: self.simFinished = False

    def getTime(self):
        if self.result is None: 
            print("Warning: returned None, call LEMA.simulate() before accessing values")
            return None
        
        return self.result.t

    def getPosition(self):
        if self.result is None: 
            print("Warning: returned None, call LEMA.simulate() before accessing values")
            return None
        
        return self.result.y[0,:]

    def getSpeed(self):
        if self.result is None: 
            print("Warning: returned None, call LEMA.simulate() before accessing values")
            return None
        
        return self.result.y[1,:]

    def getFinalVelocity(self):
        v = self.getSpeed()
        return v[-1]

    def getVoltage(self, stageNumber):
        if self.result is None: 
            print("Warning: returned None, call LEMA.simulate() before accessing values")
            return None
        
        return self.result.y[2+stageNumber,:]

    def getCurrent(self, stageNumber):
        if self.result is None: 
            print("Warning: returned None, call LEMA.simulate() before accessing values")
            return None

        return self.result.y[2+len(self.stages)+stageNumber,:]
    
    def getKineticEnergy(self):
        return 0.5 * self.proj.m * self.getSpeed()**2

    def getCapacitorEnergy(self):
        sum = np.zeros_like(self.result.t)
        for i, stage in enumerate(self.stages):
            sum += 0.5 * stage.C * self.getVoltage(i)**2
        return sum
    
    def getInductorEnergy(self):
        sum = np.zeros_like(self.result.t)
        for i, stage in enumerate(self.stages):
            sum += 0.5 * stage.L * self.getCurrent(i)**2
        return sum
    
    def getGeneratedHeat(self):
        power = np.zeros_like(self.result.t)
        for i, stage in enumerate(self.stages):
            power += stage.R * self.getCurrent(i)**2
        timesteps = np.append(0,np.diff(self.result.t))
        return np.cumsum(power * timesteps)

    def getEfficiency(self):
        v = self.getSpeed()
        # Energy accounting
        v_final = v[-1]
        v_start = v[0]
        E_k = 0.5 * self.proj.m * v_final**2 - \
              0.5 * self.proj.m * v_start**2

        E_C = 0
        for stage in self.stages:
            E_C += 0.5 * stage.C * stage.uC0**2

        if not self.simFinished:
            print("Warning: simulation did not finish in time, efficiency", 
                  "may be inaccurate. Try increasing tspan.")

        return E_k/E_C

    # State variable derivate for use with ODE solver
    def __ydot__(self, t, y):

        # Extract state variables
        n = len(self.stages)
        
        x       = y[0]                      # Current values
        xdot    = y[1]
        vec_uC  = y[2:2+n]
        vec_i   = y[2+n:2+2*n]

        dx_dt      = np.zeros(shape=(1,))   # Derivates to be calculated
        dxdot_dt   = np.zeros(shape=(1,))
        dvec_uC_dt = np.zeros(shape=(n,))
        dvec_I_dt  = np.zeros(shape=(n,))


        # Calculate derivates of state variables
        dx_dt[0] = xdot
        for i, stage in enumerate(self.stages):
            if not stage.active(x): continue

            # Derivative of mutual inductance
            dk_dx = self.proj.dk_dx(stage.n, x-stage.x)
            dM_dx = dk_dx * np.sqrt(stage.L*self.proj.L) 

            # Equations describing the system:

            # Velocity of projectile
            dxdot_dt[0] += dM_dx*self.proj.I*vec_i[i]/self.proj.m

            # Voltage of capacitor
            dvec_uC_dt[i] = - vec_i[i]/stage.C

            # Current through coil
            dvec_I_dt[i] = - self.proj.I*dM_dx*xdot/stage.L \
                           + vec_uC[i]/stage.L \
                           - stage.R*vec_i[i]/stage.L

        # Compile derivatives
        dy_dt = np.concatenate((dx_dt, dxdot_dt, dvec_uC_dt, dvec_I_dt))
        return dy_dt


# Sample experiment
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from stage import Stage
    from projectile import Projectile

    # Set up experiment
    proj = Projectile(type="50mm", x0=-2.5e-2, xdot0=2)
    stages = [Stage(n=400, gap=proj.gap(), x=0e-2,  dx=-2e-2, uC0=200),
              Stage(n=300, gap=proj.gap(), x=5e-2,  dx=-1.2e-2),
              Stage(n=200, gap=proj.gap(), x=10e-2, dx=-2e-2)]
    experiment = Experiment(stages, proj)

    # Run simulation
    experiment.simulate()

    # Get results
    t = experiment.getTime()
    v = experiment.getSpeed()

    effiency = experiment.getEfficiency()
    print("Effciency: %.1f%%" % (100*effiency))

    # Plotting
    plt.subplot(2,1,1)
    plt.plot(1000*t, v)
    plt.ylabel("Speed [m/s]")

    plt.subplot(2,1,2)
    for i, stage in enumerate(stages):
        I = experiment.getCurrent(i)
        plt.plot(1000*t, I)
    plt.ylabel("Current [A]")
    plt.xlabel("Time [ms]")

    plt.show()