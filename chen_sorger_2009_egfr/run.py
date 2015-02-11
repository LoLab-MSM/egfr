from erbb_exec import model
import numpy as np
from pysb.integrate import odesolve

t = np.linspace(0,2000, num=2000)

yout = odesolve(model, t, integrator='lsoda')