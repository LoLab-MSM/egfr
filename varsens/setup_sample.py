# Load libraries
import numpy as np
from   egfr.erbb_exec import model
import egfr.parameter_dict_A431
import pickle
from   pysb.integrate import odesolve
from   varsens import *
import time

import os, errno

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

mkdir_p("varsens/samples")

# Load fitted parameters
with open('egfr/calibration_A431_highEGF_fittedparams_unnorm_splined_RESTART', 'rb') as handle:
    fittedparams = pickle.loads(handle.read())

# Set those fits in the model
for i in range(len(model.parameters)):
    if model.parameters[i].name in fittedparams:
        model.parameters[i].value = fittedparams[model.parameters[i].name]

# Generate Reference
t = np.linspace(0, 7200, 14401)
t0=time.time()
reference = odesolve(model, t)
exec_t = time.time() - t0

k = len(fittedparams) # Dimensions (i.e. parameters)
n = 5                 # Number of samples

# Just how bad is it going to be!
print "Evaluations required is %d\n" % (n*(2*k+2))
print "Hours estimated is %f\n" % (n*(2.0*k+2)*exec_t / 60.0 / 60.0)


# Given a set of parameters, determine outcome
def objective(x):
    for i in range(len(x)):
        model.parameters[fittedparams.keys()[i]].value = x[i]
    x = odesolve(model, t)

    return [    (np.sum((reference['obsAKTPP']           - x['obsAKTPP'])          /14401.0) ** 2),
                (np.sum((reference['obsErbB1_ErbB_P_CE'] - x['obsErbB1_ErbB_P_CE'])/14401.0) ** 2),
                (np.sum((reference['obsERKPP']           - x['obsERKPP'])          /14401.0) ** 2)
           ]

# Magnitude scaling a [0..1] set of points into the +/- 3 orders of magnitude
def scaling(points): return scale.magnitude(points, np.array(fittedparams.values()))

# Too hard this way. Need to batch it up
#v = Varsens(objective, scaling, k, n, True)

sample = Sample(k, n, scaling)
sample.export("varsens/samples/egfr-batch-", ".csv", 190)