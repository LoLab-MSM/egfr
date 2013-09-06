# Load libraries
import numpy as np
from   egfr.erbb_exec import model
import egfr.parameter_dict_A431
import pickle
from   pysb.integrate import odesolve
from   varsens import *
import time
import sys
import os

# Load fitted parameters
with open('egfr/calibration_A431_highEGF_fittedparams_unnorm_splined_RESTART', 'rb') as handle:
    fittedparams = pickle.loads(handle.read())

# Set those fits in the model
for i in range(len(model.parameters)):
    if model.parameters[i].name in fittedparams:
        model.parameters[i].value = fittedparams[model.parameters[i].name]

# Generate Reference
t = np.linspace(0, 7200, 14401)
reference = odesolve(model, t)

# Given a set of parameters, determine outcome
def objective(x):
    for i in range(len(x)):
        model.parameters[fittedparams.keys()[i]].value = x[i]
    x = odesolve(model, t)

    return [    (np.sum((reference['obsAKTPP']           - x['obsAKTPP'])          /14401.0) ** 2),
                (np.sum((reference['obsErbB1_ErbB_P_CE'] - x['obsErbB1_ErbB_P_CE'])/14401.0) ** 2),
                (np.sum((reference['obsERKPP']           - x['obsERKPP'])          /14401.0) ** 2)
           ]

samplefile            = sys.argv[1]
(sampledir,filename)  = os.path.split(samplefile)
outfilename           = os.path.join(sampledir, "objective-"+(filename.split("-")[-1]))

samples = np.loadtxt(samplefile)
result = np.array([objective(sample) for sample in samples])
np.savetxt(outfilename, result)

return 0
