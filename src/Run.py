from math import *
from optparse import OptionParser
import sys

from scipy.interpolate import interp1d
from scipy.integrate import quad, dblquad, fixed_quad
from numpy import *
from scipy import stats

from Beam import Beam
from Config import Config
from Luminosity import Luminosity

#TODO: structure imports in the sections python_standard_lib, 3rd_party_libs and own_modules (vimaier)
#TODO: do not 'from x import *'
# Actually all imports except of the own classes and OptionParser are unnecessary.
# Probably you need 'tan' from math but this should be imported in Luminosity.
#--vimaier

#TODO: Do not 'pollute' the global space. Use main() function. See run_lumi_code.py (vimaier)


parser = OptionParser()
parser.add_option("-n", "--name", dest="name",
                  help="Configuration name", default="US2")

(options, args) = parser.parse_args()

name = options.name
config = Config(name)

print "Initializing Beam"
myBeam = Beam(config)
print ""
myBeam.printBeamParam()
print ""

print "Initializing Luminosity calculations"
lumi = Luminosity(config,name)
print ""
lumi.printLumiParam()
print ""


lumi.doFill(myBeam)

#print ""
#myBeam.printBeamParam()
#print ""

