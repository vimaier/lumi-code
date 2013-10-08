from scipy.interpolate import interp1d
from scipy.integrate import quad, dblquad, fixed_quad
from numpy import *
from scipy import stats
from math import *
from optparse import OptionParser
import sys
from Beam import Beam
from Config import Config
from Luminosity import Luminosity


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

