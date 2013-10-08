'''
Created on 8 Oct 2013

@author: vimaier

run_lumi_code <description>

'''

import sys
import argparse

from Beam import Beam
from Config import Config
from Luminosity import Luminosity

#===================================================================================================
# parse_args()-function
#===================================================================================================
def parse_args():
    ''' Parses the arguments and returns args '''
    parser = argparse.ArgumentParser()

    parser.add_argument("-n", "--name", dest="config_name", help="Configuration name", default="US2a")

    args = parser.parse_args()

    return args


#===================================================================================================
# main()-function
#===================================================================================================
def main(config_name):
    '''
    :Parameters:
        'config_name': string

    :Return: int
        0 if execution was successful otherwise !=0
    '''
    config = Config(config_name)

    print "Initializing Beam"
    my_beam = Beam(config)
    print ""
    my_beam.printBeamParam()
    print ""

    print "Initializing Luminosity calculations"
    lumi = Luminosity(config, config_name)
    print ""
    lumi.printLumiParam()
    print ""


    lumi.doFill(my_beam)

    #print ""
    #my_beam.printBeamParam()
    #print ""

    return 0


#===================================================================================================
# helper-functions
#===================================================================================================


#===================================================================================================
# main invocation
#===================================================================================================
def run_main():
    args = parse_args()

    return_value = main(config_name = args.config_name)

    sys.exit(return_value)

if __name__ == "__main__":
    run_main()