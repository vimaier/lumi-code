'''
    I need a module docstring. Eclipse and its templates will help me.
'''
#TODO: module description(vimaier)


import numpy as np

from Constants import Constants as cst
#TODO: take care of naming convention: https://twiki.cern.ch/twiki/bin/view/BEABP/PythonStyleGuide
#        Better 'import constants'
#        access through 'constants.Constants.x'
#        if space does matter(e.g. variable will be used in formula):
#            new local variable in small scope(small function): c = constants.Constants.clight'
#-- vimaier

class Beam:
    ''' I need a description. '''
    #TODO: docstring for class(vimaier)
    #TODO: description(vimaier)
    _momeV = 0.0
    _epsxn = 0.0
    _epsyn = 0.0
    _sigs = 0.0
    _dpp = 0.0
    _pmass = 0.0
    _intensity = 0.0
    _gamma = 0.0
    _betarel = 0.0
    _beta = []
    _betamin = []
    _sepLR = []
    _phi = []
    _circ = 0.0
    _rho = 0.0
    _frev = 0.0
    _av_beta = 0.0
    _dx = 0.0
    _alfmom = 0.0
    _tuneb = 0.0
    _tunes = 0.0
    _nip = 0
    _nbunch = []
    _xplane = []
    _wcc = 0.0
    _oncc = []
    _taux_sr = 0.0
    _tauy_sr = 0.0
    _tauz_sr = 0.0
    _dpp0 = 0.0
    _sigs0 = 0.0
    _epsx0 = 0.0
    _epsy0 = 0.0
    _kappa = 0.0
    _kappa_c = 0.0
    #TODO: the attributes are static attributes by convention. Definition of object member in constructor.
    #     See: http://stackoverflow.com/questions/1697273/differences-between-static-and-instance-variables-in-python-do-they-even-exist
    #          http://stackoverflow.com/questions/68645/static-class-variables-in-python

    def __init__(self,config):
        self._momeV = config._momeV
        self._intensity = config._intensity
        self._pmass = config._pmass
        self._gamma = config._momeV/config._pmass+1.0 #new
        self._betarel = np.sqrt(1.0-1.0/(self._gamma**2.0)) #new
        self._epsxn = config._epsxn
        self._epsyn = config._epsyn
        self._sigs = config._sigs
        self._dpp = config._dpp
        self._beta = np.array(config._beta)
        self._betamin = np.array(config._beta) #new
        self._sepLR = np.array(config._sepLR)
        self._circ=config._circ
        self._rho=config._rho
        self._alfmom = config._alfmom
        self._tuneb = config._tuneb
        self._tunes = config._tunes
        self._nip = config._nip
        self._nbunch = config._nbunch
        self._frev = cst.clight/config._circ #new
        self._dx = config._alfmom*config._circ/(2*np.pi) #new
        self._av_beta = config._circ/(config._tuneb*2*np.pi) #new
        self._xplane = config._xplane
        self._wcc = config._wcc
        self._oncc = config._oncc
        self._identicalIP01= (self._beta[0][0] == self._beta[1][0])*(
                         self._beta[0][1] == self._beta[1][1])*(
                         self._oncc[0] == self._oncc[1])*(
                           self._sepLR[0] == self._sepLR[1])*(
                           self._nbunch[0] == self._nbunch[1]) # True * True * True... is really strange. 'and' operator is more intuitiv.
        #TODO: do not access to protected members(config._x). Write getters
        #TODO: a lot of variables copied from Config. That is not necessary
        #    Better in __init__: self.__config = config
        #--vimaier

        #TODO: extract into private function: def __init_phi() (vimaier
        for i in range(self._nip):
            sigp = np.sqrt(self._epsyn/(self._gamma*self._betarel)/self._beta[i][self._xplane[i]])
            self._phi.append(sigp*self._sepLR[i])

        #TODO: extract into private function: def __init_variables_for_radiation_damping() (vimaier)
        #for radiation damping
        dEsr = cst.e**2*self._betarel**3*self._gamma**4/(3*cst.eps0*self._rho)/cst.e/self._momeV
        self._taux_sr = 2.0/(dEsr*self._frev*3600.0)
        self._tauy_sr = 2.0/(dEsr*self._frev*3600.0)
        self._tauz_sr = 1.0/(dEsr*self._frev*3600.0)

        #TODO: extract into private function(until self._epsy0) (vimaier)
        I1 = self._dx*2*np.pi
        I2 = 2*np.pi/self._rho
        I3 = 2*np.pi/self._rho**2
        I4 = self._circ*self._alfmom/self._rho**2
        I5 = 1/self._av_beta*self._dx**2/self._rho**2*2*np.pi

        cq = 55/(32*np.sqrt(3.0))*cst.hbar/(cst.m*cst.clight)

        self._dpp0 = np.sqrt(cq*self._gamma**2*I3/(2*I2))
        self._sigs0 = self._alfmom*cst.clight/(2*np.pi*self._tunes*cst.clight/self._circ)*self._dpp0
        self._epsx0= cq*self._gamma**2*I5/I2*self._betarel*self._gamma
        self._epsy0= 13.0/55.0*cq/I2*self._av_beta/self._rho**2*2*np.pi*self._betarel*self._gamma

        #TODO: extract into private function: def __init_variables_for_ibs(); what is IBS? (vimaier)
        #for IBS
        self._kappa = config._kappa
        self._kappa_c = config._kappa_c


    def getsigma(self,ip):#TODO: naming convention: get_sigma (vimaier)
        betx = self._beta[ip][0] #TODO: Difficult to remember structure of lists in lists. Better: New datastructure --> new class?
                                #could be then something like: betx = self.get_beta().get_x_for_ip(ip_str)
                                #--vimaier
        bety = self._beta[ip][1]
        sigx = np.sqrt(self._epsxn*betx/(self._gamma*self._betarel))
        sigy = np.sqrt(self._epsyn*bety/(self._gamma*self._betarel))
        return sigx,sigy


    def printBeamParam(self):#TODO: naming convention: print_beam_param (vimaier)
        print "Beam parameters:"
        print "Momentum: ",self._momeV
        print "gamma:",self._gamma
        print "betarel:",self._betarel
        print "Emittances: ",self._epsxn,self._epsyn
        print "Beta at IP: ",self._beta
        print "Bunch length:",self._sigs
        print "dp/p",self._dpp
        print "Rest mass: ",self._pmass
        print "Intensity:",self._intensity
        print ""
        print "Machine parameters:"
        print "Circumference: ",self._circ
        print "Curvature: ", self._rho
        print "Revolution frequency: ",self._frev
        print "Average beta: ",self._av_beta
        print "Average dispersion: ",self._dx
        print "Momentum compaction: ",self._alfmom
        print "Betatron tune: ",self._tuneb
        print "Synchrotron tune: ",self._tunes
        print "Number of IPs: ",self._nip
        print "Colliding pairs/IP: ",self._nbunch
        print "LR separation: ",self._sepLR
        print "Crossing angle: ",self._phi
        print "Crossing plane (H=0,V=1): ", self._xplane
        print "CC frequency: ",self._wcc
        print "CC ON:",self._oncc
        print ""
        print "Radiation damping and IBS:"
        print "Radiation damping times: ",self._taux_sr,self._tauy_sr,self._tauz_sr
        print "IBS kappas: ",self._kappa,self._kappa_c

