class Config:
  #parameters for class Beam
  _momeV = 0.0
  _intensity = 0.0
  _pmass = 0.0
  _epsxn = 0.0
  _epsyn = 0.0
  _sigs = 0.0
  _dpp = 0.0
  _beta = [] 
  #parameters for class Machine
  _circ = 0.0
  _rho = 0.0
  _frev=0.0
  _av_beta = 0.0
  _dx = 0.0
  _alfmom = 0.0
  _tuneb = 0.0
  _tunes = 0.0
  _nip = 0
  _nbunch = []
  _sepLR = []
  _xplane = []
  _wcc = 0.0
  _oncc=[]
  #parameters for luminosity calculations
  _level_Lumi=[]
  _betatol=0.0
  _xsec = 0.0
  _xsecburn = 0.0
  _days = 0.0
  _efficiency = 0.0
  _turnaround = 0.0
  _step = 0.0
  _max_length = 0.0
  #IBS
  _kappa=0.0
  _kappa_c = 0.0
  #bunch length leveling
  _lengthleveling = 0
  _constantbunchlength=0   # 0 no 1 yes
  
  def __init__(self,name):
    print "Loading configuration ", name
    self.loadcommon()
    exec("self.load"+name+"()")
    self.checkconfig()


  def loadcommon(self):
      #General
      self._momeV = 6500.0e9
      self._pmass = 0.93827231e9
      self._circ = 26658.8832
      self._rho = 3000.0
      self._alfmom = 3.225e-4
      self._tuneb = 64.31
      self._tunes = 0.002342
      self._nip = 3
      self._betatol = 0.006
      self._step = 600
      self._xsec = 85
      self._xsecburn = 100

      self._days = 160
      self._dpp = 1.2e-4  
      self._efficiency = 0.5
      self._turnaround = 3.0
      self._max_length = 15*3600 
      self._wcc = 400.0e6
      self._xplane = [0,0,1]
      self._lengthleveling = 0
      #IBS
      self._kappa=1.0
      self._kappa_c = 0.1
      self._constantbunchlength=0

  def loadUS2a(self):
    #General
    self._nbunch = [2760,2760,2760]
    self._sepLR = [12.0,12.0,12.0]
    self._intensity = 2.2e11
    self._oncc=[1,1,0]
    #Beam
    self._beta = [[0.15,0.15],[0.15,0.15],[3.5,3.5]]
    self._epsxn = 2.5e-6
    self._epsyn = 2.5e-6
    self._sigs = 0.1
    #LUminosity
    self._level_Lumi = [5.1e34,5.1e34,2.0e33]  
  



  def loadUS2b(self):
    #General
    self._nbunch = [2748,2748,2748]
    self._sepLR = [10.0,10.0,10.0]
    self._intensity = 2.2e11
    self._oncc=[1,1,0]
    #Beam
    self._beta = [[0.3,0.075],[0.30,0.075],[3.5,3.5]]
    self._epsxn = 2.5e-6
    self._epsyn = 2.5e-6
    self._sigs = 0.075
    #LUminosity
    self._level_Lumi = [5.1e34,5.1e34,2.0e33]  
    
  def loadUS2NoPic(self):
    #General
    self._nbunch = [2760,2760,2760]
    self._sepLR = [12.0,12.0,12.0]
    self._intensity = 2.2e11
    self._oncc=[0,0,0]
    self._xplane = [1,0,1]
    #Beam
    self._beta = [[0.25,0.8],[0.8,0.25],[3.5,3.5]]
    self._epsxn = 2.5e-6
    self._epsyn = 2.5e-6
    self._sigs = 0.1
    #LUminosity
    self._level_Lumi = [5.1e34,5.1e34,2.0e33]




  def load200MHzPushed(self):
    #General
    self._nbunch = [2760,2760,2760]
    self._sepLR = [15.0,15.0,15.0]
    self._intensity = 2.5e11
    self._oncc=[0,0,0]
    #Beam
    self._beta = [[0.2,1],[1,0.2],[3.5,3.5]]
    self._epsxn = 2.5e-6
    self._epsyn = 2.5e-6
    self._sigs = 0.15
    self._xplane = [1,0,1]
    #LUminosity
    self._level_Lumi = [5.1e34,5.1e34,2.0e33]

  def load200MHz400CC(self):
    self._nbunch = [2760,2760,2760]
    self._sepLR = [12.0,12.0,12.0]
    self._intensity = 2.56e11
    self._oncc=[1,1,0]
    self._beta = [[0.15,0.15],[0.15,0.15],[3.5,3.5]]
    self._epsxn = 3e-6
    self._epsyn = 3e-6
    self._sigs = 0.13
    self._xplane = [1,1,1]
    self._level_Lumi = [5.1e34,5.1e34,2.0e33]
    self._lengthleveling = -0.008 # Relative increment every step

  def load200MHz400(self):
    self._nbunch = [2760,2760,2760]
    self._sepLR = [12.0,12.0,12.0]
    self._intensity = 2.56e11
    self._oncc=[0,0,0]
    self._beta = [[0.20,0.8],[0.8,0.20],[3.5,3.5]]
    self._epsxn = 3e-6
    self._epsyn = 3e-6
    self._sigs = 0.13
    self._xplane = [1,0,1]
    self._level_Lumi = [5.1e34,5.1e34,2.0e33]
    self._lengthleveling = -0.008 # Relative increment every step

  def load200MHz(self):
    #General
    self._nbunch = [2760,2760,2760]
    self._sepLR = [12.0,12.0,12.0]
    self._intensity = 2.56e11
    self._oncc=[0,0,0]
    #Beam
    self._beta = [[0.20,0.8],[0.8,0.20],[3.5,3.5]]
    self._epsxn = 3e-6
    self._epsyn = 3e-6
    self._sigs = 0.15
    self._xplane = [1,0,1]
    #LUminosity
    self._level_Lumi = [5.1e34,5.1e34,2.0e33]

  def load200MHzBCMS(self):
    self._nbunch = [2504,2504,2504]
    self._sepLR = [12.0,12.0,12.0]
    self._intensity = 2.56e11
    self._oncc=[0,0,0]
    self._beta = [[0.20,0.8],[0.8,0.20],[3.5,3.5]]
    self._epsxn = 2.4e-6
    self._epsyn = 2.4e-6
    self._sigs = 0.15
    self._xplane = [1,0,1]
    self._level_Lumi = [5.1e34,5.1e34,2.0e33]

  def load200MHz8b4e(self):
    self._nbunch = [1840,1840,1840]
    self._sepLR = [12.0,12.0,12.0]
    self._intensity = 2.84e11
    self._oncc=[0,0,0]
    self._beta = [[0.20,0.8],[0.8,0.20],[3.5,3.5]]
    self._epsxn = 2.6e-6
    self._epsyn = 2.6e-6
    self._sigs = 0.15
    self._xplane = [1,0,1]
    self._level_Lumi = [4.6e34,4.6e34,1.2e33]




  def loadUS2NoPicMaxi(self):
    #General
    self._nbunch = [2760,2760,2760]
    self._sepLR = [12.0,12.0,12.0]
    self._intensity = 2.5e11
    self._oncc=[0,0,0]
    self._xplane = [1,0,1]
    #Beam
    self._beta = [[0.25,0.8],[0.8,0.25],[3.5,3.5]]
    self._epsxn = 2.5e-6
    self._epsyn = 2.5e-6
    self._sigs = 0.1
    #LUminosity
    self._level_Lumi = [5.1e34,5.1e34,2.0e33]
    

  def load2012(self):
      self._momeV = 4000.0e9
      self._nbunch = [1340,1340,1340]
      self._sepLR = [10.0,10.0,10.0]
      self._intensity = 1.6e11
      self._oncc=[0,0,0]
      self._beta = [[0.6,0.6],[0.6,0.6],[3.5,3.5]]
      self._epsxn = 2e-6
      self._epsyn = 2e-6
      self._sigs = 0.1
      self._level_Lumi = [1e34,1e34,2.0e33]
      self._efficiency = 0.4



  def loadLS1Nom(self):
    #General
    self._nbunch = [2760,2760,2760]
    self._sepLR = [10.0,10.0,10.0]
    self._intensity = 1.1e11
    self._oncc=[0,0,0]
    #Beam
    self._beta = [[0.4,0.4],[0.4,0.4],[3.5,3.5]]
    self._epsxn = 3e-6
    self._epsyn = 3e-6
    self._sigs = 0.1
    #LUminosity
    self._level_Lumi = [1.2e34,1.2e34,2.0e33]
    self._xplane = [1,0,1]
    
  def loadLS18b4e(self):
    #General
    self._nbunch = [1840,1840,1840]
    self._sepLR = [10.0,10.0,10.0]
    self._intensity = 1.9e11
    self._oncc=[0,0,0]
    #Beam
    self._beta = [[1.1,1.1],[1.1,1.1],[3.5,3.5]]
    self._epsxn = 2e-6
    self._epsyn = 2e-6
    self._sigs = 0.1
    #LUminosity
    self._level_Lumi = [1.e34,1.e34,1.0e33]
    self._xplane = [1,0,1]
    self._lengthleveling = 0.002   
    self._max_length = 19*3600
    
  def loadNoPic(self):
    #General
    self._nbunch = [2508,2508,2508]
    self._sepLR = [12.0,12.0,12.0]
    self._intensity = 1.38e11
    self._oncc=[0,0,0]
    #Beam
    self._beta = [[0.2,0.8],[0.80,0.20],[3.5,3.5]]
    self._epsxn = 1.8e-6
    self._epsyn = 1.8e-6
    self._sigs = 0.075
    #LUminosity
    self._level_Lumi = [5.1e34,5.1e34,2.0e33]
    self._xplane = [1,0,1]


  def loadPIC8be4(self):
    self._sepLR = [12.0,12.0,12.0]
    self._intensity = 1.85e11
    self._nbunch = [1840,1840,1840]
    self._oncc=[0,0,0]
    self._beta = [[0.2,0.8],[0.80,0.20],[3.5,3.5]]
    self._epsxn = 2.2e-6
    self._epsyn = 2.2e-6
    self._sigs = 0.1
    self._level_Lumi = [5.1e34,5.1e34,2.0e33]
    self._xplane = [1,0,1]




  def loadPIC1(self):
    #General
    self._nbunch = [2592,2592,2592]
    self._sepLR = [12.0,12.0,12.0]
    self._intensity = 1.38e11
    self._oncc=[0,0,0]
    #Beam
    self._beta = [[0.2,0.4],[0.4,0.2],[3.5,3.5]]
    self._epsxn = 1.8e-6
    self._epsyn = 1.8e-6
    self._sigs = 0.0755
    #LUminosity
    self._level_Lumi = [4.6e34,4.6e34,2.0e33]
    self._xplane = [1,0,1]
    self._constantbunchlength=1
    
  def loadPIC2(self):
    #General
    self._nbunch = [2760,2760,2760]
    self._sepLR = [12.0,12.0,12.0]
    self._intensity = 1.38e11
    self._oncc=[0,0,0]
    #Beam
    self._beta = [[0.20,0.40],[0.40,0.20],[3.5,3.5]]
    self._epsxn = 2.22e-6
    self._epsyn = 2.22e-6
    self._sigs = 0.075
    #LUminosity
    self._level_Lumi = [4.6e34,4.6e34,2.0e33]
    self._xplane = [1,0,1]


  def loadUS1ct(self):
    self._nbunch = [2508,2508,2508]
    self._sepLR = [10.0,10.0,10.0]
    self._intensity = 1.9e11
    self._oncc=[0,0,0]
    self._xplane = [1,0,1]
    self._beta = [[0.2,0.8],[0.8,0.2],[3.5,3.5]]
    self._epsxn = 2.62e-6
    self._epsyn = 2.62e-6
    self._sigs = 0.075
    self._level_Lumi = [4.64e34,4.64e34,2.0e33]

  def loadUS18b4e(self):
    self._nbunch = [1840,1840,1840]
    self._sepLR = [10.0,10.0,10.0]
    self._intensity = 2.85e11
    self._oncc=[0,0,0]
    self._xplane = [1,0,1]
    self._beta = [[0.2,0.8],[0.8,0.2],[3.5,3.5]]
    self._epsxn = 2.62e-6
    self._epsyn = 2.62e-6
    self._sigs = 0.075
    self._level_Lumi = [3.7e34,3.7e34,.8e33]



  def loadUS1a(self):
    self._nbunch = [2748,2748,2748]
    self._sepLR = [10.0,10.0,10.0]
    self._intensity = 1.9e11
    self._oncc=[0,0,0]
    self._xplane = [1,0,1]
    self._beta = [[0.2,0.4],[0.4,0.2],[3.5,3.5]]
    self._epsxn = 2.62e-6
    self._epsyn = 2.62e-6
    self._sigs = 0.0755
    self._level_Lumi = [4.64e34,4.64e34,2.0e33]

  def loadUS1b(self):
    #General
    self._nbunch = [2592,2592,2592]
    self._sepLR = [10.0,10.0,10.0]
    self._intensity = 1.9e11
    self._oncc=[0,0,0]
    #Beam
    self._beta = [[0.2,0.4],[0.4,0.2],[3.5,3.5]]
    self._epsxn = 2.62e-6
    self._epsyn = 2.62e-6
    self._sigs = 0.075
    #LUminosity
    self._level_Lumi = [4.6e34,4.6e34,2.0e33]
    self._xplane = [1,0,1]




  def loadUS2_noX_noCC(self):
    #General
    self._momeV = 6500.0e9
    self._pmass = 0.93827231e9
    self._circ = 26658.8832
    self._rho = 3000.0
    self._alfmom = 3.225e-4
    self._tuneb = 64.31
    self._tunes = 0.002342
    self._nip = 3
    self._nbunch = [2760,2760,2600]
    self._sepLR = [0.0,0.0,0.0]
    self._xplane = [1,0,1]
    self._intensity = 2.2e11
    self._wcc = 400.0e6
    self._oncc=[0,0,0]
    #Beam
    self._beta = [[0.15,0.15],[0.15,0.15],[3.5,3.5]]
    self._epsxn = 2.5e-6
    self._epsyn = 2.5e-6
    self._sigs = 0.075
    self._dpp = 1.2e-4
    #LUminosity
    self._level_Lumi = [5.1e34,5.1e34,2.0e33]
    self._betatol = 0.006
    self._xsec = 85
    self._days = 160
    self._efficiency = 0.5
    self._turnaround = 3.0
    self._step = 600
    self._max_length = 15*3600  
    #IBS
    self._kappa=1.0
    self._kappa_c = 0.1
    self.checkconfig()
    
  def checkconfig(self):
      if(len(self._nbunch)!=self._nip):
	print "ERROR: wrong bunch definition: ",self._nip,self._nbunch
      if(len(self._beta)!=self._nip):
	print "ERROR: wrong beta1 definition: ",self._nip,self._beta
      if(len(self._sepLR)!=self._nip):
	print "ERROR: wrong sepLR definition: ",self._nip,self._sepLR
      if(len(self._xplane)!=self._nip):
	print "ERROR: wrong xplane definition: ",self._nip,self._xplane
      if(len(self._level_Lumi)!=self._nip):
	print "ERROR: wrong Level_Lumi definition: ",self._nip,self._level_Lumi
  
