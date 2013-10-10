import numpy as np
from scipy.integrate import quad, dblquad
from Constants import Constants as cst
import sys

class Luminosity:
  _level_Lumi=[]
  _betatol = 0.0
  _xsec = 0.0
  _days = 0.0
  _efficiency = 0.0
  _turnaround = 0.0
  _step = 0.0
  _max_length = 0.0
  _name = ""
  
  
  def __init__(self,config,name):
    self._level_Lumi = config._level_Lumi
    self._xsec = config._xsec
    self._xsecburn = config._xsecburn
    self._days = config._days
    self._efficiency = config._efficiency
    self._turnaround = config._turnaround
    self._step = config._step
    self._max_length = config._max_length
    self._betatol = config._betatol
    self._name = name
    
  def printLumiParam(self):
    print "Leveling Luminosity: ",self._level_Lumi
    print "Beta tolerance: ",self._betatol
    print "Cross section: ",self._xsec
    print "Cross section for Burn-off: ",self._xsecburn
    print "Days in physics: ",self._days
    print "beam efficiency: ",self._efficiency
    print "Turnaround time: ",self._turnaround
    print "Time steps: ",self._step
    print "Maximum store length: ",self._max_length
    
   #kinematic factor
  def kinfact(self,beam,ip):
    phi = beam._phi[ip]
    k = 2.0*np.cos(phi/2.0)     
    return k
    
  #Luminosity without reduction  
  def getlumi0(self,beam,ip):     
    sig = beam.getsigma(ip) 
    bint = beam._intensity
    frev = beam._frev
    nbunch = beam._nbunch[ip]
    k = self.kinfact(beam,ip)  
    lumi=k*bint**2*frev*nbunch/(8*np.pi*sig[0]*sig[1])*1e-4	    
    return lumi
  
  #Reduction factor
  def getreduc(self,beam,ip):   
    sig = beam.getsigma(ip)
    bet = beam._beta[ip]
    phi = beam._phi[ip]
    xplane = beam._xplane[ip]
    oncc = beam._oncc[ip]
    wcc = beam._wcc
    sigx2 = sig[xplane]**2
    betx = bet[xplane]
    bets = bet[1-xplane]  
    sigs2 = beam._sigs**2    
    #equal beams only so far
    #betx is beta in crossing plane bets is beta in separation plane   
    reduc=1.0
    if(oncc==1):
      reduc = self.getHGXCCfact(phi,sigs2,sigx2,betx,bets,wcc)
    else:
      reduc = self.getHGXfact(sigs2,sigx2,phi,betx,bets)  
    return reduc


  
  #L, L0, R in all IPs
  def getlumi(self,beam):
    lumi=[]
    lumi0=[]
    reduc=[]
    for i in range(beam._nip):
      if i!=1 or not beam._identicalIP01:   # Trick to avoid recalculating for ip =1, which is ip5 equal to ip1, ip=0.  
        r = self.getreduc(beam,i)
        l0 = self.getlumi0(beam,i)
      lumi.append(l0*r)
      lumi0.append(l0)
      reduc.append(r)
    return lumi,lumi0,reduc
    
  def getlumiip(self,beam,ip):
    r = self.getreduc(beam,ip)
    l0 = self.getlumi0(beam,ip)
    return l0*r,l0,r    
    
    
  def integrandnoCC(self,s,sigs2,sigx2,phi,betx,bets):
    const = 1.0/np.sqrt(np.pi*sigs2)
    a = np.sin(phi/2)**2/(sigx2)
    b = np.cos(phi/2)**2/(sigs2)
    return const*np.exp(-(a/(1+s**2/betx**2)+b)*s**2)/(np.sqrt(1+s**2/betx**2)*np.sqrt(1+s**2/bets**2))
    
  def getHGXfact(self,sigs2,sigx2,phi,betx,bets):
    reduc = quad(lambda s : self.integrandnoCC(s,sigs2,sigx2,phi,betx,bets), -np.inf, np.inf)
    return reduc[0]
       
  def integrandCC(self,ct,s,phi,sigs2,sigx2,betx,bets,wcc):
    k=wcc/cst.clight*2*np.pi
    cos2=np.cos(phi/2)**2
    k2=k**2
    sin2=np.sin(phi/2)**2
    sigmax2=sigx2*(1+s**2/betx**2)
    const=1.0/(np.pi*sigs2)
    return const*np.exp(
        -ct**2/sigs2
        -s**2*cos2/sigs2
        -sin2/(4*k2*sigmax2)*(
            2
            +4*k2*s**2
            -np.cos(2*k*(s-ct))
            -np.cos(2*k*(s+ct))
            -8*k*s*np.cos(k*ct)*np.sin(k*s) 
            -4 * np.cos(k*s)**2 * np.sin(k*ct)**2)
        )/np.sqrt(1+s**2/betx**2)/np.sqrt(1+s**2/bets**2)
    

  def getHGXCCfact(self,phi,sigs2,sigx2,betx,bets,wcc):
    reduc= dblquad(self.integrandCC, -np.inf, np.inf, lambda x: -np.inf, lambda x: np.inf, 
       args=(phi,sigs2,sigx2,betx,bets,wcc))
    return reduc[0]
    
  #Luminous region

  def getPeakLumRegion(self, beam, ip):
    sig = beam.getsigma(ip)
    bet = beam._beta[ip]
    phi = beam._phi[ip]
    xplane = beam._xplane[ip]
    oncc = beam._oncc[ip]
    wcc = beam._wcc
    sigx2 = sig[xplane]**2
    betx = bet[xplane]
    bets = bet[1-xplane]
    sigs2 = beam._sigs**2
    if oncc == 1:
      i = quad(self.integrandCC, -np.inf, np.inf, args=(0,phi,sigs2,sigx2,betx,bets,wcc))
      norm = dblquad(self.integrandCC, -np.inf, np.inf, lambda x: -np.inf, lambda x: np.inf,   args=(phi,sigs2,sigx2,betx,bets,wcc))
    else:
      i = [self.integrandnoCC(0,sigs2, sigx2, phi, betx, bets),0]
      norm = quad(self.integrandnoCC, -np.inf, np.inf,    args=(sigs2,sigx2, phi,betx,bets))
    return i[0]/norm[0], 0.5*i[1]/norm[0]


  def getLuminousRegion(self,beam,ip):
    sig = beam.getsigma(ip)
    bet = beam._beta[ip]
    phi = beam._phi[ip]
    xplane = beam._xplane[ip]
    oncc = beam._oncc[ip]
    wcc = beam._wcc
    sigx2 = sig[xplane]**2
    betx = bet[xplane]
    bets = bet[1-xplane]  
    sigs2 = beam._sigs**2 
    if oncc == 1:
      s2integrand= lambda ct,s,phi,sigs2,sigx2,betx,bets,wcc: s**2*self.integrandCC(ct,s,     phi,sigs2,sigx2,betx,bets,wcc)
      i =    dblquad(s2integrand, -np.inf, np.inf, lambda x: -np.inf, lambda x: np.inf, args=(phi,sigs2,sigx2,betx,bets,wcc))
      norm = dblquad(self.integrandCC, -np.inf, np.inf, lambda x: -np.inf, lambda x: np.inf, args=(phi,sigs2,sigx2,betx,bets,wcc))
    else:     
      s2integrand= lambda  s,sigs2,sigx2,phi,betx,bets: s**2*self.integrandnoCC(s,sigs2,sigx2,phi,betx,bets)
      i =    quad(s2integrand,   -np.inf, np.inf,                           args=(sigs2,sigx2,phi,betx,bets))
      norm = quad(self.integrandnoCC, -np.inf, np.inf,                      args=(sigs2,sigx2,phi,betx,bets))
    return np.sqrt(i[0]/norm[0]), 0.5*i[1]/norm[0]

  def getPileUp(self, frev,levLumi15,nbunch15):
    return levLumi15*self._xsec*1.0e-27/nbunch15/frev


  #  beam-beam parameter 
  def getBBparam(self,beam,ip):
    r0 = 1.535e-18
    sig = beam.getsigma(ip)
    phi = beam._phi[ip]
    xplane = beam._xplane[ip]
    oncc = beam._oncc[ip]
    sigx2 = sig[0]**2
    sigy2 = sig[1]**2
    sigs2 = beam._sigs**2  
    sigx = sig[0]
    sigy = sig[1]
    bint=beam._intensity
    bet = beam._beta[ip]
    betax = bet[0]
    betay = bet[1]
    
    if(oncc==1):
      #Approximation assumes perfect crabs
      geomx=1.0
      geomy=1.0
    else:
      #Approximation, neglects hour glass and CC RF curv
      if xplane==0:
        geomx = np.sqrt(1.+(sigs2)/(sigx2)*np.tan(phi/2.)**2)
        geomy = 1.0
      else:
        geomy = np.sqrt(1.+(sigs2)/(sigy2)*np.tan(phi/2.)**2)
        geomx = 1.0
    
    xi = []
    xi.append(bint*r0*betax/(2*np.pi*beam._gamma*sigx*geomx*(sigx*geomx+sigy*geomy)))
    xi.append(bint*r0*betay/(2*np.pi*beam._gamma*sigy*geomy*(sigx*geomx+sigy*geomy)))
    return xi


    
  #Compute integrate lumi - optimum store length  
  def getIntlumi(self,lumi,t):
    l0 = 0.0
    nh = 0.0
    lint = 0.0
    linttmp = 1.0
    ltot = self._days*24*3600*self._efficiency
    f=open("intlumi_"+self._name+".dat","w")
    while linttmp > lint:
      nh=nh+600
      l0 = nh
      lint = linttmp
      linttmp=0.0
      for i in range(len(lumi)-1):
	if t[i] < l0:
	  linttmp =linttmp + (lumi[i]+lumi[i+1])/2.0*(t[i+1]-t[i])
      linttmp = linttmp*ltot/(l0+self._turnaround*3600)
      print >>f, l0/3600.0, linttmp
    f.close()
    return  l0/3600.0, linttmp
    
  def getBetaLevel(self,beam,ip): 
    levlumi = self._level_Lumi[ip]
    beta = beam._beta[ip]
    bmin = beam._betamin[ip]
    xplane = beam._xplane[ip]
    r = abs(bmin[0]-bmin[1])
    lumitot = self.getlumiip(beam,ip)
    
    #if current lumi smaller than levlumi and beta limit reached
    if lumitot[0] < levlumi and beta[1-xplane]<=bmin[1-xplane]:
      return beam._beta[ip],lumitot
      
    #first compute rough estimate (Work pretty well!!! for small time steps avoids going through the loop)
    fact = levlumi/lumitot[0]
    if(r<1e-3):
      beta = beta/fact 
    else:
      beta[1-xplane] = beta[1-xplane]/fact**2
    beam._beta[ip] = beta
    if beta[1-xplane]<bmin[1-xplane]:
      beam._beta[ip] = beam._betamin[ip]
    lumitot = self.getlumiip(beam,ip)
    
    #If good enough move on
    if lumitot[0]/levlumi<1.01 and lumitot[0]/levlumi>0.99:
      return beam._beta[ip],lumitot
      
    #Now look for exact beta with steps betatol
    while lumitot[0] > levlumi:
      if(r<1e-6):
	beta = beta + self._betatol 
      else:
	beta[1-xplane] = beta[1-xplane] + self._betatol
      beam._beta[ip] = beta
      lumitot = self.getlumiip(beam,ip)
      
    #If good enough move on
    if lumitot[0]/levlumi<1.01 and lumitot[0]/levlumi>0.99:
      return beam._beta[ip],lumitot
      
    #Try going other direction      
    while lumitot[0] < levlumi and beta[1-xplane]>bmin[1-xplane]:
      if(r<1e-6):
	beta = beta - self._betatol 
      else:
	beta[1-xplane] = beta[1-xplane] - self._betatol
      beam._beta[ip] = beta
      if beta[1-xplane]<bmin[1-xplane]:
	beam._beta[ip] = beam._betamin[ip]
      lumitot = self.getlumiip(beam,ip)
      
    return beam._beta[ip],lumitot
     
  #Burn-off  
  def getIntensity(self,t,lumi,beam):
    rate = 0.0
    for i in range(beam._nip):
      rate = rate + lumi[i]*self._xsecburn*1.0e-27/beam._nbunch[i]
    beam._intensity = beam._intensity-t*rate
  
  #IBS and radiation damping    
  def IBS_SRgrowth(self,beam):
    taux=((16*beam._tuneb*(beam._epsxn*beam._epsyn)*np.sqrt(beam._kappa)*
	      np.sqrt(beam._kappa+1.)*beam._gamma*beam._sigs*beam._dpp)/
		(2.*cst.clight*(cst.r0**2)*beam._intensity*23))/3600.
    steph = float(self._step)/3600.0
    # Horizontal
    epsx1=beam._epsxn*(1+steph/taux)-2*steph/beam._taux_sr*(beam._epsxn-beam._epsx0)
    # Vertical
    tauy=(1./beam._kappa_c)*taux
    epsy1=beam._epsyn*(1+steph/tauy)-2*steph/beam._tauy_sr*(beam._epsyn-beam._epsy0)
    # Longitudinal
    tauz=(1.+beam._kappa_c)*taux   
    epss=beam._sigs*beam._dpp*(1+steph/tauz)-2*steph/beam._tauz_sr*(beam._sigs*beam._dpp-beam._sigs0*beam._dpp0)
    sigs1 = np.sqrt(epss*beam._sigs/beam._dpp)
    dpp1 = np.sqrt(epss*beam._dpp/beam._sigs)
    beam._epsxn = epsx1
    beam._epsyn = epsy1
    beam._sigs = sigs1
    beam._dpp = dpp1
    
  def doFill(self,beam):
    intensity0=beam._intensity
    epsx0=beam._epsxn
    maxpileupdensity=0.
    maxpeakdensity=0.
    maxpileup=0.
    lumiall = self.getlumi(beam)
    
    print "Virtual Luminosity: "
    print lumiall[0]
    print "Luminosity no reduction: "
    print lumiall[1]
    print "Reduction factor: "
    print lumiall[2]
    
    fout = open('level_'+self._name+'.out','w')

    for i in range(beam._nip):
      fout.write("@ Lp"+str(i)+" = "+str(lumiall[0][i])+"\n")
      fout.write("@ F"+str(i)+" = "+str(lumiall[2][i])+"\n")  


    fout.write("hour"+"\t"+"bint"+"\t"+"\t"+"epsxn"+"\t"+"epsyn"+"\t"+"sigs"+"\t"+"dpp"+"\t")
    for i in range(beam._nip):
      fout.write("betax"+str(i)+"\t"+"betay"+str(i)+"\t"+"lumitot"+str(i)+"\t"+"reduc"+str(i)+"\t")
    fout.write("xitotx"+"\t"+"xitoty"+"\t"+"Lregion"+"PileUp"+"\n")
      
    
    xix=np.zeros(beam._nip)
    xiy=np.zeros(beam._nip)
    lumiip0 =[]
    timeinleveling=-1
    i=0
    #for i in range(0,int(self._max_length),int(self._step)): # Old loop until 15 hours
    lint = 0.0
    linttmp = 0.0000001 # small number to pass the first while evaluation
    ltot = self._days*24*3600*self._efficiency
    xitotx=0.
    xitoty=0.
    xixmax=0.
    xiymax=0.
    while (linttmp > lint and i<=int(self._max_length)):
             
      if i > 0.0:
        if beam._sigs>0.075:
          beam._sigs= beam._sigs*(1.+  beam._lengthleveling)
        if beam._sigs < beam._initsigs and    beam._constantbunchlength:
          beam._sigs = beam._initsigs  # To be done correctly we should go trough the voltage and change dpp
        lint = linttmp
	self.getIntensity(self._step,lumiall,beam)
	self.IBS_SRgrowth(beam)
	
      print "Intensity: ",beam._intensity
      lumiall = []
      reducall = []
      print "Time: ",float(i)/3600.0
      for j in range(beam._nip):
	if j!=1 or not beam._identicalIP01:   # Trick to avoid calculating for IP5 (ip 1) since it is identical to ip=0 (IP1)
          beta,lumitot =  self.getBetaLevel(beam,j)
        xix[j], xiy[j] = self.getBBparam(beam, j)
	if j ==0:
	  print beta
	  print lumitot[0]
	  lumiip0.append(lumitot[0])	
	  xplane = beam._xplane[j]
          #Check if still leveling 
	  if beam._beta[j][1-xplane]<=beam._betamin[j][1-xplane] and timeinleveling < 0:
	    timeinleveling=i/3600.0
          
          LuminousRegion0=self.getLuminousRegion(beam, 0)
          pileup0=self.getPileUp(beam._frev,lumitot[0], beam._nbunch[0])
          pileupdensity=pileup0/LuminousRegion0[0]/1000. # event/mm
          peakdensity=pileup0*self.getPeakLumRegion(beam, 0)[0]/1000.  #event/mm
          if pileupdensity > maxpileupdensity:
            maxpileupdensity = pileupdensity
          if peakdensity > maxpeakdensity:
            maxpeakdensity = peakdensity
          if pileup0 > maxpileup:
            maxpileup=pileup0
          #compute temp int lumi
          if len(lumiip0)>1:
            factnow=ltot/(i+self._turnaround*3600)
            factbefore=ltot/(i-self._step+self._turnaround*3600)
            linttmp = (linttmp/factbefore + (lumiip0[-1]+lumiip0[-2])/2.0*self._step) * factnow
          print "linttmp, lint", linttmp, lint
	lumiall.append(lumitot[0])
	reducall.append(lumitot[2])

      xitotx=sum(xix)
      xitoty=sum(xiy)
      if xitotx>xixmax:
        xixmax = xitotx
      if xitoty>xiymax:
        xiymax = xitoty
      fout.write(str(i/3600.0)+"\t"+str(beam._intensity)+"\t"+str(beam._epsxn)+"\t"+str(beam._epsyn)+"\t"+str(beam._sigs)+"\t"+str(beam._dpp)+"\t")
      for k in range(beam._nip):
	fout.write(str(beam._beta[k][0])+"\t"+str(beam._beta[k][1])+"\t"+str(lumiall[k])+"\t"+str(reducall[k])+"\t")
      fout.write(str(xitotx)+"\t"+str(xitoty)+"\t"+str(LuminousRegion0[0])+"\t"+str(pileup0)+"\t"+str(lint*1e-39)+"\t"+str(peakdensity)+"\t")
      fout.write(str(xix[0])+"\t"+str(xix[1])+"\t"+str(xix[2])+"\t"+str(xiy[0])+"\t"+str(xiy[1])+"\t"+str(xiy[2])+"\t"+  "\n")
                                                                       
      i=i+self._step
    
    if i==self._max_length:
      leveltime=i/3600.
    else:
      leveltime=(i-2*self._step)/3600.
    print leveltime,lint*1.0e-39,timeinleveling, maxpeakdensity
    print >> fout,"@", self._name, str(intensity0/1e11)[:4], str(epsx0*1e6)[:4], beam._beta[0][0],beam._beta[0][1] ,str(lint*1e-39)[0:5], str(leveltime)[0:3], str(timeinleveling)[0:3], str(maxpileup)[0:4], str(maxpeakdensity)[0:4], str(xixmax)[0:5], str(xiymax)[0:5], str(max(lumiip0)*1e-34)[0:5],str(beam._sigs)[0:5]
    
     
 

