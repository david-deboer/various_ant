import math
import matplotlib.pyplot as plt
import scipy.special as spec
import scipy.integrate as integ
import numpy as np


print """antPerf.py is a series of functions that determine antenna parameters.
             ants.py has a wrapper class that calls these functions."""
# See Balanis pg 627

def angle2fD(angle,units='degrees'):
    """assumes the full illumination beamwidth"""
    if units=='degrees':
        angle*=(math.pi/180.0)
    X = 4.0*math.tan(angle/4.0)
    fD = 1.0/X
    return fD

def fD2angle(fD,units='degrees'):
    X = 1.0/(4.0*fD)
    angle = 4.0*math.atan(X)
    if units=='degrees':
        angle*=(180.0/math.pi)
    return angle

def antennaMismatch(s11range=[-40.0,0.0,1.0]):
    s11 = np.arange(s11range[0],s11range[1],s11range[2])
    mag21 = np.power(10.0,s11/20.0)
    vswr = (1.0 + mag21)/(1.0 - mag21)
    eff21 = 1.0 - mag21**2.0
    effdB = 10.0*np.log10((1.0-eff21)/10.0)

    plt.figure('EfficiencyLinear')
    plt.plot(s11,eff21)
    plt.xlabel('S11 [dB]')
    plt.ylabel('Efficiency')
    
    plt.figure('EfficiencydB')
    plt.plot(s11,effdB)
    plt.xlabel('S11 [dB]')
    plt.ylabel('Efficiency [dB]')
    
    plt.figure('VSWR')
    plt.plot(s11,vswr)
    plt.xlabel('S11 [dB]')
    plt.ylabel('VSWR')


def freespaceTaper(fD):
    T = -20.0*math.log10(1.0 + (1.0/(4.0*fD))**2)
    return T
def feedTaper(freq, fD, FBW, dB_at_bw, feed_type):
    th = fD2angle(fD,'degrees')/2.0
    T = feedPattern(freq, th, FBW, dB_at_bw, feed_type)
    T = 10.0*math.log10(T)
    return T

def feedPattern(freq, angle, FBW, dB_at_bw, feed_type):
    if feed_type is None:
        print 'No feed_type specified'
        return
    if ':' in feed_type:
        feed_freq = float(feed_type.split(':')[1])
        feed_exponent = float(feed_freq.split(',')[1])
        feed_freq = float(feed_freq.split(',')[0])
        feed_type = feed_type.split(':')[0]
        print 'Not implemented'
    if type(angle) is not np.array:
        angle = np.array(angle,dtype=float)
    P = None
    if feed_type.lower() == 'gaussian':  # beam is full width, level is where that beam is in dB (e.g. 3)
        A = 4.0*dB_at_bw/(10.*math.log10(math.e))
        P = np.exp(-A*(angle/FBW)**2)
    elif feed_type.lower()[0:4] == 'file':
        fname = feed_type.split(':')[1]
        #print 'Opening feed file:  '+fname+' (assuming csv with header line)'
        fpat = np.loadtxt(fname,delimiter=',',skiprows=1)
        avecols = range(1,len(fpat[0]))
        avepat = []
        for i in range(len(fpat)):
            avepat.append(np.average(np.power(10.0,fpat[i,avecols]/10.0)))
        mP = max(avepat)
        theta = fpat[:,0]
        P = np.interp(angle,theta,avepat)
        P = P/mP
    else:
        print feed_type+' not supported'
    return P
def illuminationFactor(r,fD):
    # assumes r is normalized to D (r -> r/(D/2))
    # using negative r values is to allow plotFeed in ants.py to work correctly
    X = 4.0*fD
    if type(r) is float and r<0.0:
        r = X*math.tan(abs(r)*math.pi/180.0/2.0)
        print 'DOUBLE CHECK THIS IN illuminationFactor'
    if type(r) is not np.array:
        r = np.array(r,dtype=float)
    P = 1.0/(1.0 + (r/X)**2)**2
    return P


def compute_fD_given_taper(taper, feedFBW, dB_at_bw, feed_type, firstGuess=0.35):
    """Computes f/D for given taper"""
    # Need to get f/D based on total taper, including free space taper.
    #    get first guess
    print "I HAVEN'T CHANGED/CHECKED THIS YET"
    print 'Computing for ',taper,':  feed is ',feedFBW,dB_at_bw
    print 'Only does gaussian'
    llin = pow(10.0,dB_at_bw/10.0)
    X = math.sqrt( 40.0*math.log(llin)/math.log(10.0) )
    taper = float(taper)
    fD = firstGuess
    err = 1.0
    while err>0.0001:
        edge = taper - freespaceTaper(fD)
        theta0 = feedFBW*(edge/X)
        fD = angle2fD(theta0,'degrees')
        tinc = edge + freespaceTaper(fD)
        err = abs(tinc - taper)/taper
    print "f/D=%.4f (target_taper=%.4f, actual_taper=%.4f, feed_edge=%.4f, freespace_taper=%.4f)" % (fD,taper,tinc,edge,freespaceTaper(fD))
    return fD

def ruze(frequency,rms):
    if type(frequency) is list:
        frequency = np.array(frequency)
    if type(rms) is list:
        rms = np.array(rms)
    if type(frequency) is int:
        frequency = float(frequency)
    if type(rms) is int:
        rms = float(rms)
    lambd = 30./frequency
    if type(frequency) is float and type(rms) is float:
        print 'wavelength = '+str(lambd)+' cm'
        print 'rms = '+str(rms)+' cm'
    sig = 4.0*math.pi*rms/lambd
    eff = np.exp(-sig**2)
    return eff

def illuminationEff(freq,fD, FFBW, dB_at_bw, feed_type):
    """Computes illumination efficiency - Balanis p627.  Illumination taper already included so just use feed pattern."""
    theta0 = fD2angle(fD,units='degrees')
    dtt = 0.1
    theta = np.arange(0.0,theta0/2.0+dtt,dtt)
    g = feedPattern(freq, theta, FFBW, dB_at_bw, feed_type)
    theta = theta*math.pi/180.0
    kern = np.sqrt(g)*np.tan(theta/2.0)
    num = 32.0*(fD*integ.trapz(kern,dx=dtt*math.pi/180.0))**2
    kern = g*np.sin(theta)
    den = integ.trapz(kern,dx=dtt*math.pi/180.0)

    n_ill = num/den    
    return n_ill


def spilloverEff(freq,fD, FFBW, dB_at_bw, feed_type):
    """Computes spillover efficiency  - Balanis p627.  Illumination taper already included so just use feed pattern."""
    theta0 = fD2angle(fD,units='degrees')
    tt = 0.0
    dtt = 0.1
    theta = np.arange(0.0,180.0+dtt,dtt)
    g = feedPattern(freq, theta, FFBW, dB_at_bw, feed_type)
    theta = theta*math.pi/180.0

    # integrate over main beam
    gmb = np.where(theta < (theta0/2.0)*math.pi/180.0)
    kern = g[gmb]*np.sin(theta[gmb])    
    num = integ.trapz(kern,dx=dtt*math.pi/180.0)
    # integrate over full beam
    kern = g*np.sin(theta)
    den = integ.trapz(kern,dx=dtt*math.pi/180.0)
    
    n_spill = num/den
    return n_spill 

def blockageEff(freq, Dpri, dsub, n_ill):
    """Computes blockage efficiency"""
    n_block = (1.0 - (dsub/Dpri)**2/n_ill)**2
    #################print 'blockageEff ==> UPDATE USING KILDAL!!!!!!!'
    return n_block

def test_feedPattern(freq, fD, D, FFBW, dB_at_bw, feed_type):
    r = []
    t = []
    p = []
    h = []
    rint = 0.0
    while rint<=1.0:
        r.append(rint)
        tt = (2.0*180.0/math.pi)*math.atan(rint/(4.0*fD))
        t.append(tt)
        g = feedPattern(tt, FFBW, dB_at_bw, feed_type)
        p.append(g*illuminationFactor(rint,fD))
        h.append(g)
        rint+=0.001
    plt.figure(20)
    plt.plot(r,p)
    plt.figure(21)
    plt.plot(t,p)
    plt.plot(t,h)
    h = np.array(h)
    ih0 = np.where(h<=0.0)
    h[ih0] = 1e-60
    gdb = 10.0*np.log10(h)
    plt.figure(1011)
    plt.plot(t,gdb)

def getTaper(t,p,e):
    """Calculate the feed taper given basically everything:
             t:  angles
             p:  pattern
             e:  edge angle"""
    if type(t) is not np.array:
        t = np.array(t)
    ctr_i = np.where(t==0.0)[0][0]
    tap_i = np.where(t>e)[0][0]
    taper = p[ctr_i] - p[tap_i]
    return taper

def calcBW(t,g,bw):
    """compute FWHM given angle and gain at bw dB down"""
    if type(g) is not np.array:
        g = np.array(g)
    if bw>0.0:
        bw*=-1.0
    maxg = np.max(g)
    g = g-maxg
    #plt.figure('testbw')
    
    imax = np.argmax(g)
    mb = np.where(g>3.0*bw)
    # ...left
    gg = g[mb[0][0]:imax]
    tt = t[mb[0][0]:imax]
    #plt.plot(tt,gg)
    ha1 = np.interp(bw,gg,tt)
    # ...right
    gg = np.flipud(g[imax:mb[0][-1]])
    tt = np.flipud(t[imax:mb[0][-1]])
    ha2 = np.interp(bw,gg,tt)
    #plt.plot(tt,gg)
    # ...full
    FWHM = abs(ha1) + abs(ha2)
    #print 'bw:  ',ha1,ha2
    return FWHM
    
def beamPattern(freq, fD, D, efficiency, FFBW, dB_at_bw, feed_type, defects={},lw=2,plotbp=None,plotgp=102,plot_label_prefix='',plot_color=None):
    """This is the one from Baars - modified for feed pattern + free space rather than using alpha
        D and wavelength must be in same units
        efficiency:  in number, not percentage
        feedFBW:
        BWlev: (default 3.0)
        feedType:  (default 'gaussian')
        defects:  dictionary with error entries.  Allowed:
            'ruze_rms':
            'ruze_corr':
            'block':
            'defocus':
        lw:  linewidth of plotline"""

    wavelength = 300.0/freq
    # compute diffraction pattern
    dtheta = 0.1
    T = np.arange(-90.0,90.0+dtheta,dtheta)
    G = []
    F = fD*D
    dr = 0.001
    r = np.arange(0.0,1.0+dr,dr)
    tt = (2.0*180.0/math.pi)*np.arctan(r/(4.0*fD))
    g = np.sqrt(feedPattern(freq, tt, FFBW, dB_at_bw=dB_at_bw, feed_type=feed_type)*illuminationFactor(r,fD))
    for theta in T:
        u = (math.pi*D/wavelength)*math.sin(theta*math.pi/180.0)
        kern=[]
        for ii,rint in enumerate(r):
            kern.append(g[ii]*spec.jn(0,u*rint)*rint)
        fu = integ.trapz(kern,dx=dr)*math.pi*(D**2)/2.0
        G.append( 10.0*math.log10(fu**2.) )
    G = np.array(G)
    
    # compute error pattern(s) -- see Baars 86-90 (ruze/block)
    if 'ruze_rms' in defects.keys():
        if 'ruze_corr' in defects.keys():
            C = defects['ruze_corr']
        else:
            C = D/25.0  #assume fairly small correlation length
        sigma = 4.0*math.pi*defects['ruze_rms']/wavelength
        if sigma > 1.0:
            ferr = "Doesn't do anything yet"
            
    # compute normalized pattern and FWHM
    bp = G - max(G)
    FWHM = calcBW(T,bp,-3.0)
    taper = freespaceTaper(fD) + feedTaper(freq, fD, FFBW, dB_at_bw=dB_at_bw, feed_type=feed_type)
    #print 'FWHM (f/D=%.2f, taper=%.2f) = %.4f' % (fD,taper,FWHM)

    # plot beam pattern
    if plotbp is not None:
        plt.figure(plotbp)
        s = '%s%.1f-m:  %.1f$^o$' % (plot_label_prefix,D,FWHM)
        if plot_color is not None:
            plt.plot(T,bp,color=plot_color,label=s)
        else:
            plt.plot(T,bp,label=s)
        plt.grid()

    # gain pattern:  compute, plot and write
    Do = 4.0*np.pi*efficiency*(np.pi*D**2.0/4.0)/(wavelength**2)
    Do = 10.0*np.log10(Do)
    gp = bp + Do
    if plotgp is not None:
        plt.figure(plotgp)
        s = '%s:%.1f:  %.1f$^o$' % (plot_label_prefix,D,FWHM)
        if plot_color is not None:
            plt.plot(T,gp,color=plot_color,label=s,lw=lw)
        else:
            plt.plot(T,gp,label=s,lw=lw)
        plt.legend()
        plt.grid()
    bpfn = "beamPattern%.2f" % (freq)
    bpfn = bpfn.replace('.','_') + '.dat'
    print "Writing ",bpfn
    fp = open(bpfn,'w')
    for i,v in enumerate(T):
        s = '%.1f\t%f\n' % (v,gp[i])
        fp.write(s)
    fp.close()

    return FWHM

def integrate_beam(t, P, tmax, normVal=1.0, pattern_type='t_2pi_dB'):
    """pattern_type is whether theta goes from -pi to +pi (t_2pi)
                    versus 0 to pi (t_pi)
       DO OTHER CHECKS/OPTIONS ETC - BUT FOR NOW JUST GET ANSWER!!!"""
    shape_check = np.shape(P)
    step_t = 361.0/shape_check[0]
    step_p = 179.0/shape_check[1]
    #print 'steps:  ',step_t,step_p,shape_check

    if 'dB' in pattern_type:
        Pn = np.power(10.0,(P-normVal)/10.0)
    else:
        Pn = P/normVal
    #plt.figure('linear')
    #for i in range(shape_check[1]):
    #    plt.plot(Pn[:,i])

    Omega = 0.0
    for i in range(shape_check[0]-1):
        if abs(t[i]) > tmax:
            continue
        for j in range(shape_check[1]):
            Omega+= Pn[i,j]*math.fabs(math.sin(t[i]*math.pi/180.0))
    Omega*=((step_t*math.pi/180.0)*(step_p*math.pi/180.0))
    return Omega
