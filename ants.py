import math
import matplotlib.pyplot as plt
import numpy as np
import antPerf as ap
from matplotlib import colors as mcolors

#colors = ['b','r','k','g','c','m','y','b','r']
colors = mcolors.cnames.keys()

class Ants:
    def __init__(self,D=14.0,fD=0.32,feed_type='gaussian',feed_fbw=None,dB_at_bw=None):
        """
        paraboloidal antenna (primary optics) antenna parameters
        This just gets the input information and calls the various modules in antPerf.
        """
        self.D = D
        self.fD = fD
        self.f = fD*D
        if feed_fbw is None and dB_at_bw is None:
            dB_at_bw = 10.
            feed_fbw = ap.fD2angle(fD)
        self.feed_fbw = feed_fbw
        self.dB_at_bw = dB_at_bw
        self.feed_type = feed_type
        self.set_optics(D,None,fD,feed_fbw,dB_at_bw,feed_type,False)
        self.set_defects()

    def ant_par(self,freq=None,plotgp='GainPattern',plot_label_prefix='',plot_color='k'):
        self.feedtaper = ap.feedTaper(freq, self.fD, self.feed_fbw, self.dB_at_bw, self.feed_type)
        self.taper = self.feedtaper + self.freetaper
        self.theta0 = ap.fD2angle(self.fD,units='degrees')
        self.n_ill = ap.illuminationEff(freq,self.fD,self.feed_fbw,self.dB_at_bw,self.feed_type)
        self.n_spill = ap.spilloverEff(freq,self.fD,self.feed_fbw,self.dB_at_bw,self.feed_type)
        self.n_block = ap.blockageEff(freq,self.D,self.defects['block'],self.n_ill)
        self.n_total = self.n_ill*self.n_spill*self.n_block
        self.FWHM = ap.beamPattern(freq,self.fD,self.D,self.n_total,self.feed_fbw,self.dB_at_bw,
                                   self.feed_type,defects=self.defects,
                                   plotgp=plotgp,plot_label_prefix=plot_label_prefix,plot_color=plot_color)
        print '========================ANTPAR=========================='
        print 'Frequency:  %.2f MHz' % (freq)
        print 'Feed %s:  %.2f deg at %.0f dB' % (self.feed_type,self.feed_fbw,self.dB_at_bw)
        print 'Full primary angle (%.2f dB taper) = %.4f' % (self.taper,self.theta0)
        print 'f = %.2f, D=%.2f, f/D = %.4f' % (self.f,self.D,self.fD)
        print 'illumination efficiency = %f' % (self.n_ill)
        print 'spillover efficiency = %f' % (self.n_spill)
        print 'blockage efficiency = %f' % (self.n_block)
        print 'total efficiency = %f' % (self.n_total)
        print 'FWHM = %f' % (self.FWHM)
        print '========================================================'

    def set_optics(self,D=None,f=None,fD=None,feed_fbw=None,dB_at_bw=None,feed_type=None,fit_ant_to_feed=False):
        """Set parameters for feed.
                  feed_type:  only gaussian right now
                  feed_fbw:  full beamwidth of feed at
                  dB_at_bw
                  fD:  f/D
                  if one of ..."""

        #Check/set feed parameters
        if fit_ant_to_feed:
            fD = ae.compute_fD_given_taper(fit_ant_to_feed,feed_fbw,dB_at_bw,feed_type)
            if D is None:
                D = f/fD
            else:
                f = fD*D
        self.freetaper = ap.freespaceTaper(fD)
        
        #Check/set primary parameters
        prim_pars = [D,f,fD]
        if prim_pars.count(None)>1:
            print 'Need to specify 2/3:  D,f,f/D'
            return
        elif prim_pars.count(None)==1:
            if D is None:
                D = f/fD
            elif f is None:
                f = D*fD
            elif fD is None:
                fD = f/D
        else:
            if abs(fD - f/D)/fD > .01:
                print "you specified all f, D, and f/D, but they don't match"
                print "reset f/D to match f and D"
                fD = f/D
        self.fD = fD
        self.D = D
        self.f = f

    def set_defects(self,sub=1.5,ruze=0.0,defocus=0.0):
        self.defects = {'ruze_rms':ruze,'block':sub,'defocus':defocus}

    def plot_feed(self,feedlabel=None,plot_color=None):
        PdB = []
        gdB = []
        xstop = ap.fD2angle(self.fD,units='degrees')
        xstart = -1.0*xstop
        xstep = 0.2
        degrees = np.arange(xstart,xstop+xstep,xstep)
        r = 4.0*self.fD*np.tan(degrees*math.pi/180.0/2.0)
        P = ap.feedPattern(degrees,self.feed_fbw,self.dB_at_bw,self.feed_type)
        g = P*ap.illuminationFactor(r,self.fD)
        PdB = 10.0*np.log10(P)
        gdB = 10.0*np.log10(g)
        if plot_color:
            plt.plot(degrees,PdB,plot_color,label=feedlabel)
            plt.plot(degrees,gdB,plot_color,linestyle='--')
        else:
            plt.plot(degrees,PdB)
            plt.plot(degrees,gdB,'--')
        if feedlabel:
            plt.legend()
        plt.title('Feed Pattern')
        plt.xlabel('Deg')
        plt.ylabel('dB')
        axval = plt.axis()
        plt.grid(True)
        #plt.axis(xmin=-100.0,xmax=100.0,ymin=axval[2],ymax=axval[3])
