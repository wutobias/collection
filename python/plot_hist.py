import numpy as np
from scipy  import stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate

k_b  = 1.38064852 * 10**-23
R_g  = 8.3144
Temp = 300

class hist(object):

   def __init__(self, x, y, xlim=None, ylim=None, num_bins=200):

      if xlim == None:

        self.xmin, self.xmax = x.min(), x.max()

      else:

        self.xmin, self.xmax = xlim[0], xlim[1]

      if ylim == None:

        self.ymin, self.ymax = y.min(), y.max()

      else:

        self.ymin, self.ymax = ylim[0], ylim[1]

      self.x           = x
      self.y           = y
      self.num_bins    = num_bins

   def plot2d(self, xlab="X0", ylab="Y0", title="Analysis", name="output.png", inline=False, smooth=True, norm=True,
              weights=None, contour=True, hist_label="Counts", annotate=None, hist_lim=None, log=False, exp=False):

      plt.clf()

      #xnew, ynew = np.mgrid[self.xmin:self.xmax:200j, self.ymin:self.ymax:200j]
      xx, yy                      = np.linspace(self.xmin, self.xmax, self.num_bins),\
                                    np.linspace(self.ymin, self.ymax, self.num_bins)
      hist, xedges, yedges, image = plt.hist2d(self.x, self.y, 
                                               bins=(xx, yy), 
                                               weights=weights, 
                                               normed=norm)
      plt.clf()

      hist = hist.T

      if log and exp:
        raise ValueError("Connot specify both exp and log!")
      if log:
        hist = np.log(hist)
      elif exp:
        hist = np.exp(hist)

      if smooth:
      	#interpolation="spline16"
      	interpolation="nearest"
      else:
      	interpolation="none"

      if hist_lim == None:
      	hist_lim = np.max(hist)

      fig = plt.figure()
      ax  = fig.gca()

      plot2d = plt.imshow(hist, interpolation=interpolation, 
                                origin='lower',
                                extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], 
                                aspect='auto',
                                cmap=cm.coolwarm)

      #cfset = plt.contourf(xedges[:-1], yedges[:-1], hist, levels=np.linspace(0,hist_lim,100))
      cbar  = fig.colorbar(plot2d)
      cbar.ax.set_ylabel(hist_label)

      if contour:
        
        cset = ax.contour( xedges[:-1], yedges[:-1], hist, colors='k')
        # Add the contour line levels to the colorbar
        cbar.add_lines(cset)
        
        if inline:
        
           # Contour-label
           ax.clabel(cset, inline=1, fontsize=10)

      if annotate != None:

      	# annotate must be of type list.
      	# annoate = [ [ xy1_point, xy1_text, text1 ], ...]
      	# xy is x-y-coordinates as array, text is string
      	
      	for xy, xy_text, s in annotate:

      		ax.annotate(xy=xy, xytext=xy_text, s=s,
      			arrowprops=dict(facecolor='black', shrink=0.05))

      ax.set_xlabel(xlab)
      ax.set_ylabel(ylab)
      ax.set_title(title)
      ax.set_xlim((self.xmin, self.xmax))
      ax.set_ylim((self.ymin, self.ymax))

      fig.savefig(name, dpi=1000)
      fig.clf()


class hist_1d(object):

   def __init__(self, x, xlim=None, num_bins=200):

      if xlim == None:

        self.xmin, self.xmax = x.min(), x.max()

      else:

        self.xmin, self.xmax = xlim[0], xlim[1]

      self.x           = x
      self.num_bins    = num_bins
      self.kernel      = stats.gaussian_kde(self.x, bw_method="scott")
      self.xx          = np.linspace(self.xmin, self.xmax, num_bins)
      self.f           = self.kernel.evaluate(self.xx)


   def plot(self, xlab="X0", ylab="Normalized Density", ylim=None, title="Normalized Density Plot", name="output.png", color=None, kde=True, norm=True, weights=None, vertical=None, return_plt=False, label=None):

      fig       = plt.figure()
      ax        = fig.gca()

      if not kde:

        if label != None:

          n, bins, patches = plt.hist(self.x, self.xx, color=color, weights=weights, normed=norm, facecolor='green', alpha=0.5, label=label)

        else:

          n, bins, patches = plt.hist(self.x, self.xx, color=color, weights=weights, normed=norm, facecolor='green', alpha=0.5)

        ax.set_xlim(self.xmin, self.xmax)

      else:

        if label != None:

          ax.plot(self.xx, self.f, color=color, label=label)

        else:

          ax.plot(self.xx, self.f, color=color)

        ax.set_xlim(self.xmin, self.xmax)
        ax.set_ylim(0.0, np.max(self.f)+0.0001)
      
      ax.set_xlabel(r"%s"%xlab)
      ax.set_ylabel(r"%s"%ylab)

      if ylim != None:

        ax.set_ylim(ylim[0], ylim[1])

      if vertical != None:

        ax.axvline(x = vertical, linewidth=2, linestyle="dashed", color="r")

      # Tweak spacing to prevent clipping of ylabel
      fig.subplots_adjust(left=0.15)
      #ax.title(title)

      if return_plt:

        return fig, ax

      fig.savefig(name)
      plt.clf()

# alias
hist_2d = hist