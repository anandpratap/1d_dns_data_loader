import os
from collections import OrderedDict
import numpy as np
from scipy.interpolate import interp1d

_current_file_path = os.path.dirname(os.path.abspath(__file__))

_dns_solution_list = np.array([180, 395, 550, 950, 2000, 4200, 5200])
_dnssol_keys = ["y", "y+", "u+", "ub+", "vb+", "wb+", "-om_z+", "om_xb+", "om_yb+", "om_zb+", "uvb+", "uwb+", "vwb+", "prb+", "psb+", "pstob+", "pb"]
_kbal_keys = ["y", "y+", "dissip", "prodoc", "p-strain", "p-diff", "t-diff", "v-diff", "bal", "tp-kbal"]

# wilcox data is not dimensionalized in the file
_wilcox_keys = ["y", "unknown_1", "unknown_2", "u+", "uvb+", "k+", "dissip", "produc"]

class DNSDataLoader(object):
    def __init__(self, Retau, y = None):
        self.Retau = Retau
        self.dnssol = _dns_solution_list[abs(_dns_solution_list - self.Retau).argmin()]
        print "Selecting nearest Re = %d as proxy for Re = %.2f" %(self.dnssol, self.Retau)
        if y is not None:
            self.y = y
            self.if_interpolate = True
            print "will interpolate."
            assert abs(y[-1] - 0.5) < 1e-13, "y should be between 0.0 and 0.5"
            assert abs(y[0] - 0.0) < 1e-13, "y should be between 0.0 and 0.5" 
        self.data = OrderedDict()
        self.load_dnssol()
        try:
            self.load_kbal()
        except:
            print "Warning: kbal not found!"
        self.load_wilcox()
        if self.if_interpolate:
            self.interpolate()

    def load_dnssol(self):
        dnssol = np.loadtxt(os.path.join(_current_file_path, "data/Re_%d/DNSsol.dat"%self.dnssol))
        nvar = dnssol.shape[1]
        for i in range(nvar):
            fac = 1.0
            if i == 0:
                fac = 0.5
            self.data[_dnssol_keys[i]] = dnssol[:,i]*fac
        self.data['k+'] = 0.5*(self.data["ub+"]**2+self.data["vb+"]**2+self.data["wb+"]**2)

    def load_kbal(self):
        kbal = np.loadtxt(os.path.join(_current_file_path, "data/Re_%d/kbal.dat"%self.dnssol))
        nvar = kbal.shape[1]
        for i in range(2, nvar):
            self.data[_kbal_keys[i]] = kbal[:,i]

    def load_wilcox(self):
        wilcox = np.loadtxt(os.path.join(_current_file_path, "data/Re_%d/wilcox_%d.dat"%(self.dnssol, self.dnssol)))
        nvar = wilcox.shape[1]
        nu = 1e-4
        utau = self.Retau*nu/0.5;
        for i in range(nvar):
            fac = 1.0
            if i == 0:
                fac = 0.5
            elif i == 3:
                fac = 1
            elif i == 4:
                fac = -1
            self.data["wilcox_" + _wilcox_keys[i]] = wilcox[:,i]*fac
        self.data["wilcox_y+"] = wilcox[:,0]*utau/nu

    def interpolate(self):
        print "Interpolating..."
        ytarget = self.y
        ysource = self.data['y']
        wilcox_ysource = self.data['wilcox_y']
        for key in self.data.keys():
            if key not in ["y", "wilcox_y"]:
                if key[0:6] == "wilcox":
                    f = interp1d(wilcox_ysource, self.data[key])
                else:
                    f = interp1d(ysource, self.data[key])
                self.data[key] = f(ytarget)
        self.data['y'] = ytarget
        self.data['wilcox_y'] = ytarget
        
if __name__ == "__main__":
    for Re in _dns_solution_list:
        y = np.linspace(0.0, 0.5, 201);
        loader = DNSDataLoader(Re, y)
        data = loader.data
        print data.keys()
        from matplotlib.pyplot import *
        figure(1)
        title("- DNS, -- Wilcox")
        c = semilogx(data['y+'], data['u+'], label="Re = %.2f"%Re)
        semilogx(data['y+'], data['wilcox_u+'], "--", color=c[0].get_color())
        legend(loc="best")
        xlabel('y')
        ylabel('u+')
        figure(2)
        title("- DNS, -- Wilcox")
        semilogx(data['y+'], data['uvb+'], color=c[0].get_color(), label="Re = %.2f"%Re)
        semilogx(data['y+'], data['wilcox_uvb+'], "--", color=c[0].get_color())
        legend(loc="best")
        xlabel('y')
        ylabel('uv+')
        figure(3)
        title("- DNS, -- Wilcox")
        semilogx(data['y+'], data['k+'],  color=c[0].get_color(), label="Re = %.2f"%Re)
        semilogx(data['y+'], data['wilcox_k+'], "--", color=c[0].get_color())
        legend(loc="best")
        xlabel('y')
        ylabel('k+')
    show()
        
