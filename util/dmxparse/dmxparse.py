#! /usr/bin/env python
import sys, struct, math, numpy

class tempo_cov:

    def __init__(self,filename='matrix.tmp'):
        self.read_matrix(filename)

    def reset(self,npar):
        self.npar = npar
        self.params = [''] * npar
        self.errs = [0.0] * npar
        self.cov = [[]] * npar

    def add_row(self,raw):
        stuff = struct.unpack("=ii5sdd" + "d"*self.npar, raw)
        idx = stuff[1] - 1
        self.params[idx] = stuff[2]
        self.errs[idx] = stuff[4]
        self.cov[idx] = list(stuff[5:])

    def get_cov(self,par1,par2):
        idx1 = self.params.index(par1)
        idx2 = self.params.index(par2)
        return self.cov[idx1][idx2]

    def read_matrix(self,filename='matrix.tmp'):
        """Read a tempo matrix.tmp file."""
        f = open(filename,'r')
        run = True
        first = True
        nbhdr = 4 + 4 + 5 + 8 + 8
        while run:
            try:
                l1 = struct.unpack("=i", f.read(4))[0]
                ncov = (l1 - nbhdr)/8
                if first:
                    self.reset(ncov)
                    first = False
                self.add_row(f.read(l1))
                l2 = struct.unpack("=i", f.read(4))[0]
                # TODO could put more error handling here
            except:
                run = False

class dmxparse:
    def __init__(self,parfile):
        self.val = {}
        self.err = {}
        self.epoch = {}
        self.r1 = {}
        self.r2 = {}
        self.f1 = {}
        self.f2 = {}
        for l in open(parfile).readlines():
            if not l.startswith('DMX'): continue
            fields = l.split()
            key = fields[0]
            val = fields[1].replace('D','e')
            if len(fields)==4: 
                flag = int(fields[2])
                err = fields[3].replace('D','e')
            pfx=None
            idx=None
            (pfx,idx,newkey) = (None,None,None)
            if '_' in key: 
                (pfx, idx) = key.split('_')
                newkey = "DX%03d" % int(idx)
            if l.startswith('DMX_') and flag==1:
                self.val[newkey] = float(val)
                self.err[newkey] = float(err)
            if l.startswith('DMXEP_'): self.epoch[newkey] = float(val)
            if l.startswith('DMXR1_'): self.r1[newkey] = float(val)
            if l.startswith('DMXR2_'): self.r2[newkey] = float(val)
            if l.startswith('DMXF1_'): self.f1[newkey] = float(val)
            if l.startswith('DMXF2_'): self.f2[newkey] = float(val)

    def fix_errs(self,tmp_cov):
        self.verr = {}
        n = len(self.err)
        cc = numpy.zeros((n,n))
        k = sorted(self.err.keys())
        for i in range(n):
            for j in range(n):
                cc[i,j] = tmp_cov.get_cov(k[i],k[j]) \
                        * self.err[k[i]] * self.err[k[j]]

        # Find error in mean DM
        self.mean = numpy.mean(self.val.values())
        self.mean_err = numpy.sqrt(cc.sum())/float(n)

        # Do the correction for varying DM
        m = numpy.identity(n) - numpy.ones((n,n))/float(n)
        cc = numpy.dot(numpy.dot(m,cc),m)
        for i in range(n):
            self.verr[k[i]] = math.sqrt(cc[i,i])


if __name__ == '__main__':

    if len(sys.argv)<2:
        print >>sys.stderr, "Usage: dmxparse.py par_file (cov_matrix_file)"
        print >>sys.stderr, ""
        print >>sys.stderr, "Note, if matrix file name is not given, 'matrix.tmp' is assumed."
        sys.exit(0)

    parfile = sys.argv[1]
    if len(sys.argv)>2: matrixfile = sys.argv[2]
    else: matrixfile = 'matrix.tmp'

    d = dmxparse(parfile)
    c = tempo_cov(matrixfile)

    # Check that DMX's match up
    ndmx1 = ['DX' in p for p in c.params].count(True)
    ndmx2 = len(d.err)
    if ndmx1 != ndmx2:
        print >>sys.stderr, "Error: DMX entries do not match up"
        print >>sys.stderr, "  parfile NDMX = %d, matrix NDMX = %d" % (ndmx1, ndmx2)
        sys.exit(1)

    d.fix_errs(c)
    print "# Mean DMX value = %+.6e" % d.mean
    print "# Uncertainty in average DM = %.5e" % d.mean_err
    print "# Columns: DMXEP DMX_value DMX_var_err DMXR1 DMXR2 DMXF1 DMXF2 DMX_bin"
    for k in sorted(d.val.keys()):
        #print d.epoch[k], d.val[k], d.err[k], d.verr[k], k
        print "%.4f %+.7e %.3e %.4f %.4f %7.2f %7.2f %s" % (d.epoch[k], 
                d.val[k]-d.mean, d.verr[k], 
                d.r1[k], d.r2[k], d.f1[k], d.f2[k], k)
