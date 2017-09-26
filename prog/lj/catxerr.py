#!/usr/bin/env python

''' combine xerr.dat from different directories
Example:
  python catxerr.py drun14/xerr.dat drun14a/xerr.dat drun14b/xerr.dat -o xerr_c.dat
'''

import sys, os

class XErr:
    def __init__(self, fn):
        s = open(fn).readlines()
        self.tag = s[0]
        nsamp = int( s[1][1:].strip() )
        dat = s[2:]
        self.n = n = len(dat)
        self.sx = [[0, 0, 0, 0] for i in range(n)]
        self.sxx = [[0, 0, 0, 0] for i in range(n)]
        self.ref = [0]*n
        self.efref = [0]*n
        self.eiref = [0]*n
        for i in range(n):
            ln = dat[i].strip()
            a = [float(x) for x in ln.split()]
            self.ref[i] = ref = a[9]
            self.efref[i] = a[10]
            self.eiref[i] = a[11]
            for k in range(4):
                a[2*k+2] += ref # average
                self.sx[i][k] = a[2*k+2] * nsamp
                # a[2*k+1] is the variance
                self.sxx[i][k] = (a[2*k+1] + a[2*k+2]**2) * nsamp
        print "loaded %s samples from %s" % (nsamp, fn)
        self.nsamp = nsamp

    def add(self, fn):
        s = open(fn).readlines()
        tag = s[0]
        if tag != self.tag:
            print "file %s has incompatible tagline\n%s%s" % (fn, tag, self.tag)
        nsamp = int( s[1][1:].strip() )
        dat = s[2:]
        n = len(dat)
        if n != self.n:
            print "file %s has incompatible n %d != %d" % (fn, n, self.n)
        for i in range(n):
            ln = dat[i].strip()
            a = [float(x) for x in ln.split()]
            ref = a[9]
            for k in range(4):
                a[2*k+2] += ref
                self.sx[i][k] += a[2*k+2] * nsamp
                self.sxx[i][k] += (a[2*k+1] + a[2*k+2]**2) * nsamp
        self.nsamp += nsamp
        print "added %s samples (now %s) from %s" % (nsamp, self.nsamp, fn)

    def save(self, fn):
        fp = open(fn, "w")
        fp.write(self.tag)
        fp.write("# %d\n" % self.nsamp)
        ave = [0, 0, 0, 0]
        var = [0, 0, 0, 0]
        nsamp = self.nsamp
        for i in range(self.n):
            s = "%d" % (i + 1)
            for k in range(4):
                ave[k] = self.sx[i][k] / nsamp
                var[k] = self.sxx[i][k] / nsamp - ave[k] * ave[k]
                ave[k] -= self.ref[i]
                s += " %g %g" % (var[k], ave[k])
            s += " %g %g %g\n" % (self.ref[i], self.efref[i], self.eiref[i])
            fp.write(s)
        print "saved the combined data to %s with %d samples" % (fn, nsamp)

fns = []
fnout = None
i = 1
while i < len(sys.argv):
    fn = sys.argv[i]
    if fn == "-o":
        i += 1
        fnout = sys.argv[i]
    else:
        fns += [fn,]
    i += 1
print "input:", fns, "output:", fnout


xerr = XErr(fns[0])
for i in range(1, len(fns)):
    xerr.add(fns[i])
if not fnout:
    fn = "xerr.dat"
    if os.path.exists(fn):
        fn = "xerr_c.dat"
xerr.save(fnout)
