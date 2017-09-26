#!/usr/bin/env python

import sys

fn = "verr.log"
if len(sys.argv) > 1:
    fn = sys.argv[1]

class Ave:
    def __init__(self):
        self.cnt = 0
        self.sx  = 0
        self.sxx = 0

    def add(self, x):
        self.cnt += 1
        self.sx += x
        self.sxx += x * x

    def getave(self):
        return self.sx/(self.cnt+1e-300)

avf = Ave()
avi = Ave()

i = 0
xiref = xfref = 0
for ln in open(fn).xreadlines():
    i += 1
    if ln.startswith("#"):
        if i == 2:
            arr = ln.split()
            xiref = float( arr[-2] )
            xfref = float( arr[-1] )
        continue
    arr = ln.split()
    xf = float( arr[1] )
    xi = float( arr[2] )
    avf.add(xf)
    avi.add(xi)

print "%s: count %s, init. %g(%g) final %g(%g)" % (fn,
        avf.cnt, avi.getave(), xiref, avf.getave(), xfref)
