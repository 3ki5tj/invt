#!/usr/bin/env python

fn = "verr.log"

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

av = Ave()

for ln in open(fn).xreadlines():
    if ln.startswith("#"): continue
    x = float( ln.split()[1] )
    av.add(x)

print av.cnt, av.getave()
