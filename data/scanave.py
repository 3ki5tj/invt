#!/usr/bin/env python


'''
averaging the error from the output files of invtscan.py
'''


import sys, os, getopt, shutil, re, math
import zcom



fnins = []
fnout = None
nsamp = 1000
verbose = 0



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS] input_files""" % sys.argv[0]

  print """
  Averaging the results from a few input files

  OPTIONS:

    -o                          set the output file
    -n                          number of samples in each data file
    -v                          be verbose
    --verbose=                  set verbosity
    -h, --help                  help
  """
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "hvo:n:",
        [
          "output=",
          "nsamp=",
          "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global fnins, fnout, nsamp, verbose

  for o, a in opts:
    if o in ("-n", "--nsamp"):
      nsamp = atoi(a)
    elif o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-o", "--output"):
      fnout = a
    elif o in ("-h", "--help"):
      usage()

  if len(args) == 0:
    usage()
  fnins = args


class Record:
  def __init__(self, arr):
    self.n = len(arr)
    self.arr = [float( arr[i] ) for i in range(self.n)]
    self.w = nsamp

  def update(self, arr):
    w1 = self.w + nsamp

    oldarr = self.arr[:]
    # update the average
    for i in range(1, 4):
      x = float( arr[i] )
      self.arr[i] = (self.arr[i] * self.w + x * nsamp) / w1

    # update the variance
    for i in range(4, 5):
      # old sum x^2
      vara = self.arr[i]**2 * (self.w - 1)
      x2a = vara + oldarr[i-3]**2
      # new sum x^2
      x = float( arr[i-3] ) # new average
      y = float( arr[i] ) # new error
      varb = y * y * (nsamp - 1)
      x2b = varb + x**2
      x2 = (x2a * self.w + x2b * nsamp) / w1
      var = x2 - self.arr[i-3]**2
      #print x2a, x2b, x2, var, self.w, nsamp, w1
      self.arr[i] = math.sqrt( var / (w1 - 1) )
    self.w = w1



def loadf(fn, dic):
  ''' update the dictionary according to fn '''
  lns = open(fn).readlines()
  n = len(lns)
  tag = lns[0]
  for i in range(1, n):
    s = lns[i]
    arr = s.split()
    key = arr[0]
    if dic[i] != 0:
      if math.fabs( float(arr[0]) - dic[i].arr[0] ) > 1e-10:
        print "key mismatch for %s, line %s, %s vs %s" % (fn, i, arr[0], dic[i].arr[0])
        raise Exception

      #arr1 = dic[i].arr[:]
      dic[i].update(arr)
      #print key, dic[i].w, "\n", dic[i].arr, "\n", arr1, "\n", arr
      #raw_input()
    else:
      dic[ i ] = Record(arr)
  return tag



def ave(fnins, fnout):
  if not fnins: return
  # get the number of lines
  n = len(open(fnins[0]).readlines())
  # allocate the dictionary
  mydic = [0]*n
  tag = ""
  for fnin in fnins:
    tag = loadf(fnin, mydic)
  
  s = tag
  for i in range(1, n):
    a = mydic[i].arr
    s += "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
        a[0], a[1], a[2], a[3], a[4], a[5], mydic[i].w)
  if not fnout:
    print s
  else:
    open(fnout, "w").write(s)


if __name__ == "__main__":
  doargs()
  ave(fnins, fnout)
