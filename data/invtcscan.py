#!/usr/bin/env python



'''
Compute error over a range of c values
for the inverse-time schedule, where

  alpha(t) = c / (t + t0),

with t0 assuming the default value
'''



import sys, os, getopt, shutil, re
import zcom



fncfg = None
fnout = None
fnprd = None
cmdopt = ""
verbose = 0

# scan lambda = 1 / c instead of c
lscan = 0

# parameters for c scan
cmin = 0.5
cdel = 0.1
cmax = 1.5

# parameters for lambda scan
lmin = 0.5
ldel = 0.1
lmax = 1.5



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS] input.cfg""" % sys.argv[0]

  print """
  Compute the error over a c range

  OPTIONS:

    --c=cmin:dc:cmax        set the c range in a c-scan
    --l=lmin:dl:lmax        set the lambda range in a lambda = 1/c scan
    -o                      set the output file
    --prd=                  set the prediction file
    --opt=                  set options to be passed to the command line
    -v                      be verbose
    --verbose=              set verbosity
    -h, --help              help
  """
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "hvo:",
        [ "c=", "crange=",
          "l=", "lrange=", "lambda=",
          "output=", "opt=",
          "prd=", "predict=",
          "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global fncfg, fnout, fnprd, cmdopt, verbose
  global lscan
  global cmin, cdel, cmax
  global lmin, ldel, lmax

  for o, a in opts:
    if o in ("--c",):
      carr = a.split(':')
      cmin = float( carr[0] )
      cdel = float( carr[1] )
      cmax = float( carr[2] )
      lscan = 0
    elif o in ("--l", "--lambda"):
      larr = a.split(':')
      lmin = float( larr[0] )
      ldel = float( larr[1] )
      lmax = float( larr[2] )
      lscan = 1
    elif o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("--opt",):
      cmdopt = a
    elif o in ("-o", "--output"):
      fnout = a
    elif o in ("--prd", "--predict"):
      fnprd = a
    elif o in ("-h", "--help"):
      usage()

  if len(args) == 0:
    usage()
  fncfg = args[0]

  if not fnout:
    fnout = os.path.splitext(fncfg)[0] + "_err.dat"

  if not fnprd:
    fnprd = os.path.splitext(fncfg)[0] + "_prd.dat"



def geterror(out):
  ''' get the initial and final errors from the output '''

  # construct the pattern
  numpat = "([0-9][\.0-9e+-]+)"
  pat = "average error: .*sqr %s -> %s,.*stdsqr %s -> %s" % (
      numpat, numpat, numpat, numpat)

  s = out.strip().split("\n")

  m = None
  i = -1
  while i >= -10:
    ln = s[i].strip()

    m = re.search(pat, ln)
    if m: break

    i -= 1

  if not m:
    print "cannot find error information"
    print '\n'.join( s[-10:] )
    exit(1)

  ei    = float( m.group(1) )
  e     = float( m.group(2) )
  stdei = float( m.group(3) )
  stde  = float( m.group(4) )

  return e, ei, stde, stdei



def c_scan():
  global fncfg, fnout, fnprd, cmdopt, fnout

  progdir = "../prog"
  if not os.path.isdir(progdir):
    progdir = "../" + progdir

  # build the program
  zcom.runcmd("make -C %s" % progdir)

  prog = "invt"

  try:
    # make a copy of the program in case
    # it gets modified or deleted later
    shutil.copy("%s/%s" % (progdir, prog), "./%s" % prog)
  except:
    pass

  cmd0 = "./%s %s %s" % (prog, fncfg, cmdopt)
  cmd0 = cmd0.strip()

  # create a table of c-values
  cval = []
  if not lscan:
    # linear c values
    c = cmin
    while c < cmax + cdel * 0.01:
      cval += [ c, ]
      c += cdel
    srange = "c=%s:%s:%s" % (cmin, cdel, cmax)

    # construct the larger range for the analytical prediction
    if cmin > cdel:
      cmin1 = cmin - cdel
    else:
      cmin1 = 0
    cmax1 = cmax + cdel
  else:
    # linear l values, c = 1/l
    l = lmin
    while l < lmax + ldel * 0.01:
      cval += [ 1.0 / l, ]
      l += ldel
    srange = "lambda=%s:%s:%s" % (lmin, ldel, lmax)

    # construct the larger range for the analytical prediction
    cmin1 = 1 / (lmax + ldel)
    if lmin > ldel:
      cmax1 = 1 / (lmin - ldel)
    else:
      cmax1 = 1 / lmin

  dc1 = 0.01
  if cmin1 <= 0:
    cmin1 = dc1
  # generate the prediction result
  os.system("%s/predict %s %s --cmin=%s --cdel=%s --cmax=%s > %s"
      % (progdir, fncfg, cmdopt, cmin1, dc1, cmax1, fnprd) )

  # we save data to a temporary file
  # and overwrite the actual file only at the end
  if os.path.exists(fnout):
    fnouttmp = fnout + ".tmp"
  else:
    fnouttmp = fnout

  # text of the output
  txt = ""

  # print out an information line
  sinfo = "# %s %s %s\n" % (fncfg, cmdopt, srange)
  txt += sinfo

  # go over the predefined c-values
  for i in range(len(cval)):
    print "%d: testing for c-value %s..." % (i, cval[i])

    c = cval[i]

    # override the c value
    cmd = "%s -c %s" % (cmd0, c)
    ret, out, err = zcom.runcmd(cmd.strip(), capture = True)
    e, ei, stde, stdei = geterror(out)

    # construct a line of output
    s = "%s\t%s\t%s\t%s\t%s\n" % (
        c, e, ei, stde, stdei)
    txt += s

    # save to the output file
    open(fnouttmp, "w").write(txt)

  # write the entire text at the end
  # we do not trust the content of fnouttmp,
  # which is only for backup, and may have been changed
  open(fnout, "w").write(txt)

  # remove the temporary file
  if fnouttmp != fnout:
    os.remove(fnouttmp)



if __name__ == "__main__":
  doargs()
  c_scan()

