#!/usr/bin/env python



'''
Compute error over a range of values of
the width sigma of the Gaussian updating scheme, or
the cutoff K of the bandpass updating scheme
'''



import sys, os, getopt, shutil, re
import zcom



fncfg = None
fnout = None
cmdopt = "--opta"
verbose = 0

# parameters for sig scan
sigscan = 0
sigmin = 1
sigdel = 1
sigmax = 10

# parameters for ok scan
okscan = 0
okmin = 1
okdel = 1
okmax = 20



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS] input.cfg""" % sys.argv[0]

  print """
  Compute the error over a range of width sigma, or cutoff K,

  OPTIONS:

    --sig=sigmin:sigdel:sigmax  set the range of width sigma for the Gaussian updating scheme
    --ok=okmin:okdel:okmax      set the range of the cutoff K for the bandpass (sinc) updating scheme
    -o                          set the output file
    --opt=                      set options to be passed to the command line
    -v                          be verbose
    --verbose=                  set verbosity
    -h, --help                  help
  """
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "hvo:",
        [ "sig=",
          "ok=", "okmax=", "K=",
          "output=", "opt=",
          "prd=", "predict=",
          "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global fncfg, fnout, cmdopt, verbose
  global sigscan, sigmin, sigdel, sigmax
  global okscan, okmin, okdel, okmax

  for o, a in opts:
    if o in ("--sig",):
      sigarr = a.split(':')
      sigmin = float( sigarr[0] )
      sigdel = float( sigarr[1] )
      sigmax = float( sigarr[2] )
      sigscan = 1
    elif o in ("--ok", "--okmax", "--K"):
      okarr = a.split(':')
      okmin = float( okarr[0] )
      okdel = float( okarr[1] )
      okmax = float( okarr[2] )
      okscan = 1
    elif o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("--opt",):
      cmdopt = a
    elif o in ("-o", "--output"):
      fnout = a
    elif o in ("-h", "--help"):
      usage()

  if len(args) == 0:
    usage()
  fncfg = args[0]

  if not fnout:
    fnout = os.path.splitext(fncfg)[0] + "_err.dat"



def geterror(out):
  ''' get the initial and final errors from the output '''

  # construct the pattern
  numpat = "([0-9][\.0-9e+-]+)"
  pat = "average error: .*sqr %s -> %s, norm. sqr %s, stdsqr %s -> %s" % (
      numpat, numpat, numpat, numpat, numpat)

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
  enorm = float( m.group(3) )
  stdei = float( m.group(4) )
  stde  = float( m.group(5) )

  return e, ei, enorm, stde, stdei



def scan():
  global fncfg, fnout, cmdopt, fnout

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
  if sigscan:
    sig = sigmin
    while sig < sigmax + sigdel * 0.01:
      cval += [ sig, ]
      sig += sigdel
    srange = "sig=%s:%s:%s" % (sigmin, sigdel, sigmax)
    varname = "Gaussian sigma"

  else:
    ok = okmin
    while ok < okmax + okdel * 0.01:
      cval += [ ok, ]
      ok += okdel
    srange = "ok=%s:%s:%s" % (okmin, okdel, okmax)
    varname = "bandpass cutoff K"

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
    print "%d: testing for %s value %s..." % (
        i, varname, cval[i])

    c = cval[i]

    if sigscan:
      cmd = "%s --sig=%s" % (cmd0, c)
    else:
      cmd = "%s --okmax=%s" % (cmd0, c)

    ret, out, err = zcom.runcmd(cmd.strip(), capture = True)
    e, ei, enorm, stde, stdei = geterror(out)

    # construct a line of output
    s = "%s\t%s\t%s\t%s\t%s\t%s\n" % (
        c, e, ei, enorm, stde, stdei)
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
  scan()

