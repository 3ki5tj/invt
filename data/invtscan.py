#!/usr/bin/env python



'''
Compute error over a range of values of
the width sigma of the Gaussian updating scheme, or
the cutoff K of the bandpass updating scheme
'''



import sys, os, getopt, shutil, re, math
import zcom



fncfg = None
fnout = None
cmdopt = "--opta"
verbose = 0

# parameters for sig scan
sigscan = 0
sigmin = 1
sigdel = 1
sigmax = 12

# parameters for ok scan
okscan = 0
okmin = 1
okdel = 1
okmax = 20

# parameters for inverse ok scan
ikscan = 0
ikmin = 1
ikdel = 1
ikmax = 12
ikconst = 100 / math.sqrt(2 * math.pi)  # n / sqrt(2 * pi)

# parameters for the initial updating magnitude scan
iascan = 0
iamin = 5e-7
iamax = 2e-2
iadel = 0.2


def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS] input.cfg""" % sys.argv[0]

  print """
  Compute the error over a range of width sigma, or cutoff K

  OPTIONS:

    --sig=sigmin:sigdel:sigmax  set the range of width sigma for the Gaussian updating scheme
    --ok=okmin:okdel:okmax      set the range of the cutoff K for the bandpass (sinc) updating scheme
    --ik=ikmin:ikdel:ikmax      set the range of 1/K for the bandpass (sinc) updating scheme
    --ia=iamin:iadel:iamax      set the range of the initial updating magnitude, increase by a factor of 10^iadel
    --ikC                       set the constant of proportionality of the above scan
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
          "ik=", "invok=", "iK=", "invK=",
          "ikC=", "iKC=",
          "ia=", "inita=",
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
  global ikscan, ikmin, ikdel, ikmax, ikconst
  global iascan, iamin, iadel, iamax

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
    elif o in ("--ik", "--invok", "--iK", "--invK"):
      ikarr = a.split(':')
      ikmin = float( ikarr[0] )
      ikdel = float( ikarr[1] )
      ikmax = float( ikarr[2] )
      ikscan = 1
    elif o in ("--ikC", "--iKC"):
      ikC = float(a)
    elif o in ("--ia", "--inita",):
      iaarr = a.split(':')
      iamin = float( iaarr[0] )
      iadel = float( iaarr[1] )
      iamax = float( iaarr[2] )
      iascan = 1
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
  global fncfg, fnout, cmdopt 

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

  elif okscan:
    ok = okmin
    while ok < okmax + okdel * 0.01:
      cval += [ ok, ]
      ok += okdel
    srange = "ok=%s:%s:%s" % (okmin, okdel, okmax)
    varname = "bandpass cutoff K"

  elif ikscan:
    ik = ikmin
    okold = -1
    while ik < ikmax + ikdel * 0.01:
      ok = int(ikconst / ik + 0.5)
      if ok != okold: # avoid redundancy
        cval += [ ok, ]
      ik += ikdel
      okold = ok
    srange = "ik=%s:%s:%s" % (ikmin, ikdel, ikmax)
    varname = "inverse bandpass cutoff K"

  elif iascan:
    ia = iamin
    iafac = 10.0 ** iadel
    while ia < iamax * iafac:
      cval += [ ia, ]
      ia *= iafac
    srange = "ia=%s:%s:%s" % (iamin, iadel, iamax)
    varname = "initial updating magnitude, alpha(0)"

  else:
    raise

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
    elif okscan or ikscan:
      cmd = "%s --okmax=%s" % (cmd0, c)
    elif iascan:
      cmd = "%s --opta --inita=%s" % (cmd0, c)
    else:
      raise

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

