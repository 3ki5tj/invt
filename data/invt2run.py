#!/usr/bin/env python



''' Compute error over a range of c values '''



import sys, os, getopt, shutil, re
import zcom



fncfg = None
fnout = None
cmdopt = ""
verbose = 0

# scan c instead of lambda = 1 / c
cscan = 1

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
          "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global fncfg, fnout, cmdopt, verbose
  global cmin, cdel, cmax
  global lmin, ldel, lmax

  for o, a in opts:
    if o in ("--c",):
      carr = a.split(':')
      cmin = float( carr[0] )
      cdel = float( carr[1] )
      cmax = float( carr[2] )
      cscan = 1
    elif o in ("--l", "--lambda"):
      larr = a.split(':')
      lmin = float( larr[0] )
      ldel = float( larr[1] )
      lmax = float( larr[2] )
      cscan = 0
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

  ln = out.strip().split("\n")[-3].strip()

  m = re.search("sqr.* ([0-9][\.0-9e+-]+) -> ([0-9][\.0-9e+-]+)", ln)
  if not m:
    print "line [%s]: no error information" % ln
    raise Exception
  err0 = float( m.group(1) )
  err1 = float( m.group(2) )

  return err0, err1



def main():
  global fncfg, fnout, cmdopt, fnout

  progdir = "../prog"
  if not os.path.isdir(progdir):
    progdir = "../" + progdir

  # build the program
  zcom.runcmd("make -C %s" % progdir)

  prog = "invt"

  # make a copy of the program in case
  # it gets modified or deleted later
  shutil.copy("%s/%s" % (progdir, prog), "./%s" % prog)

  cmd0 = "./%s %s %s" % (prog, fncfg, cmdopt)
  cmd0 = cmd0.strip()

  # construct an array of c-values
  cval = []
  if cscan:
    c = cmin
    while c < cmax + cdel * 0.01:
      cval += [ c, ]
      c += cdel
    srange = "c=%s:%s:%s" % (cmin, cdel, cmax)
  else:
    l = lmin
    while l < lmax + ldel * 0.01:
      cval += [ 1.0 / l, ]
      l += ldel
    srange = "lambda=%s:%s:%s" % (lmin, ldel, lmax)

  # print out an information line
  sinfo = "# %s %s %s\n" % (fncfg, cmdopt, srange)
  open(fnout, "a").write(sinfo)

  for i in range(len(cval)):
    print "%d: testing for c-value %s..." % (i, cval[i])

    c = cval[i]

    # override the c value
    cmd = "%s -c %s" % (cmd0, c)
    ret, out, err = zcom.runcmd(cmd.strip(), capture = True)
    err0, err1 = geterror(out)

    # construct a line of output
    s = "%s %s %s\n" % (c, err0, err1)

    # save to the log files
    open(fnout, "a").write(s)



if __name__ == "__main__":
  doargs()
  main()

