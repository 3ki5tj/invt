Overview
=========

This project aims at testing the effectiveness
of the inverse-time formula for the updating magnitude
in adaptive sampling scheme, such as
the Wang-Landau algorithm and metadynamics.



Programs
========

 Program  |  Description
----------|------------------------------------
invt      | Main program
predict   | Predict error


Test programs
-------------

 Program  |  Description
----------|------------------------------------
modetest  | Test program for mode decomposition
langtest  | Test for Langevin equations
outest    | Test for the Ornstein-Uhlenbeck process


Source code
===========


invt
----

File      | Description
----------|-------------------
invt.c    | driver for a 1D model
invtpar.h | handling input parameters
invt.h    | utility routines
invtsamp.h| sampling routines
intq.h    | optimal schedule by integrating the equation of motion
cosmodes.h| eigenmode decomposition for translationally invariant updating schemes
cmvar.h   | computing the autocorrelation integrals from the variance of the bias potential
corr.h    | correlation of eigenmodes
invt.cfg  | sample configuration file
metad.h   | wrapper for using the invt.h, etc. in a realistic simulation
is2.c     | driver for MC simulation


Compilation
===========

The main program is `invt`, to compile it, use
```
make
```
Or manual compilation,
```
gcc invt.c -o invt -lm
```



Running invt
============

To run `invt`, type
```
./invt
```

You can add options, for a list of options, type
```
./invt --help
```

To load parameters from configuration file
(for example, `invt.cfg` in this directory)
```
./invt invt.cfg
```


Computing the error
-------------------

By default the program computes the terminal error
according to the number the schedule
```
alpha(t) =  c / (t + t0).
```

Several key parameters are the following.

### Updating scheme

The default updating scheme is the Wang-Landau one,
which updates a single bin only.
However, the program can update several nearby bins.
To update the nearest neighbors with a relative weight of 0.1 each,
the next nearest neighbors with a relative weight of 0.05 each, type
```
./invt --nb=0.1,0.05
```
This represents a symmetric updating window

  i - 2 | i - 1 |   i   | i + 1 | i + 2
  ------|-------|-------|-------|-------
  0.05  |  0.1  |  0.7  |  0.1  |  0.05

The value at bin i is automatically filled such that
the total is equal to 1.0.
In the case of reflective boundary condition,
the out-of-boundary values are wrapped around, for example,
for i = 0 (first bin), i - 1 means 0, i - 2 means 1, etc.

Note that the window function has to be stable.
That is, the Fourier transform must be nonnegative.
To check if a window function is stable,
run `invt` with `-vv` or `--verbose=2`,
it will print out a table of `gamma_i` in the beginning.
Make sure that all values are nonnegative.

### Sampling scheme

The default sampling scheme is a global Metropolis one,
which is close to the perfect sampling in this case.
Another option is the local (nearest-neighbor) Metropolis sampling.
This option is turned on by `--samp=l`.
The last option is the heatbath sampling (strict perfect sampling),
which is turned on by `--samp=h`.

Example
```
./invt --samp=l
```


### Number of bins

The number of bins can be specified by the `-n` option
```
./invt -n100
```
or equivalently
```
./invt --n=100
```


Computing autocorrelation function
-----------------------------------

It is possible to compute the autocorrelation functions
under a fixed updating magnitude alpha,
determined by the option `--a0`.
For example
```
./invt --a0=0.0001 --corr --fixa
```
Such a run will by default compute the autocorrelation functions
of different updating modes, which are save in the file given "corr.dat".

For these simulations, `ntrials` is automatically set to 1.
So it helps to increase the number of steps by manually
setting `--nsteps=`.
On a 64-bit machine, a value greater than 2e9 is acceptable.

Interesting figures can be observed by changing the sampling scheme
from global Metropolis move (default) to the local one, `--samp=l`.
Then, difference modes will show dramatically different magnitudes.




Options of invt
===============


Regular options
----------------


 o  Sampling scheme

 o  Updating scheme
    1. The default updating scheme is the single-bin (Wang-Landau) scheme.

    2. To explicitly specify a multiple-bin scheme,
    use the option `--nb`
    For example, `--nb=0.15,0.05` means the following updating window
    ```
    i - 2   i - 1   i     i + 1   i + 2
    0.05    0.15    0.6   0.15    0.05
    ```

    3. For a Gaussian updating window, it is entirely determined
    by the width
    ```
    --sig=3
    ```
    means a window with half width being 3.0


Correlation functions
----------------------

 o  Compute correlation functions, `--corr`.

 o  Name of output autocorrelation function,
    specified by `--fncorr=`.
    If this option is explicitly specified,
    `--corr` is unnecessary.

 o  Interval of saving frames for computing autocorrelation functions
    This is determined by the `--nstcorr=`.
    By default (value 0), it is set to `1/a0`,
    where `a0` is the value given by `--alpha0`.
    If this option is explicitly specified,
    and the value is different from 0,
    `--corr` is unnecessary.

 o  Tolerance to truncate autocorrelation functions, `tol`
    On output of the autocorrelation functions
    will be truncated when they drop under `tol` times the peak value.
    Usually, `tol` is positive number, which means that
    we assume the negative values of the autocorrelation functions
    are errors.
    Since the sampling process is Markovian,
    this assumption is usually good.


Running predict
===============


The program `predict` can predict the error
from the analytical result.

Note, the error reported by predict
are the square root values

The program scans over a range of an input parameters
and outputs the error there.
There are three scan modes.

1. c-scan

Scan the proportionality constant c in the schedule
alpha = c / (t + t0),
and outputs the error


2. nb-scan

For the nearest-neighbor updating scheme, scan the
magnitude of the updating the nearest neighbor,
compute the optimal c, and outputs the error.


3. sig-scan

For the Gaussian updating scheme, scan the width of
the Gaussian, compute the optimal c,
and outputs the error.


