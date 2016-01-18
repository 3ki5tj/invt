Overview
=========

This project aims at testing the effectiveness
of the inverse-t scheme for the updating magnitude
in adaptive sampling scheme, such as
Wang-Landau and metadynamics.



Compilation and running
=======================

The main program is `invt`, to compile it, use
```
make
```
Or manual compilation,
```
gcc invt.c -o invt -lm
```

Then to run the program, type
```
./invt
```


Computing autocorrelation function
-----------------------------------

It is possible to compute the autocorrelation functions
under a fixed updating magnitude alpha,
determined by the option `--a0`.
For example
```
./invt --a0=0.0001 --fixa
```
Such a run will by default compute the autocorrelation functions
of different updating modes, which are save in
the file given "corr.dat".

For these simulations, `ntrials` is automatically set to 1.
So it helps to increase the number of steps by manually
setting `--nsteps=`.
On a 64-bit machine, a value greater than 2e9 is acceptable.

Interesting figures can be observed by changing the sampling scheme
from global Metropolis move (default) to the local one, `--samp=l`.
Then, difference modes will show dramatically different magnitudes.




Options
=======


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


 o  Name of output autocorrelation function,
    specified by `--fncorr=`.

 o  Interval of saving frames for computing autocorrelation functions
    This is determined by the `--nstcorr=`.
    By default, it is set to `1/a0`,
    where `a0` is the value given by `--alpha0`.

 o  Tolerance to truncate autocorrelation functions, `tol`
    On output of the autocorrelation functions
    will be truncated when they drop under `tol` times the peak value.
    Usually, `tol` is positive number, which means that
    we assume the negative values of the autocorrelation functions
    are errors.
    Since the sampling process is Markovian,
    this assumption is usually good.

