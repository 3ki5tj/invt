Overview
=========

This project aims at testing the effectiveness
of the inverse-t scheme for the updating magnitude
in adapative sampling scheme, such as
Wang-Landau and metadynamics.



Compilation and running
=======================

The main program is `invt`, to compile it, use
```
make
```
Or manual compilation,
```
gcc invt.c -lm
```


Options
=======

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
-
