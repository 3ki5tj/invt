### Notes

* Try "--rho=0.1", default density is 0.8.
* Remember to use "--opta" for the optimal schedule for the Gaussian updating scheme.


### List of truncation error

Data for bin width (0.01).
Number of bins is n = 250 (rho = 0.8), n = 500 (rho = 0.1).

 rho  | kc | sig  |trunc. error
------|----|------|--------
 0.1  | 20 |      | 1.37e-6
 0.1  | 15 |      | 3.55e-6
 0.1  | 10 |      | 2.00e-5
 0.1  |    | 0.05 |
 0.1  |    | 0.1  | ~1e-6
 0.8  | 20 |      | 5.18e-7
 0.8  | 15 |      | 1.42e-6
 0.8  | 10 |      | 5.36e-6
 0.8  |    | 0.05 | ~5e-7
 0.8  |    | 0.1  | ~1e-4

### Refine the bias potential

The following simulations only need to run for 1-2 seconds and be stopped.

For the rho = 0.8 simulation,
```
cp run0/vf.dat vref.dat
make lj && ./lj --rho=0.8 --kc=100
cp vtrunc.dat dr0.01/rho0.8/vref.dat
```

For the rho = 0.1 simulation,
```
cp run4/vf.dat vref.dat
make lj && ./lj --rho=0.1 --kc=200
cp vtrunc.dat dr0.01/rho0.1/vref.dat
```

### Compute the gamma value

```
./lj --try=0 --gamnsteps=100000000 --gam=varv
```

### Compute the reference value, `vref.dat`

```
./lj --try=1 --nsteps=100000000
cp vf.dat vref.dat
```

### Comparison efficiency


Wang-Landau:
```
./lj --nsteps=10000000
```
Bandpass
```
./lj --nsteps=10000000 --kc=15
```
Gaussian, inverse-time
```
./lj --nsteps=10000000 --sig=0.1
```
Gaussian, optimal
```
./lj --nsteps=10000000 --sig=0.1 --opta --fngamma=gamma.dat --gam=load
```



arun1-8: rho = 0.8
```
make -C .. && ../lj --gam=load --fngamma=gamma.dat
make -C .. && ../lj --gam=load --fngamma=gamma.dat --kc=20
make -C .. && ../lj --gam=load --fngamma=gamma.dat --sig=0.1
make -C .. && ../lj --gam=load --fngamma=gamma.dat --sig=0.1 --opta
make -C .. && ../lj --gam=load --fngamma=gamma.dat --sig=0.2
make -C .. && ../lj --gam=load --fngamma=gamma.dat --sig=0.2 --opta
make -C .. && ../lj --gam=load --fngamma=gamma.dat --sig=0.5
make -C .. && ../lj --gam=load --fngamma=gamma.dat --sig=0.5 --opta
```

brun1-8: rho = 0.1
```
make -C .. && ../lj --rho=0.1 --gam=load --fngamma=gamma.dat
make -C .. && ../lj --rho=0.1 --gam=load --fngamma=gamma.dat --kc=20
make -C .. && ../lj --rho=0.1 --gam=load --fngamma=gamma.dat --sig=0.1
make -C .. && ../lj --rho=0.1 --gam=load --fngamma=gamma.dat --sig=0.1 --opta
make -C .. && ../lj --rho=0.1 --gam=load --fngamma=gamma.dat --sig=0.2
make -C .. && ../lj --rho=0.1 --gam=load --fngamma=gamma.dat --sig=0.2 --opta
make -C .. && ../lj --rho=0.1 --gam=load --fngamma=gamma.dat --sig=0.5
make -C .. && ../lj --rho=0.1 --gam=load --fngamma=gamma.dat --sig=0.5 --opta
```

