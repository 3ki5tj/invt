
cp ../../prog/lj/brun4/alpha.dat rho0.1/alpha_sig0.1.dat
cp ../../prog/lj/brun6/alpha.dat rho0.1/alpha_sig0.2.dat
cp ../../prog/lj/brun8/alpha.dat rho0.1/alpha_sig0.5.dat

cp ../../prog/lj/drun6/alpha.dat rho0.1/alpha_sig0.2_t1e8.dat
cp ../../prog/lj/drun6/verr.log  rho0.1/verr_sig0.2_t1e8.log

cp ../../prog/lj/brun1/xerr.dat rho0.1/xerr_sbin.dat  
cp ../../prog/lj/brun2/xerr.dat rho0.1/xerr_kc20.dat
cp ../../prog/lj/brun3/xerr.dat rho0.1/xerr_sig0.1_invt.dat
cp ../../prog/lj/brun4/xerr.dat rho0.1/xerr_sig0.1_opt.dat
cp ../../prog/lj/brun5/xerr.dat rho0.1/xerr_sig0.2_invt.dat
cp ../../prog/lj/brun6/xerr.dat rho0.1/xerr_sig0.2_opt.dat
cp ../../prog/lj/brun7/xerr.dat rho0.1/xerr_sig0.5_invt.dat
cp ../../prog/lj/brun8/xerr.dat rho0.1/xerr_sig0.5_opt.dat

cp ../../prog/lj/drun5/xerr.dat rho0.1/xerr_sig0.2_t1e8_invt.dat
cp ../../prog/lj/drun6/xerr.dat rho0.1/xerr_sig0.2_t1e8_opt.dat

cp ../../prog/lj/arun4/alpha.dat rho0.8/alpha_sig0.1.dat
cp ../../prog/lj/arun6/alpha.dat rho0.8/alpha_sig0.2.dat
cp ../../prog/lj/arun8/alpha.dat rho0.8/alpha_sig0.5.dat

cp ../../prog/lj/arun1/xerr.dat rho0.8/xerr_sbin.dat  
cp ../../prog/lj/arun2/xerr.dat rho0.8/xerr_kc20.dat
cp ../../prog/lj/arun3/xerr.dat rho0.8/xerr_sig0.1_invt.dat
cp ../../prog/lj/arun4/xerr.dat rho0.8/xerr_sig0.1_opt.dat
cp ../../prog/lj/arun5/xerr.dat rho0.8/xerr_sig0.2_invt.dat
cp ../../prog/lj/arun6/xerr.dat rho0.8/xerr_sig0.2_opt.dat
cp ../../prog/lj/arun7/xerr.dat rho0.8/xerr_sig0.5_invt.dat
cp ../../prog/lj/arun8/xerr.dat rho0.8/xerr_sig0.5_opt.dat

make -B -C ../../doc/fig lj_alpha.pdf lj_xerr.pdf
