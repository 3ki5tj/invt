/* Generalized ensemble of adaptive umbrella sampling
 * */


"use strict";

/* values for lnzmethod */
//enum { LNZ_WL, LNZ_AVE };
//const char *lnz_names[2] = {"WL", "Ave"};

var SQRT2 = 1.4142135623730951;

function AGE(xcmin, xcmax, delx, sig, lnzmethod,
  c1, alpha0, xmin, xmax, dx, pbc)
{
  var i, xn, n;

  n = Math.floor((xcmax - xcmin) / delx) + 1;
  this.n = n;
  this.lnzmethod = lnzmethod;
  this.alpha0 = alpha0;
  this.alphamm0 = alpha0;
  this.ave = new Array(n);
  this.sig = new Array(n);
  this.c1 = new Array(n);
  this.c2 = new Array(n);
  this.lnz = new Array(n);
  this.mmwl = new Array(n);
  this.cnt = new Array(n);
  this.acc = new Array(n);
  this.hfl = new Array(n);
  for ( i = 0; i < n; i++ ) {
    this.ave[i] = xcmin + i * delx;
    this.sig[i] = sig;
    this.c1[i] = c1;
    this.c2[i] = 0;
    this.lnz[i] = (i - (n - 1)*0.5) * c1 / sig * delx;
    this.mmwl[i] = new MMWL(this.alphamm0);
    this.cnt[i] = 0;
    this.acc[i] = 0;
    this.hfl[i] = 0;
  }
  this.xmin = xmin;
  this.dx = dx;
  this.xn = xn = (xmax - xmin) / this.dx + 1;
  this.xmax = this.xmin + this.dx * xn;
  this.hist = new Array(n);
  this.htot = new Array(xn);
  for ( i = 0; i < n; i++ ) {
    this.hist[i] = new Array(xn);
    for ( var j = 0; j < xn; j++ ) {
      this.hist[i][j] = 0.0;
    }
  }
  for ( i = 0; i < xn; i++ ) {
    this.htot[i] = 0.0;
  }
  this.alphawl = this.alpha0;
  this.t = 0;
  this.t0 = 1;
  this.invt = false;
  this.costab = mkcostab(n, pbc);
}




AGE.prototype.trimv = function(v)
{
  var i, n = this.n, v0;

  for ( i = 0; i < n; i++ ) v0 += v[i];
  v0 /= n;
  for ( i = 0; i < n; i++ ) v[i] -= v0;
}



/* compute difference of the partition function from the average method
 * lnz[j] - lnz[i] */
AGE.prototype.getlnzaveij = function(j, i, corr)
{
  var sigi = this.sig[i], sigj = this.sig[j], dx, dv, dy;

  dx = this.ave[j] - this.ave[i];
  dv = 0.5 * (this.c1[i]/sigi + this.c1[j]/sigj) * dx;
  dv += Math.log(sigj/sigi);
  if ( corr ) {
    dy = (this.c2[j]*SQRT2 - 1)/(sigj*sigj)
       - (this.c2[i]*SQRT2 - 1)/(sigi*sigi);
    dv -= 0.25*dy*(dx*dx/3 + sigi*sigi + sigj*sigj);
    dv += (this.c2[j] - this.c2[i]) / SQRT2;
  }
  return dv;
}



/* compute the partition function from the average method */
AGE.prototype.getlnzave = function()
{
  var i, n = this.n;

  this.lnz[0] = 0;
  for ( i = 1; i < n; i++ )
    this.lnz[i] = this.lnz[i-1] + this.getlnzaveij(i, i-1, 1);
  this.trimv(this.lnz);
}



/* retrieve the updating magnitude */
AGE.prototype.getalpha = function(i)
{
  var alpha = this.invt ? this.n / (this.t + this.t0) : this.alphawl;

  if ( this.lnzmethod == "WL" ) {
    this.alphamm = alpha;
  } else {
    this.alphamm = this.mmwl[i].getalpha();
  }
  if ( this.alphamm > this.alphamm0 ) {
    this.alphamm = this.alphamm0;
  }
  return alpha;
}



/* update the parameters of the bias potential */
AGE.prototype.add = function(i, x, acc)
{
  var ave = this.ave[i], sig = this.sig[i], y1, y2;
  var j;
  var mm = this.mmwl[i];

  var alpha = this.getalpha(i);

  y1 = (x - ave) / sig;
  y2 = (y1 * y1 - 1) / SQRT2;
  this.c1[i] += y1 * this.alphamm;
  this.c2[i] += y2 * this.alphamm;
  mm.add(y1, y2);

  this.t += 1;
  this.cnt[i] += 1;
  this.acc[i] += acc;
  j = (x - this.xmin) / this.dx;
  this.hist[i][j] += 1;

  if ( this.lnzmethod == "WL" ) {
    this.lnz[i] += alpha;
  }
}



/* calculate the fluctuation of histogram modes (from the variance) */
AGE.prototype.calcfl_std = function()
{
  var i, n = this.n, mm, nhh;
  var tot = 0, fl0 = 0, fl1 = 0, fl2 = 0, fl, hi;

  /* compute the total sample size */
  for ( i = 0; i < n; i++ )
    tot += this.mmwl[i].mm[0];
  if ( tot <= 0 ) return 99.0;

  for ( i = 0; i < n; i++ ) {
    mm = this.mmwl[i];
    mm.calcfl(1);
    hi = mm.mm[0]/tot;
    nhh = n * hi * hi;
    fl0 += nhh; /* inter-umbrella */
    fl1 += nhh * mm.fl[1] * mm.fl[1]; /* intra */
    fl2 += nhh * mm.fl[2] * mm.fl[2]; /* intra */
  }
  fl0 -= 1;
  fl = fl0 + fl1 + fl2;
  this.flfr = [fl0/fl, fl1/fl, fl2/fl];
  return Math.sqrt(fl);
}


/* switch a stage */
AGE.prototype.switch = function(magred)
{
  var i, j, n = this.n, xn = this.xn;

  this.alphawl *= magred;
  if ( this.alphawl < this.n/(this.t + this.t0) ) {
    this.invt = true;
    this.t0 = this.t;
  }
  console.log("alpha " + this.alpha.wl + ", " + this.n/this.t + ", "
    + "fluc " + this.hfluc + ", invt " + this.invt);
  this.t = 0;
  for ( i = 0; i < n; i++ ) {
    /* renew the accumulators */
    this.mmwl[i] = new MMWL(this.alphawl);
    this.mmwl[i].invt = this.invt;
    this.mmwl[i].t0 = this.t0;
  }
  for ( i = 0; i < n; i++ )
    for ( j = 0; j < xn; j++ )
      this.hist[i][j] = 0;
}



/* extensive check of histogram fluctuation */
AGE.prototype.wlcheckx = function(fl, magred)
{
  this.hfluc = this.calcfl_std();
  if ( this.lnzmethod == "WL" && !this.invt && this.hfluc < fl ) {
    this.switch(magred);
    return 1;
  }
  return 0;
}




/* transition to a neighboring umbrella */
AGE.prototype.move = function(x, id, local)
{
  var acc = 0, n = this.n;
  var jd, mm;
  var xi, xj, vi, vj, dv, sigi, sigj;

  if ( local ) {
    // jump to a nearest neighbor
    jd = (rand01() < 0.5) ? id - 1 : id + 1;
    if ( jd < 0 || jd >= n ) return 0;
  } else {
    // jump any other umbrella
    jd = ( id + 1 + (int) (rand01() * (n - 1)) ) % n;
  }
  sigi = this.sig[id];
  sigj = this.sig[jd];
  if ( this.lnzmethod == "Ave" ) {
    /* disable transition during WL stage */
    mm = this.mmwl[id];
    if ( mm.getalpha() < 1e-7 || !mm.invt ) return 0;
    /* approximate change of lnz */
    dv = this.getlnzaveij(jd, id, 1);
  } else {
    dv = this.lnz[jd] - this.lnz[id];
  }
  /* compute the acceptance probability */
  xi = (x - this.ave[ id]) / sigi;
  xj = (x - this.ave[ jd]) / sigj;
  vi = this.c1[ id] * xi + this.c2[ id] * (xi * xi - 1) / SQRT2;
  vj = this.c1[ jd] * xj + this.c2[ jd] * (xj * xj - 1) / SQRT2;
  dv += vj - vi;
  this.tacc = ( dv <= 0 || rand01() < Math.exp(-dv) );
  if ( this.tacc ) {
    /* copy parameters to the new umbrella */
    if ( this.cnt[jd] <= 0 && this.lnzmethod == "Ave" ) {
      this.c1[jd] = this.c1[id];
      this.c2[jd] = this.c2[id];
    }
    this.id = jd;
  }
  return this.id;
}

