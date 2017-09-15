/* Generalized ensemble of adaptive umbrella sampling
 * */


"use strict";

var SQRT2 = 1.4142135623730951;

function AGE(xcmin, xcmax, delx, sig,
  c1, alpha0, xmin, xmax, dx)
{
  var i, xn, n;

  n = Math.floor((xcmax - xcmin) / delx) + 1;
  this.n = n;
  this.alpha0 = alpha0;
  this.ave = new Array(n);
  this.sig = new Array(n);
  this.c0 = new Array(n);
  this.c1 = new Array(n);
  this.c2 = new Array(n);
  this.mm = new Array(n);
  this.cnt = new Array(n);
  this.acc = new Array(n);
  this.hfl = new Array(n);
  this.delx = delx;
  this.bc = c1/sig;
  for ( i = 0; i < n; i++ ) {
    this.ave[i] = xcmin + i * delx;
    this.sig[i] = sig;
    this.c0[i] = (i - (n - 1)*0.5) * c1 / sig * delx;
    this.c1[i] = c1;
    this.c2[i] = SQRT2*0.5;
    this.mm[i] = [0, 0, 0];
    this.cnt[i] = 0;
    this.acc[i] = 0;
    this.hfl[i] = 0;
  }
  this.flfr = [0, 0, 0];
  this.xmin = xmin;
  this.dx = dx;
  this.xn = xn = (xmax - xmin) / this.dx + 1;
  this.xmax = this.xmin + this.dx * xn;
  this.hist = new Array(n);
  for ( i = 0; i < n; i++ ) {
    this.hist[i] = new Array(xn);
    for ( var j = 0; j < xn; j++ ) {
      this.hist[i][j] = 0.0;
    }
  }
  this.htot = new Array(xn);
  this.lng = new Array(xn);
  for ( i = 0; i < xn; i++ ) {
    this.htot[i] = 0.0;
    this.lng[i] = 0.0;
  }
  this.alphawl = this.alpha0;
  this.t = 0;
  this.t0 = 1;
  this.invt = false;
}



/* retrieve the updating magnitude */
AGE.prototype.getalpha = function()
{
  var alpha = this.invt ? this.n / (this.t + this.t0) : this.alphawl;
  return alpha;
}



/* update the parameters of the bias potential */
AGE.prototype.add = function(i, x, acc)
{
  var ave = this.ave[i], sig = this.sig[i];
  var alpha = this.getalpha();

  var y1 = (x - ave) / sig;
  var y2 = (y1 * y1 - 1) / SQRT2;
  this.c0[i] += alpha;
  this.c1[i] += y1 * alpha;
  this.c2[i] += y2 * alpha;
  this.mm[i][0] += 1;
  this.mm[i][1] += y1;
  this.mm[i][2] += y2;

  this.t += 1;
  this.cnt[i] += 1;
  this.acc[i] += acc;
  var j = Math.floor( (x - this.xmin) / this.dx );
  this.hist[i][j] += 1;
}



/* calculate the fluctuation of histogram modes (from the variance) */
AGE.prototype.calcfl_rms = function()
{
  var i, k, n = this.n;
  var tot = 0, fl = [0, 0, 0], flsum = 0, s;

  for ( i = 0; i < n; i++ ) {
    for ( k = 0; k < 3; k++ )
      fl[k] += this.mm[i][k] * this.mm[i][k];
    tot += this.mm[i][0];
  }
  // normalize the fluctuations
  var s = n / (tot * tot);
  for ( k = 0; k < 3; k++ ) {
    fl[k] *= s;
    if ( k == 0 ) fl[0] -= 1;
    flsum += fl[k];
  }
  for ( k = 0; k < 3; k++ )
    this.flfr[k] = fl[k]/flsum;
  return Math.sqrt(flsum);
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
  console.log("alpha " + this.alphawl + ", " + this.n/this.t + ", "
    + "fluc " + this.hfluc + ", invt " + this.invt);
  this.t = 0;
  for ( i = 0; i < n; i++ ) {
    // renew the accumulators
    this.mm[i] = [0, 0, 0];
  }
  for ( i = 0; i < n; i++ )
    for ( j = 0; j < xn; j++ )
      this.hist[i][j] = 0;
}



/* extensive check of histogram fluctuation */
AGE.prototype.wlcheckx = function(fl, magred)
{
  this.hfluc = this.calcfl_rms();
  if ( !this.invt && this.hfluc < fl ) {
    this.switch(magred);
    return 1;
  }
  return 0;
}




/* transition to a neighboring umbrella */
AGE.prototype.move = function(x, id, local)
{
  var acc = 0, n = this.n, jd;
  var xi, xj, vi, vj, dv, sigi, sigj;

  if ( local ) {
    // jump to a nearest neighbor
    jd = (rand01() < 0.5) ? id - 1 : id + 1;
    if ( jd < 0 || jd >= n ) {
      this.tacc = false;
      return id;
    }
  } else {
    // jump any other umbrella
    jd = ( id + 1 + Math.floor(rand01() * (n - 1)) ) % n;
  }
  sigi = this.sig[id];
  sigj = this.sig[jd];
  /* compute the acceptance probability */
  xi = (x - this.ave[ id]) / sigi;
  xj = (x - this.ave[ jd]) / sigj;
  vi = this.c0[id] + this.c1[id] * xi + this.c2[id] * (xi * xi - 1) / SQRT2;
  vj = this.c0[jd] + this.c1[jd] * xj + this.c2[jd] * (xj * xj - 1) / SQRT2;
  dv = vj - vi;
  //console.log(id, jd, x, vi, this.c0[id], vj, this.c0[jd], dv, local);
  this.tacc = ( dv <= 0 || rand01() < Math.exp(-dv) );
  if ( this.tacc ) {
    this.id = jd;
  }
  return this.id;
}

var LN0 = -1000000;

// ln(exp(x) + exp(y))
function lnadd(x, y)
{
  if ( x < y ) { // swap x and y
    var z = y; y = x; x = z;
  }
  y -= x;
  return x + Math.log(1 + Math.exp(y));
}

/* non-iterative histogram reweighting */
AGE.prototype.reweight = function()
{
  var ix, i, x, dxi, ui, lnden, htot, all = 0;
  var xn = this.xn, n = this.n, xmin = this.xmin, dx = this.dx;

  for ( ix = 0; ix < xn; ix++ ) {
    // compute the density of states
    // g[ix] = htot[ix] / Sum_i cnt[i] exp(-ui[ix]-c0[i])
    this.lng[ix] = lnden = LN0;
    x = xmin + ix * dx;
    htot = 0;
    for ( i = 0; i < n; i++ ) {
      htot += this.hist[i][ix];
    }
    this.htot[ix] = htot;
    all += htot;
    if ( htot <= 0 ) continue;
    for ( i = 0; i < n; i++ ) {
      if ( this.cnt[i] <= 0 ) continue;
      dxi = (x - this.ave[i])/this.sig[i];
      ui = this.c0[i] + this.c1[i] * dxi + this.c2[i] * (dxi*dxi - 1) / SQRT2;
      lnden = lnadd(lnden, Math.log(this.cnt[i]) - ui);
    }
    this.lng[ix] = Math.log(htot) - lnden;
    //console.log(ix, this.lng[ix]);
  }
  return all;
}

/* find the peak of arr[i] - x_i * bc in [ileft, iright) */
function findpeak(arr, bc, ileft, iright, xmin, dx)
{
  var x, y, i, imax = ileft, ym = -1e100;

  for ( i = ileft; i < iright; i++ ) {
    if ( arr[i] <= LN0 ) continue;
    x = xmin + i * dx;
    y = arr[i] - x * bc;
    if ( y > ym ) {
      ym = y;
      imax = i;
    }
  }
  return [imax, ym];
}

/* find the critical inverse temperature */
AGE.prototype.seekcrit = function()
{
  var bc, beta, y1, y2, y, lns, lns2;
  var i, ix, ix1, ix2, im, t, x, ret;
  var xn = this.xn, n = this.n, dx = this.dx, xmin = this.xmin;
  var arr = this.lng;

  // since c1/sigma is roughly beta, we will estimate the critical temperature
  // at a point where c1/sigma rises quickest
  var dbmax = 0;
  for ( i = 0; i < n - 1; i++ ) {
    var b0 = this.c1[i]/this.sig[i];
    var b1 = this.c1[i+1]/this.sig[i+1];
    var db = (b1 - b0)/(this.ave[i+1] - this.ave[i]);
    if ( db > dbmax ) {
      dbmax = db;
      bc = b0;
      im = Math.floor((this.ave[i] - this.xmin)/this.dx);
    }
  }

  // adjust the critical temperature util the heights are equal
  for ( t = 0; t < 10; t++ ) {
    ret = findpeak(arr, bc, 0,  im, xmin, dx);
    ix1 = ret[0]; y1 = ret[1];
    ret = findpeak(arr, bc, im, xn, xmin, dx);
    ix2 = ret[0]; y2 = ret[1];
    console.log("iter", t, "bc", bc, ix1, y1, ix2, y2, im);
    if ( Math.abs(y2-y1) < 1e-10 ) break;
    bc += (y2 - y1) / ((ix2 - ix1)*dx);
  }
  this.bc = bc;
  this.x1 = xmin + ix1*dx;
  this.x2 = xmin + ix2*dx;

  // normalize lng at the critical point
  lns = LN0;
  for ( ix = 0; ix < xn; ix++ ) {
    if ( arr[ix] <= LN0 ) continue;
    x = xmin + ix * dx;
    y = arr[ix] - x * bc;
    lns = lnadd(lns, y);
  }
  lns += Math.log(dx);
  // normalize the density of states
  for ( ix = 0; ix < xn; ix++ )
    arr[ix] -= lns;

  // normalize c0 at the critical temperature
  lns2 = LN0;
  for ( i = 0; i < n; i++ )
    lns2 = lnadd(lns2, this.c0[i] - this.c2[i]/SQRT2 - bc * this.ave[i]);
  lns2 += Math.log(this.delx);
  this.shift = lns;
  this.shiftc0hat = lns2;

  return bc;
}



