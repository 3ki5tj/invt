/* Generalized ensemble of adaptive umbrella sampling
 * */


"use strict";

var SQRT2 = 1.4142135623730951;

function MMWL()
{
  this.mm = [0, 0, 0];
  this.fl = [0, 0, 0];
}

MMWL.prototype.add = function(y1, y2)
{
  this.mm[0] += 1;
  this.mm[1] += y1;
  this.mm[2] += y2;
}

/* compute the fluctuation of the statistical moments
 * return the maximum fluctuation */
MMWL.prototype.calcfl = function()
{
  var t = this.mm[0];
  if ( t <= 0 ) return 99.0;
  this.fl[1] = this.mm[1] / t;
  this.fl[2] = this.mm[2] / t;
  var fl1 = Math.abs(this.fl[1]);
  var fl2 = Math.abs(this.fl[2]);
  return Math.max(fl1, fl2);
}


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
  this.mmwl = new Array(n);
  this.cnt = new Array(n);
  this.acc = new Array(n);
  this.hfl = new Array(n);
  for ( i = 0; i < n; i++ ) {
    this.ave[i] = xcmin + i * delx;
    this.sig[i] = sig;
    this.c0[i] = (i - (n - 1)*0.5) * c1 / sig * delx;
    this.c1[i] = c1;
    this.c2[i] = 0;
    this.mmwl[i] = new MMWL();
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
}



/* retrieve the updating magnitude */
AGE.prototype.getalpha = function(i)
{
  var alpha = this.invt ? this.n / (this.t + this.t0) : this.alphawl;
  return alpha;
}



/* update the parameters of the bias potential */
AGE.prototype.add = function(i, x, acc)
{
  var ave = this.ave[i], sig = this.sig[i];
  var alpha = this.getalpha(i);

  var y1 = (x - ave) / sig;
  var y2 = (y1 * y1 - 1) / SQRT2;
  this.c1[i] += y1 * alpha;
  this.c2[i] += y2 * alpha;
  this.mmwl[i].add(y1, y2);

  this.t += 1;
  this.cnt[i] += 1;
  this.acc[i] += acc;
  var j = Math.floor( (x - this.xmin) / this.dx );
  this.hist[i][j] += 1;
  this.c0[i] += alpha;
}



/* calculate the fluctuation of histogram modes (from the variance) */
AGE.prototype.calcfl = function()
{
  var i, n = this.n, mm, nhh;
  var tot = 0, fl0 = 0, fl1 = 0, fl2 = 0, fl, hi;

  /* compute the total sample size */
  for ( i = 0; i < n; i++ )
    tot += this.mmwl[i].mm[0];
  if ( tot <= 0 ) return 99.0;

  for ( i = 0; i < n; i++ ) {
    mm = this.mmwl[i];
    mm.calcfl();
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
  console.log("alpha " + this.alphawl + ", " + this.n/this.t + ", "
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
  this.hfluc = this.calcfl();
  if ( !this.invt && this.hfluc < fl ) {
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
  dv = this.c0[jd] - this.c0[id];
  /* compute the acceptance probability */
  xi = (x - this.ave[ id]) / sigi;
  xj = (x - this.ave[ jd]) / sigj;
  vi = this.c1[ id] * xi + this.c2[ id] * (xi * xi - 1) / SQRT2;
  vj = this.c1[ jd] * xj + this.c2[ jd] * (xj * xj - 1) / SQRT2;
  dv += vj - vi;
  //console.log(id, jd, x, vi, this.c0[id], vj, this.c0[jd], dv, local);
  this.tacc = ( dv <= 0 || rand01() < Math.exp(-dv) );
  if ( this.tacc ) {
    this.id = jd;
  }
  return this.id;
}

