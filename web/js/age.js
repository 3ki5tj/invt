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
  this.delx = delx;
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
AGE.prototype.calcfl_rms = function()
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

var LN0 = -1000000;

// ln(exp(x) + exp(y))
function lnadd(x, y)
{
  if ( x < y ) { // swap x and y
    var z = y;
    y = x;
    x = z;
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
  var bc, beta, w, sw, y1, y2, y, lns, lns2;
  var ix, ix1, ix2, il = -1, ir, im, t, x, ret;
  var xn = this.xn, n = this.n, dx = this.dx, xmin = this.xmin;
  var arr = this.lng;

  bc = sw = 0;
  for ( ix = 0; ix < xn - 1; ix++ ) {
    if ( arr[ix] <= LN0 || arr[ix+1] < LN0 ) {
      continue;
    } else if ( il < 0 ) {
      il = ix;
    } else {
      ir = ix;
    }
    if ( this.htot[ix] > 0 && this.htot[ix+1] > 0 ) {
      beta = (arr[ix + 1] - arr[ix]) / dx;
      w = (this.htot[ix] + this.htot[ix+1]) / 2;
      bc += beta * w;
      sw += w;
    }
  }
  bc /= sw;
  im = Math.floor((il + ir)/2);
  console.log(bc, il, ir, im);
  //return bc;

  // find the two density peaks
  // util the heights are equal
  ix1 = il;
  ix2 = ir;
  for ( t = 0; t < 10; t++ ) {
    ret = findpeak(arr, bc, 0,  im, xmin, dx);
    ix1 = ret[0]; y1 = ret[1];
    ret = findpeak(arr, bc, im, xn, xmin, dx);
    ix2 = ret[0]; y2 = ret[1];
    //printf("%d: bc %g, x %d(%g) %d(%g) | %d(%d) %g\n", t, bc, xmin + ix1*dx, y1, xmin + ix2*dx, y2, xmin + im*dx, im, fabs(y2-y1));
    console.log(t, bc, ix1, y1, ix2, y2, im);
    if ( Math.abs(y2-y1) < 1e-14 ) break;
    bc += (y2 - y1) / ((ix2 - ix1)*dx);
  }
  this.bc = bc;
  this.x1 = xmin + ix1*dx;
  this.x2 = xmin + ix2*dx;
  //printf("bc %.8f, T %.8f, x %d %d\n", bc, 1/bc, xmin + ix1*dx, xmin + ix2*dx);

  /* normalize lng at the critical point */
  lns = LN0;
  for ( ix = 0; ix < xn; ix++ ) {
    if ( arr[ix] <= LN0 ) continue;
    x = xmin + ix * dx;
    y = arr[ix] - x * bc;
    lns = lnadd(lns, y);
  }
  lns += Math.log(dx);
  // normalize the density of states
  for ( ix = 0; ix < xn; ix++ ) {
    arr[ix] -= lns;
  }

  // normalize c0 at the critical temperature
  {
    lns2 = LN0;
    for ( var i = 0; i < n; i++ ) {
      y = this.c0[i] - this.c2[i]/SQRT2 - bc * this.ave[i];
      lns2 = lnadd(lns2, y);
      console.log(i, y, lns2);
      //printf("i %d, c0 %g, bc %g, ave %g, %g\n", i, c0[i], bc, ave[i], y);
    }
    lns2 += Math.log(this.delx);
  }
  //printf("lns for normalization %g (lng), %g (c0hat)\n", lns, lns + lns2);
  this.shift = lns;
  this.shiftc0hat = lns2;

  return bc;
}



