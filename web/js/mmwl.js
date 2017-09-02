/* Wang-Landau scheme for moments */



"use strict";



function MMWL(a)
{
  this.mm = [0, 0, 0];
  this.fl = [0, 0, 0];
  this.alphawl = a;
  this.invt = false;
  this.t0 = 1;
}



MMWL.prototype.add = function(y1, y2)
{
  this.mm[0] += 1;
  this.mm[1] += y1;
  this.mm[2] += y2;
}



MMWL.prototype.getalpha = function()
{
  return this.invt ? 1./(this.mm[0] + this.t0) : this.alphawl;
}



/* compute the fluctuation of the statistical moments
 * return the maximum fluctuation */
MMWL.prototype.calcfl = function(dovar)
{
  var t = this.mm[0], fl, fl2;
  if ( t <= 0 ) return 99.0;
  this.fl[1] = this.mm[1] / t;
  this.fl[2] = this.mm[2] / t;
  fl = Math.abs(this.fl[1]);
  if ( dovar && (fl2 = Math.abs(this.fl[2])) > fl )
    fl = fl2;
  return fl;
}



/* switch to a stage with a reduced updating magnitude */
MMWL.prototype.switch = function(magred)
{
  var t = this.mm[0];
  this.alphawl *= magred;
  if ( this.alphawl < 1./t ) {
    this.t0 = t;
    this.invt = true;
  }
  this.mm = [0, 0, 0];
}


MMWL.prototype.check = function(dovar, fl, magred)
{
  var flmax = this.calcfl(dovar);
  if ( !this.invt && flmax < fl ) {
    this.switch(magred);
    return 1;
  } else {
    return 0;
  }
}

