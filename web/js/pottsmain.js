/* Handle web interface */



"use strict";



var potts = null;

var age = null;

var iage = 0; // current umbrella index
//var tempfreq = 0.5; // tempering frequency
//var hist = null; // energy histogram
//var tpmctot = 1e-16, tpmcacc = 0;

var esig = 16;
var espacing = 16;
var betac = 1.4;

var flthreshold = 0.5;
var alpha0 = 0.01;
var magred = 0.5;

var lnzmethod = "WL";
var localmove_prob = 1.0;

var timer_interval = 100; // in milliseconds
var potts_timer = null;
var mc_algorithm = "Metropolis";

var blocksizemc = 10; // number of steps for a block in the Metropolis algorithm
var nstepspsmc = 1000; // number of steps per second for MC
var nstepspfmc = 100;  // number of steps per frame for MC

var histplot = null;
var vplot = null;

function getparams()
{
  var l = get_int("L", 16);
  var q = get_int("q", 10);
  potts = new Potts(l, q);

  var ecmin = -Math.floor(1.8 * potts.n + 0.5);
  var ecmax = -Math.floor(0.8 * potts.n + 0.5);
  esig = potts.l * get_float("sigma", 1.0);
  espacing = potts.l * get_float("spacing", 1.0);
  potts_equil(ecmin);

  alpha0 = get_float("alpha0", 0.01);
  flthreshold = get_float("flthreshold", 0.5);
  magred = get_float("magred", 0.5);
  betac = get_float("betac", 1.4);
  age = new AGE(ecmin, ecmax, espacing, esig, lnzmethod,
      betac * esig, alpha0, -2*potts.n, 0, 1, 0);
  iage = 0;

  mc_algorithm = grab("mc_algorithm").value;

  blocksizemc = get_int("blocksizemc", 10);
  nstepspsmc = get_int("nstepspersecmc", 10000);
  nstepspfmc = nstepspsmc * timer_interval / 1000;
}



function changescale()
{
  paint();
}


/* equilibrate the system to raise the energy above ene */
function potts_equil(ene)
{
  var id, sn, h, t;

  for ( t = 1; ; t++ ) {
    id = potts.pick();
    potts.flip(id, potts.sn, potts.h);
    if ( potts.E > ene ) break;
  }
}


/* modified Metropolis algorithm
 * p(E) ~ exp(-b1*(E-Eave)-b2*(E-Eave)^2/2) */
function potts_metro_mod(b1, b2, Eave)
{
  var id, sn, h, enew, edev, dv, acc;

  id = potts.pick();
  enew = potts.E + potts.h;
  edev  = (enew + potts.E) * 0.5 - Eave;
  dv = potts.h * (b1 + b2 * edev);
  if ( dv <= 0 || rand01() < Math.exp(-dv) ) {
    potts.flip(id, potts.sn, potts.h);
    return 1;
  } else {
    return 0;
  }
}


function potts_wolff_mod(b1, b2, Eave)
{
  var q = potts.q;
  var l = potts.l, n = potts.n, i, ix, iy, id, so, sn, cnt = 0, h = 0;
  var padd, enew, dv;

  padd = 1 - Math.exp(-b1);
  // randomly selected a seed */
  id = Math.floor( rand01() * n );
  so = potts.s[id];
  sn = (so + 1 + Math.floor(rand01() * (q - 1))) % q;
  potts.queue[ cnt++ ] = id;
  for ( i = 0; i < n; i++ )
    potts.used[i] = false;
  potts.used[id] = true;

  /* go through spins in the queue */
  for ( i = 0; i < cnt; i++ ) {
    id = potts.queue[i];
    potts.s[id] = sn;
    /* add neighbors of i with the same spins */
    ix = id % l;
    iy = id - ix;
    h += potts.addtoqueue(iy + (ix + 1) % l,     so, sn, padd);
    h += potts.addtoqueue(iy + (ix + l - 1) % l, so, sn, padd);
    h += potts.addtoqueue((iy + l) % n + ix,     so, sn, padd);
    h += potts.addtoqueue((iy + n - l) % n + ix, so, sn, padd);
  }

  enew = potts.E + h;
  dv = h * b2 * ((enew + potts.E)*.5 - Eave);
  if ( dv <= 0 || rand01() < Math.exp(-dv) ) {
    potts.E += h;
    return 1;
  } else { /* reject */
    for ( i = 0; i < cnt; i++ ) {
      id = potts.queue[i];
      potts.s[id] = so;
    }
    return 0;
  }
}


function paint()
{
  if ( !potts ) {
    return;
  }
  var c = grab("animationbox");
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;

  // the system dimension is L + 1
  var dx = 1.0 * width / (potts.l);
  var dy = 1.0 * height / (potts.l);
  var dx1 = Math.floor(dx - 1), dy1 = Math.floor(dy - 1);
  if ( dx1 < 1 ) dx1 = 1;
  if ( dy1 < 1 ) dy1 = 1;

  // draw the background
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  // draw each spin
  var l = potts.l, q = potts.q, id = 0;
  var s, colors = new Array(potts.q);
  for ( s = 0; s < potts.q; s++ ) {
    colors[s] = getHueColor(1.0*s/potts.q);
  }

  for ( var i = 0; i < l; i++ ) {
    for ( var j = 0; j < l; j++, id++ ) {
      var x = Math.floor( (i - 0.5*l) * dx + width * 0.5 );
      var y = Math.floor( (j - 0.5*l) * dy + height * 0.5 );
      s = potts.s[id];
      ctx.fillStyle = colors[s];
      ctx.fillRect(x, y, dx1, dy1)
    }
  }
}



function updatehistplot()
{
  var i, j, ntp = age.n, xn = age.xn;
  var dat = "Energy";

  // prepare the header
  for ( j = 0; j < ntp; j++ ) {
    dat += ",Histogram " + (j+1);
  }
  dat += "\n";

  // refine the energy range
  var imin = 0;
  for ( ; imin <= xn; imin++ ) {
    for ( j = 0; j < ntp; j++ )
      if ( age.hist[j][imin] > 0 ) break;
    if ( j < ntp ) break;
  }

  var imax = xn - 1;
  for ( ; imax >= imin; imax-- ) {
    for ( j = 0; j < ntp; j++ )
      if ( age.hist[j][imax] > 0 ) break;
    if ( j < ntp ) break;
  }
  if ( imin > imax ) return;

  // fill in the energy histogram data
  for ( i = imin; i <= imax; i++ ) {
    var ep = -2 * potts.n + i;
    dat += "" + ep;
    for ( j = 0; j < ntp; j++ ) {
      dat += "," + age.hist[j][i];
    }
    dat += "\n";
  }

  if ( histplot === null || 1 ) {
    var h = grab("animationbox").height / 2 - 5;
    var w = h * 3 / 2;
    var options = {
      xlabel: "<small>Energy</small>",
      ylabel: "<small>Histogram</small>",
      includeZero:true,
      axisLabelFontSize: 10,
      xRangePad: 2,
      plotter: barChartPlotter,
      width: w,
      height: h
    };
    histplot = new Dygraph(grab("histplot"), dat, options);
  } else {
    histplot.updateOptions({file: dat});
  }

  dat = null;
}



function updatevplot()
{
  var i, j, ntp = age.n;
  // prepare the header
  var dat = "Energy,c1/sig\n";

  // fill in the energy histogram data
  for ( j = 0; j < ntp; j++ ) {
    var ene = age.ave[j];
    dat += "" + ene + "," + (age.c1[j]/esig) + "\n";
  }

  if ( vplot === null || 1 ) {
    var h = grab("animationbox").height / 2 - 5;
    var w = h * 3 / 2;
    var options = {
      xlabel: "<small>Energy</small>",
      ylabel: "<small><i>c</i><sub>1</sub>/<i>&sigma;</i></small>",
      //includeZero:true,
      axisLabelFontSize: 10,
      xRangePad: 2,
      drawPoints: true,
      pointSize: 2,
      width: w,
      height: h
    };
    vplot = new Dygraph(grab("vplot"), dat, options);
  } else {
    vplot.updateOptions({file: dat});
  }
}



/* draw a color bar */
function drawcolorbar(target)
{
  var c = grab(target);
  var ctx = c.getContext("2d");
  var width = c.width;
  var height = c.height;

  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);

  var grd = ctx.createLinearGradient(0, 0, width, 0);
  grd.addColorStop(0, "#0000cc");
  grd.addColorStop(1, "#cc0000");

  ctx.fillStyle = grd;
  ctx.fillRect(0, 4, width, height-4);

  // draw grids
  if ( age ) {
    var i, x;
    for ( i = 0; i < age.n; i++ ) {
      x = width * i / (age.n - 1.0);
      drawLine(ctx, x, 4, x, height, "#808080", 1);
    }
    x = width * iage / (age.n - 1.0);
    drawLine(ctx, x, 0, x, height, "#000000", 3);
  }
}



function pulse()
{
  var sinfo, istep, jstep;

  for ( istep = 0; istep <= nstepspfmc; istep++ ) {
    var beta1 = age.c1[iage] / esig;
    //TODO:
    //age.c2[iage] = 1.0;
    var beta2 = age.c2[iage] / (esig * esig);
    var Eave = age.ave[iage];
    var acc = 0;
    if ( mc_algorithm === "Metropolis" ) {
      for ( jstep = 0; jstep < blocksizemc; jstep++ ) {
        acc += potts_metro_mod(beta1, beta2, Eave);
      }
      acc = 1.0*acc/blocksizemc;
    } else if ( mc_algorithm === "Wolff" ) {
      acc = potts_wolff_mod(beta1, beta2, Eave);
    }
    iage = age.move(potts.E, iage, rand01() <= localmove_prob);
    age.add(iage, potts.E, acc);
  }

  //var e = potts.E; potts.energy(); console.log("energy " + e + " vs. " + potts.E, esig);

  age.wlcheckx(flthreshold, magred);
  drawcolorbar("tpscale");
  sinfo = "Umbrella " + iage + "/" + age.n
    + ", ave " + age.ave[iage] + ", sig " + esig
    + ", c<sub>1</sub>/&sigma; " + roundto(age.c1[iage]/esig, 2)
    + ", c<sub>2</sub> " + roundto(age.c2[iage], 2)
    + ", E/N " + roundto(1.0*potts.E/potts.n, 2) + ", data " + age.cnt[iage]
    + ", &alpha; " + age.getalpha().toPrecision(3) + ", 1/t " + age.invt
    + ", fl " + roundto(age.hfluc, 2)
    + " (" + roundto(age.flfr[0]*100.0,2) + "%/" + roundto(age.flfr[1]*100.0,2) + "%/" + roundto(age.flfr[2]*100.0,2) + "%)";
  grab("sinfo").innerHTML = sinfo;
  //if ( age.cnt[iage] > 10000 ) age.invt = true;

  paint();
  updatehistplot();
  updatevplot();
}



function stopsimul()
{
  if ( potts_timer !== null ) {
    clearInterval(potts_timer);
    potts_timer = null;
  }
  grab("pause").innerHTML = "&#9724;";
  //tpmctot = 1e-16;
  //tpmcacc = 0;
}



function pausesimul()
{
  if ( !potts ) {
    return;
  }
  if ( potts_timer !== null ) {
    clearInterval(potts_timer);
    potts_timer = null;
    grab("pause").innerHTML = "&#10704";
  } else {
    potts_timer = setInterval(
        function() { pulse(); },
        timer_interval);
    grab("pause").innerHTML = "&#9724;";
  }
}



function startsimul()
{
  stopsimul();
  getparams();
  //installmouse("animationbox", "animationboxscale");
  potts_timer = setInterval(
    function(){ pulse(); },
    timer_interval );
}



function pausesimul2()
{
  // skip a mouse-move
  //if ( mousemoved > 0 ) {
  //  return;
  //}
  if ( !potts ) {
    startsimul();
  //} else if ( mousemoved === 0 ) {
  } else {
    pausesimul();
  }
}


function runwham()
{
  if ( potts_timer ) {
    pausesimul();
  }
}


/* respond to critical parameter changes: restart simulation */
function changeparams()
{
  if ( potts_timer !== null ) {
    startsimul();
  }
}



function showtab(who)
{
  who = grab(who);
  var par = who.parentNode;
  var c = par.childNodes;
  var i, iwho, k = 0;

  // arrange the tabs
  for ( i = 0; i < c.length; i++ ) {
    if ( c[i].className === "params-panel" ) {
      if ( c[i] !== who ) {
        c[i].style.zIndex = k;
        k += 1;
      } else {
        iwho = k;
      }
    }
  }
  who.style.zIndex = k;

  // arrange the clickable tab titles
  k += 1;
  var pt = grab("tabsrow");
  pt.style.zIndex = k;
  var ct = pt.childNodes, ik = 0;
  for ( i = 0; i < ct.length; i++ ) {
    if ( ct[i].tagName ) {
      if ( ik === iwho ) {
        ct[i].style.fontWeight = "bold";
        ct[i].style.borderTop = "2px solid #c0c0d0";
      } else {
        ct[i].style.fontWeight = "normal";
        ct[i].style.borderTop = "0px solid #e0e0f0";
      }
      ik += 1;
    }
  }
}



function resizecontainer(a)
{
  var canvas = grab("animationbox");
  var ctx = canvas.getContext("2d");
  var w, h;
  if ( a === null || a === undefined ) {
    w = canvas.width;
    h = canvas.height;
  } else {
    a = parseInt( grab(a).value );
    w = h = a;
    canvas.width = w;
    canvas.height = h;
  }
  ctx.font = "24px Verdana";
  ctx.fillText("Click to start", w/2-40, h/2-10);

  var hsbar = 30; // height of the global scaling bar
  var hcbar = 40; // height of the control bar
  var htbar = 30; // height of the tabs bar
  var wr = h*3/4; // width of the plots
  var wtab = 560; // width of the tabs
  var htab = 280;

  grab("simulbox").style.width = "" + w + "px";
  grab("simulbox").style.height = "" + h + "px";
  grab("simulbox").style.top = "" + hsbar + "px";
  grab("controlbox").style.top = "" + (h + hsbar) + "px";
  //grab("animationboxscale").style.width = "" + (w - 100) + "px";
  grab("tpscale").style.width = "" + (w - 220) + "px";
  drawcolorbar("tpscale");
  histplot = null;
  grab("histplot").style.left = "" + w + "px";
  grab("histplot").style.width = "" + wr + "px";
  grab("vplot").style.top = "" + hcbar + "px";
  grab("histplot").style.height = "" + h/2 + "px";
  vplot = null;
  grab("vplot").style.left = "" + w + "px";
  grab("vplot").style.width = "" + wr + "px";
  grab("vplot").style.top = "" + (h/2 + hcbar) + "px";
  grab("vplot").style.height = "" + h/2 + "px";
  grab("tabsrow").style.top = "" + (h + hsbar + hcbar) + "px";
  grab("tabsrow").style.width = "" + wtab + "px";

  var c = grab("container").childNodes;
  var i;
  // tabs
  for ( i = 0; i < c.length; i++ ) {
    if ( c[i].className === "params-panel" ) {
      c[i].style.top = "" + (h + hsbar + hcbar + htbar) + "px";
      c[i].style.width = "" + (wtab - 20) + "px";
      c[i].style.height = "" + htab + "px";
    }
  }
  grab("sinfo").style.top = "" + (h + hsbar + hcbar + htbar) + "px";
  grab("sinfo").style.left = "" + (wtab + 10) + "px";
  grab("sinfo").style.width = "" + (w + wr - wtab - 20) + "px";
  grab("sinfo").innerHTML = "simulation information";

  grab("container").style.height = "" + (h + hsbar + hcbar + htbar + htab) + "px";
  grab("container").style.width = "" + (w + wr) + "px";
}


function init()
{
  resizecontainer();
  showtab("system-params");
}
