var nztm2000 = (function () {
  /* Define the parameters for the International Ellipsoid
     used for the NZGD2000 datum (and hence for NZTM) */

  var NZTM_A = 6378137;
  var NZTM_RF = 298.257222101;
  var NZTM_CM = 173.0;
  var NZTM_OLAT = 0.0;
  var NZTM_SF = 0.9996;
  var NZTM_FE = 1600000.0;
  var NZTM_FN = 10000000.0;

  /* Defines PI (from Abramowitz and Stegun Table 1.1) */

  //var  PI = 3.1415926535898
  var TWOPI = 2 * Math.PI;
  var rad2deg = (180 / Math.PI);

  /* Structure used to define a TM projection */

  /* Initiallize the TM structure  */

  function define_tmprojection(a, rf, cm, sf, lto, fe, fn, utom) {

    var f;

    tm = {};

    tm.meridian = cm;
    tm.scalef = sf;
    tm.orglat = lto;
    tm.falsee = fe;
    tm.falsen = fn;
    tm.utom = utom;
    if (rf != 0.0) f = 1.0 / rf; else f = 0.0;
    tm.a = a;
    tm.rf = rf;
    tm.f = f;
    tm.e2 = 2.0 * f - f * f;
    tm.ep2 = tm.e2 / (1.0 - tm.e2);

    tm.om = meridian_arc(tm, tm.orglat);
    return tm;
  }


  /***************************************************************************/
  /*                                                                         */
  /*  meridian_arc                                                           */
  /*                                                                         */
  /*  Returns the length of meridional arc (Helmert formula)                 */
  /*  Method based on Redfearn's formulation as expressed in GDA technical   */
  /*  manual at http://www.anzlic.org.au/icsm/gdatm/index.html               */
  /*                                                                         */
  /*  Parameters are                                                         */
  /*    projection                                                           */
  /*    latitude (radians)                                                   */
  /*                                                                         */
  /*  Return value is the arc length in metres                               */
  /*                                                                         */

  /***************************************************************************/


  function meridian_arc(tm, lt) {
    var e2 = tm.e2;
    var a = tm.a;
    var e4;
    var e6;
    var A0;
    var A2;
    var A4;
    var A6;

    e4 = e2 * e2;
    e6 = e4 * e2;

    A0 = 1 - (e2 / 4.0) - (3.0 * e4 / 64.0) - (5.0 * e6 / 256.0);
    A2 = (3.0 / 8.0) * (e2 + e4 / 4.0 + 15.0 * e6 / 128.0);
    A4 = (15.0 / 256.0) * (e4 + 3.0 * e6 / 4.0);
    A6 = 35.0 * e6 / 3072.0;

    return a * (A0 * lt - A2 * Math.sin(2 * lt) + A4 * Math.sin(4 * lt) - A6 * Math.sin(6 * lt));
  }

  /*************************************************************************/
  /*                                                                       */
  /*   foot_point_lat                                                      */
  /*                                                                       */
  /*   Calculates the foot point latitude from the meridional arc          */
  /*   Method based on Redfearn's formulation as expressed in GDA technical*/
  /*   manual at http://www.anzlic.org.au/icsm/gdatm/index.html            */
  /*                                                                       */
  /*   Takes parameters                                                    */
  /*      tm definition (for scale factor)                                 */
  /*      meridional arc (metres)                                          */
  /*                                                                       */
  /*   Returns the foot point latitude (radians)                           */
  /*                                                                       */

  /*************************************************************************/


  function foot_point_lat(tm, m) {
    var f = tm.f;
    var a = tm.a;
    var n;
    var n2;
    var n3;
    var n4;
    var g;
    var sig;
    var phio;

    n = f / (2.0 - f);
    n2 = n * n;
    n3 = n2 * n;
    n4 = n2 * n2;

    g = a * (1.0 - n) * (1.0 - n2) * (1 + 9.0 * n2 / 4.0 + 225.0 * n4 / 64.0);
    sig = m / g;

    phio = sig + (3.0 * n / 2.0 - 27.0 * n3 / 32.0) * Math.sin(2.0 * sig)
      + (21.0 * n2 / 16.0 - 55.0 * n4 / 32.0) * Math.sin(4.0 * sig)
      + (151.0 * n3 / 96.0) * Math.sin(6.0 * sig)
      + (1097.0 * n4 / 512.0) * Math.sin(8.0 * sig);

    return phio;
  }


  /***************************************************************************/
  /*                                                                         */
  /*   tmgeod                                                                */
  /*                                                                         */
  /*   Routine to convert from Tranverse Mercator to latitude and longitude. */
  /*   Method based on Redfearn's formulation as expressed in GDA technical  */
  /*   manual at http://www.anzlic.org.au/icsm/gdatm/index.html              */
  /*                                                                         */
  /*   Takes parameters                                                      */
  /*      input easting (metres)                                             */
  /*      input northing (metres)                                            */
  /*      output latitude (radians)                                          */
  /*      output longitude (radians)                                         */
  /*                                                                         */

  /***************************************************************************/

  function tm_geod(tm, ce, cn) {
    var fn = tm.falsen;
    var fe = tm.falsee;
    var sf = tm.scalef;
    var e2 = tm.e2;
    var a = tm.a;
    var cm = tm.meridian;
    var om = tm.om;
    var utom = tm.utom;
    var cn1;
    var fphi;
    var slt;
    var clt;
    var eslt;
    var eta;
    var rho;
    var psi;
    var E;
    var x;
    var x2;
    var t;
    var t2;
    var t4;
    var trm1;
    var trm2;
    var trm3;
    var trm4;

    cn1 = (cn - fn) * utom / sf + om;
    fphi = foot_point_lat(tm, cn1);
    slt = Math.sin(fphi);
    clt = Math.cos(fphi);

    eslt = (1.0 - e2 * slt * slt);
    eta = a / Math.sqrt(eslt);
    rho = eta * (1.0 - e2) / eslt;
    psi = eta / rho;

    E = (ce - fe) * utom;
    x = E / (eta * sf);
    x2 = x * x;


    t = slt / clt;
    t2 = t * t;
    t4 = t2 * t2;

    trm1 = 1.0 / 2.0;

    trm2 = ((-4.0 * psi
      + 9.0 * (1 - t2)) * psi
      + 12.0 * t2) / 24.0;

    trm3 = ((((8.0 * (11.0 - 24.0 * t2) * psi
      - 12.0 * (21.0 - 71.0 * t2)) * psi
      + 15.0 * ((15.0 * t2 - 98.0) * t2 + 15)) * psi
      + 180.0 * ((-3.0 * t2 + 5.0) * t2)) * psi + 360.0 * t4) / 720.0;

    trm4 = (((1575.0 * t2 + 4095.0) * t2 + 3633.0) * t2 + 1385.0) / 40320.0;

    var lt = fphi + (t * x * E / (sf * rho)) * (((trm4 * x2 - trm3) * x2 + trm2) * x2 - trm1);

    trm1 = 1.0;

    trm2 = (psi + 2.0 * t2) / 6.0;

    trm3 = (((-4.0 * (1.0 - 6.0 * t2) * psi
      + (9.0 - 68.0 * t2)) * psi
      + 72.0 * t2) * psi
      + 24.0 * t4) / 120.0;

    trm4 = (((720.0 * t2 + 1320.0) * t2 + 662.0) * t2 + 61.0) / 5040.0;

    var ln = cm - (x / clt) * (((trm4 * x2 - trm3) * x2 + trm2) * x2 - trm1);
    return {'lt': lt, 'ln': ln};
  }


  /***************************************************************************/
  /*                                                                         */
  /*   geodtm                                                                */
  /*                                                                         */
  /*   Routine to convert from latitude and longitude to Transverse Mercator.*/
  /*   Method based on Redfearn's formulation as expressed in GDA technical  */
  /*   manual at http://www.anzlic.org.au/icsm/gdatm/index.html              */
  /*   Loosely based on FORTRAN source code by J.Hannah and A.Broadhurst.    */
  /*                                                                         */
  /*   Takes parameters                                                      */
  /*      input latitude (radians)                                           */
  /*      input longitude (radians)                                          */
  /*      output easting  (metres)                                           */
  /*      output northing (metres)                                           */
  /*                                                                         */

  /***************************************************************************/


  function geod_tm(tm, ln, lt) {
    var fn = tm.falsen;
    var fe = tm.falsee;
    var sf = tm.scalef;
    var e2 = tm.e2;
    var a = tm.a;
    var cm = tm.meridian;
    var om = tm.om;
    var utom = tm.utom;
    var dlon;
    var m;
    var slt;
    var eslt;
    var eta;
    var rho;
    var psi;
    var clt;
    var w;
    var wc;
    var wc2;
    var t;
    var t2;
    var t4;
    var t6;
    var trm1;
    var trm2;
    var trm3;
    var gce;
    var trm4;
    var gcn;

    dlon = ln - cm;
    while (dlon > Math.PI) dlon -= TWOPI;
    while (dlon < -Math.PI) dlon += TWOPI;

    m = meridian_arc(tm, lt);

    slt = Math.sin(lt);

    eslt = (1.0 - e2 * slt * slt);
    eta = a / Math.sqrt(eslt);
    rho = eta * (1.0 - e2) / eslt;
    psi = eta / rho;

    clt = Math.cos(lt);
    w = dlon;

    wc = clt * w;
    wc2 = wc * wc;

    t = slt / clt;
    t2 = t * t;
    t4 = t2 * t2;
    t6 = t2 * t4;

    trm1 = (psi - t2) / 6.0;

    trm2 = (((4.0 * (1.0 - 6.0 * t2) * psi
      + (1.0 + 8.0 * t2)) * psi
      - 2.0 * t2) * psi + t4) / 120.0;

    trm3 = (61 - 479.0 * t2 + 179.0 * t4 - t6) / 5040.0;

    gce = (sf * eta * dlon * clt) * (((trm3 * wc2 + trm2) * wc2 + trm1) * wc2 + 1.0);
    var ce = gce / utom + fe;

    trm1 = 1.0 / 2.0;

    trm2 = ((4.0 * psi + 1) * psi - t2) / 24.0;

    trm3 = ((((8.0 * (11.0 - 24.0 * t2) * psi
      - 28.0 * (1.0 - 6.0 * t2)) * psi
      + (1.0 - 32.0 * t2)) * psi
      - 2.0 * t2) * psi
      + t4) / 720.0;

    trm4 = (1385.0 - 3111.0 * t2 + 543.0 * t4 - t6) / 40320.0;

    gcn = (eta * t) * ((((trm4 * wc2 + trm3) * wc2 + trm2) * wc2 + trm1) * wc2);
    var cn = (gcn + m - om) * sf / utom + fn;

    return {'ce': ce, 'cn': cn};
  }


  /* Define a static implementation of tmprojection */
  /* Note: for some implementations it may be better to create this
     dynamically and develop modified versions of the transformation
     functions to take this as a parameter */

  var nztm_projection;

  function get_nztm_projection() {
    if (nztm_projection == null) {
      nztm_projection = define_tmprojection(NZTM_A, NZTM_RF,
        NZTM_CM / rad2deg, NZTM_SF, NZTM_OLAT / rad2deg, NZTM_FE, NZTM_FN,
        1.0);
    }
    return nztm_projection;
  }

  /* Functions implementation the TM projection specifically for the
     NZTM coordinate system
  */

 return {
  nztm_geod: function (e, n) {
    var nztm = get_nztm_projection();
    var wgs84 = tm_geod(nztm, e, n);
    wgs84.lt *= rad2deg;
    wgs84.ln *= rad2deg;
    return wgs84;
  },
  geod_nztm: function (lt, ln) {
    var nztm = get_nztm_projection();
    return geod_tm(nztm, ln / rad2deg, lt / rad2deg);
  },
  test: function () {
    // Example:
    // NZTM2000: E1817224, N5675344
    // Latitude: -39.04398599, Longitude: 175.50998658
    var ltlg = nztm2000.nztm_geod(1817224, 5675344);
    console.log(ltlg);
    // output:
    // {lt: -39.04398599563426, ln: 175.50998657532008}

    var lglt = nztm2000.geod_nztm(-39.04398599, 175.50998658);
    console.log(lglt);
    // output:
    // {ce: 1817224.0004226866, cn: 5675344.000490918}
  },
};


})();

module.exports = nztm2000;
