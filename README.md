# nztm2000
nztm2000 coordinates to wgs84. This version is translated from [the official c language](https://www.linz.govt.nz/data/geodetic-services/download-geodetic-software) and [leighghunt's c# version](https://github.com/leighghunt/nztm).

# Usage

```
// NZTM2000: E1817224, N5675344
// Latitude: -39.04398599, Longitude: 175.50998658

var ltlg = nztm2000.nztm_geod('1817224', '5675344');

console.log(ltlg); // {lt: -39.04398599563426, ln: 175.50998657532008}

var lglt = nztm2000.geod_nztm('-39.04398599','175.50998658');

console.log(lglt); // {ce: 1817224.0004226866, cn: 5675344.000490918}
```
