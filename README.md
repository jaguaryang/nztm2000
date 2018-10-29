# nztm2000
nztm2000 coordinates to wgs84. This version is translated from [the official c language](https://www.linz.govt.nz/data/geodetic-services/download-geodetic-software) and [leighghunt's c# version](https://github.com/leighghunt/nztm).

# Usage

// Example:
// NZTM2000: E1817224, N5675344
// Latitude: -39.04398599, Longitude: 175.50998658
var ltlg = nztm2000.nztm_geod('1817224', '5675344');
ltlg.lt *= rad2deg;
ltlg.ln *= rad2deg;
console.log(ltlg);
// output:
// {lt: -39.04398599563426, ln: 175.50998657532008}

var lglt = nztm2000.geod_nztm(
'-39.04398599' / rad2deg,
'175.50998658' / rad2deg
);
console.log(lglt);
// output:
// {ce: 1817224.0004226866, cn: 5675344.000490918}

## npm

npm install nztm

import nztm from "nztm";

nztm.test();

## Js Non-modular usage

// module.exports = nztm2000;
nztm2000.test();
