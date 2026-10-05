import fs from 'node:fs';
const d=JSON.parse(fs.readFileSync('data/ogle-2023-blg-0100-fit.json','utf8'));const fail=[];
if(d.provenance.photometry_url!=='https://www.astrouw.edu.pl/ogle/ogle4/ews/2023/blg-0100/phot.dat')fail.push('canonical source URL changed');
if(d.comparison.delta_bic_constant_minus_pspl<10000)fail.push('constant-model rejection regression failed');
if(!(d.fits.pspl_blended.reduced_chi2>1&&d.fits.pspl_blended.reduced_chi2<2.5))fail.push('fit-quality regression failed');
if(d.comparison.delta_bic_unblended_minus_blended>=0)fail.push('BIC should prefer the simpler unblended model');
if(!d.scope.not_inferred.includes('unique lens mass'))fail.push('scientific boundary missing');
if(fail.length){console.error(fail.join('\n'));process.exit(1)}
console.log('Research validation passed: provenance, fit regression, model comparison and inference limits verified.');
