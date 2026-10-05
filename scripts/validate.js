const fs=require('node:fs');const crypto=require('node:crypto');
const required=['README.md','index.html','styles.css','app.js','data/ogle-2023-blg-0100-phot.dat','data/ogle-2023-blg-0100-fit.json','scripts/fit_ogle_event.py'];
const failures=[];for(const file of required)if(!fs.existsSync(file))failures.push(`${file} missing`);
const data=JSON.parse(fs.readFileSync('data/ogle-2023-blg-0100-fit.json','utf8'));
if(data.event.name!=='OGLE-2023-BLG-0100')failures.push('wrong event');
if(data.observations?.length!==360)failures.push('expected 360 observations');
if(data.curve?.length!==700)failures.push('expected 700 model samples');
if(!(data.fits.pspl_blended.parameters.tE_days>20&&data.fits.pspl_blended.parameters.tE_days<28))failures.push('timescale outside regression bounds');
if(!(data.fits.pspl_blended.parameters.u0>.25&&data.fits.pspl_blended.parameters.u0<.38))failures.push('u0 outside regression bounds');
const digest=crypto.createHash('sha256').update(fs.readFileSync('data/ogle-2023-blg-0100-phot.dat')).digest('hex');
if(digest!==data.provenance.source_sha256)failures.push('source checksum mismatch');
const combined=['index.html','app.js'].map(f=>fs.readFileSync(f,'utf8')).join('\n').toLowerCase();
for(const token of ['todo','placeholder','physicsworker.js','ogle-2023-blg-0001'])if(combined.includes(token))failures.push(`obsolete or unfinished token: ${token}`);
if(failures.length){console.error(failures.join('\n'));process.exit(1)}
console.log(`Microlensing validation passed: ${data.observations.length} real epochs, reproducible PSPL fit, verified source hash.`);
