pro dr1_astrom_performance

; Astrometric accuracy and precision of final object coordinates
; by comparing to GAIA values

files = file_search('/data/smash/cp/red/photred/catalogs/final/v4/Field*_combined_allobj_bright.fits',count=nfiles)

fitstr = replicate({field:'',nmatch:0L,ramed:0.0,rasig:0.0,decmed:0.0,decsig:0.0},nfiles)

for i=0,nfiles-1 do begin

  field = file_basename(files[i],'_combined_allobj_bright.fits')

  str = mrdfits(files[i],1,/silent)
  gaiafile = '/data/smash/cp/red/photred/gaia/'+field+'_gaia.fits.gz'
  gaia = mrdfits(gaiafile,1,/silent)
  srcmatch,gaia.ra_icrs,gaia.de_icrs,str.ra,str.dec,0.5,ind1,ind2,/sph,count=nmatch
  radiff = (gaia[ind1].ra_icrs-str[ind2].ra)*cos(gaia[ind1].de_icrs/!radeg)*3600*1e3
  decdiff = (gaia[ind1].de_icrs-str[ind2].dec)*3600.*1e3
  ;gd = where(str[ind2].g lt 19,ngd)
  ;ramed = median(radiff[gd])
  ;rasig = mad(radiff[gd])
  ;decmed = median(decdiff[gd])
  ;decsig = mad(decdiff[gd])
  ramed = median(radiff)
  rasig = mad(radiff)
  decmed = median(decdiff)
  decsig = mad(decdiff)

  fitstr[i].field = field
  fitstr[i].nmatch = nmatch
  fitstr[i].ramed = ramed
  fitstr[i].rasig = rasig
  fitstr[i].decmed = decmed
  fitstr[i].decsig = decsig

  form = '(I-5,A-10,I-8,4F-10.4)'
  print,format=form,i+1,field,nmatch,ramed,rasig,decmed,decsig

  ;stop

endfor

; Final values
print,'RA med:  ',median(fitstr.ramed),' +/- ',mad(fitstr.ramed),' mas'
print,'RA sig:  ',median(fitstr.rasig),' +/- ',mad(fitstr.rasig),' mas'
print,'DEC med: ',median(fitstr.decmed),' +/- ',mad(fitstr.decmed),' mas'
print,'DEC sig: ',median(fitstr.decsig),' +/- ',mad(fitstr.decsig),' mas'

;mwrfits,fitstr,'dr1_astrom_performance.fits',/create

stop

end
