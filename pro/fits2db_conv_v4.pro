pro fits2db_conv_v4,ifield,version

; Convert the SMASH binary fits files to the fits
; binary format needed to load the database

if n_elements(version) eq 0 then begin
  print,'Please enter the version number, e.g. "v2"'
  return
endif
dir = '/data/smash/cp/red/photred/catalogs/final/'+version+'/'
outdir = dir+'db/'
if file_test(outdir,/directory) eq 0 then file_mkdir,outdir

fieldstr = importascii('/data/smash/cp/red/photred/catalogs/pro/smash_fields_final.txt',/header)

fieldid = fix(strmid(strtrim(ifield,2),5))
print,ifield


; --- chips ---
; field, file, expnum, chip, base, filter, exptime, utdate,
; uttime, airmass, gain, rdnoise, nx, ny, wcstype, pixscale,
; ra, dec, wcsrms, fwhm, skymode, skysig, dao_nsources, dao_depth
; dao_npsfstars, dao_psftype, dao_psfboxsize, dao_psfvarorder,
; dao_psfchi, alf_nsources, alf_depth, calib_depth, calib_color,
; calib_zpterm, calib_amterm, calib_colorterm, calib_magname,
; apcor, ebv, refexpnum, vertices_ra, vertices_dec, nsrck
; allsrcindx, mjd, calibrated, ubercal_magoffset, ubercal_flag,
; photometric, badsoln, band, colband, colsign, zpterm, zptermsig,
; amterm, amtermsig, colterm, coltermsig, amcolterm, amcoltermsig,
; colsqterm, colsqtermsig
;  -strip out CALIB_XXX columns since we didn't calibrate using PHOTRED
;  -strip out ALLSRCINDX, AMCOL*, COLSQ* terms
;  -add JD, DATEOBS, break out vertices_RA/DEC?
;  -maybe add CALIB_ or CAL_ prefix for the calibration columns
;  -change NANs and 99.99 to something else?? 

; NEED TO ADD:
; -alftiletype, zpcalib_magoffset, zpcalib_magofferr, zpcalibflag
; leaving out alftiletype
; REMOVE:
; -gaiacal_magoffset, gaiacal_magofferr

chips = mrdfits(dir+ifield+'_combined_chips.fits.gz',1)
nchips = n_elements(chips)
schema_chips = {fieldid:0,photred_field:'',file:'',expnum:'',chip:0,base:'',$
                airmass:0.0,gain:0.0,rdnoise:0.0,nx:0L,ny:0L,wcstype:'',$
                pixscale:0.0,ra:0.0d0,dec:0.0d0,wcsrms:0.0,gaiawcsrms:0.0,$
                gaiawcsnmatch:0L,fwhm:0.0,$
                skymode:0.0,skysig:0.0,dao_nsources:0L,dao_depth:0.0,$
                dao_npsfstars:0L,dao_psftype:'',dao_psfboxsize:0L,$
                dao_psfvarorder:0L,dao_psfchi:0.0,alf_nsources:0L,$
                alf_depth:0.0,apcor:0.0,ebv:0.0,refexpnum:'',$
                vertex_ra1:0.0d0,vertex_ra2:0.0d0,vertex_ra3:0.0d0,$
                vertex_ra4:0.0d0,vertex_dec1:0.0d0,vertex_dec2:0.0d0,$
                vertex_dec3:0.0d0,vertex_dec4:0.0d0,$
                nsrc:0L,calibrated:0B,ubercal_magoffset:0.0,$
                ubercal_flag:0,zpcalib_magoffset:0.0,zpcalib_magofferr:0.0,$
                zpcalibflag:0,photometric:0B,badsoln:0B,calib_band:'',$
                calib_colband:'',calib_colsign:0,calib_zpterm:0.0,$
                calib_zptermsig:0.0,calib_amterm:0.0,calib_amtermsig:0.0,$
                calib_colterm:0.0,calib_coltermsig:0.0}
newchips = replicate(schema_chips,nchips)
STRUCT_ASSIGN,chips,newchips
newchips.fieldid = fieldid
newchips.photred_field = chips.field
newchips.gaiawcsrms = chips.gaiarms
newchips.gaiawcsnmatch = chips.gaianmatch
newchips.vertex_ra1 = chips.vertices_ra[0]
newchips.vertex_ra2 = chips.vertices_ra[1]
newchips.vertex_ra3 = chips.vertices_ra[2]
newchips.vertex_ra4 = chips.vertices_ra[3]
newchips.vertex_dec1 = chips.vertices_dec[0]
newchips.vertex_dec2 = chips.vertices_dec[1]
newchips.vertex_dec3 = chips.vertices_dec[2]
newchips.vertex_dec4 = chips.vertices_dec[3]
newchips.calib_band = chips.band
newchips.calib_colband = chips.colband
newchips.calib_colsign = chips.colsign
newchips.calib_zpterm = chips.zpterm
newchips.calib_zptermsig = chips.zptermsig
newchips.calib_amterm = chips.amterm
newchips.calib_amtermsig = chips.amtermsig
newchips.calib_colterm = chips.colterm
newchips.calib_coltermsig = chips.coltermsig
chips_mjd = dblarr(nchips)  ; need this for newsrc
for j=0,nchips-1 do chips_mjd[j]=date2jd(chips[j].utdate+'T'+chips[j].uttime,/mjd)
MWRFITS,newchips,outdir+ifield+'_chip.fits',/create
; Add the units
head = headfits(outdir+ifield+'_chip.fits',exten=1)
sxaddhist,'',head,/comment
sxaddhist,'  ***  Column units ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TUNIT1','None'
sxaddpar,head,'TUNIT2','None'
sxaddpar,head,'TUNIT3','None'
sxaddpar,head,'TUNIT4','None'
sxaddpar,head,'TUNIT5','None'
sxaddpar,head,'TUNIT6','None'
sxaddpar,head,'TUNIT7','None'
sxaddpar,head,'TUNIT8','Electrons/ADU'
sxaddpar,head,'TUNIT9','Electrons'
sxaddpar,head,'TUNIT10','None'
sxaddpar,head,'TUNIT11','None'
sxaddpar,head,'TUNIT12','None'
sxaddpar,head,'TUNIT13','Arcseconds'
sxaddpar,head,'TUNIT14','Degrees'
sxaddpar,head,'TUNIT15','Degrees'
sxaddpar,head,'TUNIT16','Arcseconds'
sxaddpar,head,'TUNIT17','Arcseconds'
sxaddpar,head,'TUNIT18','None'
sxaddpar,head,'TUNIT19','Pixels'
sxaddpar,head,'TUNIT20','Counts'
sxaddpar,head,'TUNIT21','Counts'
sxaddpar,head,'TUNIT22','None'
sxaddpar,head,'TUNIT23','Magnitude'
sxaddpar,head,'TUNIT24','None'
sxaddpar,head,'TUNIT25','None'
sxaddpar,head,'TUNIT26','Pixels'
sxaddpar,head,'TUNIT27','None'
sxaddpar,head,'TUNIT28','None'
sxaddpar,head,'TUNIT29','None'
sxaddpar,head,'TUNIT30','Magnitude'
sxaddpar,head,'TUNIT31','Magnitude'
sxaddpar,head,'TUNIT32','Magnitude'
sxaddpar,head,'TUNIT33','None'
sxaddpar,head,'TUNIT34','Degrees'
sxaddpar,head,'TUNIT35','Degrees'
sxaddpar,head,'TUNIT36','Degrees'
sxaddpar,head,'TUNIT37','Degrees'
sxaddpar,head,'TUNIT38','Degrees'
sxaddpar,head,'TUNIT39','Degrees'
sxaddpar,head,'TUNIT40','Degrees'
sxaddpar,head,'TUNIT41','Degrees'
sxaddpar,head,'TUNIT42','None'
sxaddpar,head,'TUNIT43','None'
sxaddpar,head,'TUNIT44','Magnitude'
sxaddpar,head,'TUNIT45','None'
sxaddpar,head,'TUNIT46','Magnitude'
sxaddpar,head,'TUNIT47','Magnitude'
sxaddpar,head,'TUNIT48','None'
sxaddpar,head,'TUNIT49','None'
sxaddpar,head,'TUNIT50','None'
sxaddpar,head,'TUNIT51','None'
sxaddpar,head,'TUNIT52','None'
sxaddpar,head,'TUNIT53','None'
sxaddpar,head,'TUNIT54','Magnitude'
sxaddpar,head,'TUNIT55','Magnitude'
sxaddpar,head,'TUNIT56','Magnitude'
sxaddpar,head,'TUNIT57','Magnitude'
sxaddpar,head,'TUNIT58','None'
sxaddpar,head,'TUNIT59','None'

; Column UCD
sxaddhist,'',head,/comment
sxaddhist,'  ***  Column UCD ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TUCD1','meta.id'
sxaddpar,head,'TUCD2','meta.id'
sxaddpar,head,'TUCD3','meta.file'
sxaddpar,head,'TUCD4','meta.id'
sxaddpar,head,'TUCD5','meta.id'
sxaddpar,head,'TUCD6','meta.id;meta.main'
sxaddpar,head,'TUCD7','obs.airMass'
sxaddpar,head,'TUCD8','instr.det.gain'
sxaddpar,head,'TUCD9','instr.det.noise'
sxaddpar,head,'TUCD10','meta.number'
sxaddpar,head,'TUCD11','meta.number'
sxaddpar,head,'TUCD12','meta.code'
sxaddpar,head,'TUCD13','instr.pixel'
sxaddpar,head,'TUCD14','pos.eq.ra;meta.main'
sxaddpar,head,'TUCD15','pos.eq.dec;meta.main'
sxaddpar,head,'TUCD16','stat.stdev'
sxaddpar,head,'TUCD17','stat.stdev'
sxaddpar,head,'TUCD18','meta.number'
sxaddpar,head,'TUCD19','instr.det.psf'
sxaddpar,head,'TUCD20',''
sxaddpar,head,'TUCD21','stat.error'
sxaddpar,head,'TUCD22','meta.number'
sxaddpar,head,'TUCD23',''
sxaddpar,head,'TUCD24','meta.number'
sxaddpar,head,'TUCD25','meta.code'
sxaddpar,head,'TUCD26','meta.number'
sxaddpar,head,'TUCD27','meta.code'
sxaddpar,head,'TUCD28','stat.fit.chi2'
sxaddpar,head,'TUCD29','meta.number'
sxaddpar,head,'TUCD30',''
sxaddpar,head,'TUCD31',''
sxaddpar,head,'TUCD32','phys.absorption'
sxaddpar,head,'TUCD33','meta.id'
sxaddpar,head,'TUCD34','pos.eq.ra'
sxaddpar,head,'TUCD35','pos.eq.ra'
sxaddpar,head,'TUCD36','pos.eq.ra'
sxaddpar,head,'TUCD37','pos.eq.ra'
sxaddpar,head,'TUCD38','pos.eq.dec'
sxaddpar,head,'TUCD39','pos.eq.dec'
sxaddpar,head,'TUCD40','pos.eq.dec'
sxaddpar,head,'TUCD41','pos.eq.dec'
sxaddpar,head,'TUCD42','meta.number'
sxaddpar,head,'TUCD43','meta.code'
sxaddpar,head,'TUCD44','phot.mag'
sxaddpar,head,'TUCD45','stat.error;phot.mag'
sxaddpar,head,'TUCD46','arith.zp;phot.mag'
sxaddpar,head,'TUCD47','stat.error;phot.mag'
sxaddpar,head,'TUCD48','meta.code;phot.mag'
sxaddpar,head,'TUCD49','meta.code'
sxaddpar,head,'TUCD50','meta.code'
sxaddpar,head,'TUCD51','instr.filter'
sxaddpar,head,'TUCD52','instr.filter'
sxaddpar,head,'TUCD53','meta.code'
sxaddpar,head,'TUCD54','arith.zp;phot.mag'
sxaddpar,head,'TUCD55','stat.error;phot.mag'
sxaddpar,head,'TUCD56','phys.absorption'
sxaddpar,head,'TUCD57','stat.error'
sxaddpar,head,'TUCD58',''
sxaddpar,head,'TUCD59','stat.error'

; Column descriptions
sxaddhist,'',head,/comment
sxaddhist,'  ***  Column descriptions ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TCOMM1','SMASH Field ID'
sxaddpar,head,'TCOMM2','Short "field" name for the PHOTRED reduction'
sxaddpar,head,'TCOMM3','Filename for this chip'
sxaddpar,head,'TCOMM4','Eight-digit exposure number'
sxaddpar,head,'TCOMM5','Chip number (1-62)'
sxaddpar,head,'TCOMM6','Base of all files for this chip'
sxaddpar,head,'TCOMM7','Airmass'
sxaddpar,head,'TCOMM8','Gain, electrons/ADU'
sxaddpar,head,'TCOMM9','Read noise, electrons'
sxaddpar,head,'TCOMM10','Number of pixels in X-dimension'
sxaddpar,head,'TCOMM11','Number of pixels in Y-dimension'
sxaddpar,head,'TCOMM12','WCS type, always TPV'
sxaddpar,head,'TCOMM13','Arcsec per pixel'
sxaddpar,head,'TCOMM14','Right ascension at center of chip (in degrees)'
sxaddpar,head,'TCOMM15','Declination at center of chip (in degrees)'
sxaddpar,head,'TCOMM16','RMS in WCS fit with respect to reference catalog (normally USNO-B1) in arcsec'
sxaddpar,head,'TCOMM17','The RMS of the residuals using the GAIA reference stars.'
sxaddpar,head,'TCOMM18','GAIA matches used for the astrometric fitting.'
sxaddpar,head,'TCOMM19','FWHM of PSF (i.e. seeing) in pixels (multiply by 0.26 to get arcsec)'
sxaddpar,head,'TCOMM20','Median sky background level'
sxaddpar,head,'TCOMM21','Sigma in sky background'
sxaddpar,head,'TCOMM22','Total number of DAOPHOT/ALLSTAR sources in chip'
sxaddpar,head,'TCOMM23','ALLSTAR "depth" (peak of the source histogram), instrumental mags'
sxaddpar,head,'TCOMM24','Number of stars used to construct PSF'
sxaddpar,head,'TCOMM25','Analytical function type for DAOPHOT PSF'
sxaddpar,head,'TCOMM26','PSF box size in pixels'
sxaddpar,head,'TCOMM27','Spatial PSF variations, -1 to 2'
sxaddpar,head,'TCOMM28','DAOPHOT PSF "chi" "goodness-of-fit" value'
sxaddpar,head,'TCOMM29','Total number of DAOPHOT/ALLFRAMEs sources. -1 if ALLFRAME was not run on this exposure'
sxaddpar,head,'TCOMM30','Median ALLFRAME depth (peak of the source histogram), instrumental mags, NaN if ALLFRAME not run'
sxaddpar,head,'TCOMM31','Aperture correction'
sxaddpar,head,'TCOMM32','Median SFD E(B-V) for this chip'
sxaddpar,head,'TCOMM33','Exposure number for the reference frame.'
sxaddpar,head,'TCOMM34','First RA value for the four corners of the chip'
sxaddpar,head,'TCOMM35','Second RA value for the four corners of the chip'
sxaddpar,head,'TCOMM36','Third RA value for the four corners of the chip'
sxaddpar,head,'TCOMM37','Fourth RA value for the four corners of the chip'
sxaddpar,head,'TCOMM38','First DEC value for the four corners of the chip'
sxaddpar,head,'TCOMM39','Second DEC value for the four corners of the chip'
sxaddpar,head,'TCOMM40','Third DEC value for the four corners of the chip'
sxaddpar,head,'TCOMM41','Fourth DEC value for the four corners of the chip'
sxaddpar,head,'TCOMM42','Number of sources in SOURCE table for this chip'
sxaddpar,head,'TCOMM43','Is this chip calibrated or not.'
sxaddpar,head,'TCOMM44','The magnitude offset applied from the ubercal algorithm'
sxaddpar,head,'TCOMM45','Ubercal flag: 0-none, 1-good, 2-used med exp+chip offset, 3-used med exp+std chip offset.'
sxaddpar,head,'TCOMM46','The zero-point offset used for this chip.'
sxaddpar,head,'TCOMM47','The uncertainty in the zero-point offset.'
sxaddpar,head,'TCOMM48','ZP calib type: 1-phot DECam data, 2-overlap, 3-0.9m data, 4-SMASH-GAIA color-color relations.'
sxaddpar,head,'TCOMM49','Is this from a photometric night'
sxaddpar,head,'TCOMM50','Did we have a good transformation equation/solution for these data'
sxaddpar,head,'TCOMM51','Transformation equation band name'
sxaddpar,head,'TCOMM52','Transformation equation band name to construct the color'
sxaddpar,head,'TCOMM53','Trans. eqn. sign: 1 - color=band-colband; -1 - color=colband-band.'
sxaddpar,head,'TCOMM54','Transformation equation zeropoint term'
sxaddpar,head,'TCOMM55','Transformation equation uncertainty in zeropoint term'
sxaddpar,head,'TCOMM56','Transformation equation airmass/extinction term'
sxaddpar,head,'TCOMM57','Transformation equation uncertainty in airmass/extinction term'
sxaddpar,head,'TCOMM58','Transformation equation color term'
sxaddpar,head,'TCOMM59','Transformation equation uncertainty in color term'

MODFITS,outdir+ifield+'_chip.fits',0,head,exten_no=1

;TTYPE1  = 'FIELDID '           /                                                
;TTYPE2  = 'PHOTRED_FIELD'      /                                                
;TTYPE3  = 'FILE    '           /                                                
;TTYPE4  = 'EXPNUM  '           /                                                
;TTYPE5  = 'CHIP    '           /                                                
;TTYPE6  = 'BASE    '           /                                                
;TTYPE7  = 'AIRMASS '           /                                                
;TTYPE8  = 'GAIN    '           /                                                
;TTYPE9  = 'RDNOISE '           /                                                
;TTYPE10 = 'NX      '           /                                                
;TTYPE11 = 'NY      '           /                                                
;TTYPE12 = 'WCSTYPE '           /                                                
;TTYPE13 = 'PIXSCALE'           /                                                
;TTYPE14 = 'RA      '           /                                                
;TTYPE15 = 'DEC     '           /                                                
;TTYPE16 = 'WCSRMS  '           /                                                
;TTYPE17 = 'GAIAWCSRMS'         /                                                
;TTYPE18 = 'GAIAWCSNMATCH'      /                                                
;TTYPE19 = 'FWHM    '           /                                                
;TTYPE20 = 'SKYMODE '           /                                                
;TTYPE21 = 'SKYSIG  '           /                                                
;TTYPE22 = 'DAO_NSOURCES'       /                                                
;TTYPE23 = 'DAO_DEPTH'          /                                                
;TTYPE24 = 'DAO_NPSFSTARS'      /                                                
;TTYPE25 = 'DAO_PSFTYPE'        /                                                
;TTYPE26 = 'DAO_PSFBOXSIZE'     /                                                
;TTYPE27 = 'DAO_PSFVARORDER'    /                                                
;TTYPE28 = 'DAO_PSFCHI'         /                                                
;TTYPE29 = 'ALF_NSOURCES'       /                                                
;TTYPE30 = 'ALF_DEPTH'          /                                                
;TTYPE31 = 'APCOR   '           /                                                
;TTYPE32 = 'EBV     '           /                                                
;TTYPE33 = 'REFEXPNUM'          /                                                
;TTYPE34 = 'VERTEX_RA1'         /                                                
;TTYPE35 = 'VERTEX_RA2'         /                                                
;TTYPE36 = 'VERTEX_RA3'         /                                                
;TTYPE37 = 'VERTEX_RA4'         /                                                
;TTYPE38 = 'VERTEX_DEC1'        /                                                
;TTYPE39 = 'VERTEX_DEC2'        /                                                
;TTYPE40 = 'VERTEX_DEC3'        /                                                
;TTYPE41 = 'VERTEX_DEC4'        /                                                
;TTYPE42 = 'NSRC    '           /                                                
;TTYPE43 = 'CALIBRATED'         /                                                
;TTYPE44 = 'UBERCAL_MAGOFFSET'  /                                                
;TTYPE45 = 'UBERCAL_FLAG'       /                                                
;TTYPE46 = 'ZPCALIB_MAGOFFSET'  /                                                
;TTYPE47 = 'ZPCALIB_MAGOFFERR'  /                                                
;TTYPE48 = 'ZPCALIBFLAG'        /                                                
;TTYPE49 = 'PHOTOMETRIC'        /                                                
;TTYPE50 = 'BADSOLN '           /                                                
;TTYPE51 = 'CALIB_BAND'         /                                                
;TTYPE52 = 'CALIB_COLBAND'      /                                                
;TTYPE53 = 'CALIB_COLSIGN'      /                                                
;TTYPE54 = 'CALIB_ZPTERM'       /                                                
;TTYPE55 = 'CALIB_ZPTERMSIG'    /                                                
;TTYPE56 = 'CALIB_AMTERM'       /                                                
;TTYPE57 = 'CALIB_AMTERMSIG'    /                                                
;TTYPE58 = 'CALIB_COLTERM'      /                                                
;TTYPE59 = 'CALIB_COLTERMSIG'   /  

; ADD COMMENTS KEYWORDS!!!

; --- exposures ---
; expnum, nchips, filter, exptime, utdate, uttime, airmass,
; wcstype, ra, dec, wcsrms, fwhm, skymode, skysig, dao_nsources,
; dao_depth, dao_psfchi, alf_nsources, alf_depth, apcor,
; ebv, magname, photometric, badsoln
;  -strip out MAGNAME?
;  -add DATEOBS, JD
;  -change NANs and 99.99 to something else??
exp = mrdfits(dir+ifield+'_combined_exposures.fits.gz',1)
nexp = n_elements(exp)
schema_exp = {expnum:'',nchips:0L,fieldid:0,filter:'',exptime:0.0,$
              dateobs:'',mjd:0.0d0,night_mjd:0L,airmass:0.0,wcstype:'',ra:0.0d0,$
              dec:0.0d0,wcsrms:0.0,fwhm:0.0,skymode:0.0,skysig:0.0,$
              dao_nsources:0L,dao_depth:0.0,dao_psfchi:0.0,alf_nsources:0L,$
              alf_depth:0.0,apcor:0.0,ebv:0.0,photometric:0B,badsoln:0B}
newexp = replicate(schema_exp,nexp)
STRUCT_ASSIGN,exp,newexp
newexp.fieldid = fieldid
newexp.dateobs = exp.utdate+'T'+exp.uttime
for j=0,nexp-1 do newexp[j].mjd=date2jd(newexp[j].dateobs,/mjd)
MATCH,newexp.expnum,chips.expnum,ind1,ind2,/sort
newexp[ind1].night_mjd = chips[ind2].mjd
MWRFITS,newexp,outdir+ifield+'_exposure.fits',/create
; Add the units
head = headfits(outdir+ifield+'_exposure.fits',exten=1)
sxaddhist,'',head,/comment
sxaddhist,'  ***  Column units ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TUNIT1','None'
sxaddpar,head,'TUNIT2','None'
sxaddpar,head,'TUNIT3','None'
sxaddpar,head,'TUNIT4','None'
sxaddpar,head,'TUNIT5','Seconds'
sxaddpar,head,'TUNIT6','None'
sxaddpar,head,'TUNIT7','Days'
sxaddpar,head,'TUNIT8','Days'
sxaddpar,head,'TUNIT9','None'
sxaddpar,head,'TUNIT10','None'
sxaddpar,head,'TUNIT11','Degrees'
sxaddpar,head,'TUNIT12','Degrees'
sxaddpar,head,'TUNIT13','Arcseconds'
sxaddpar,head,'TUNIT14','Pixels'
sxaddpar,head,'TUNIT15','Counts'
sxaddpar,head,'TUNIT16','Counts'
sxaddpar,head,'TUNIT17','None'
sxaddpar,head,'TUNIT18','Magnitude'
sxaddpar,head,'TUNIT19','None'
sxaddpar,head,'TUNIT20','None'
sxaddpar,head,'TUNIT21','Magnitude'
sxaddpar,head,'TUNIT22','Magnitude'
sxaddpar,head,'TUNIT23','Magnitude'
sxaddpar,head,'TUNIT24','None'
sxaddpar,head,'TUNIT25','None'

; Add UCD
sxaddhist,'',head,/comment
sxaddhist,'  ***  Column UCD ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TUCD1','meta.id;meta.main'
sxaddpar,head,'TUCD2','meta.number'
sxaddpar,head,'TUCD3','meta.id'
sxaddpar,head,'TUCD4','instr.filter'
sxaddpar,head,'TUCD5',''
sxaddpar,head,'TUCD6','time.start'
sxaddpar,head,'TUCD7','time.start'
sxaddpar,head,'TUCD8',''
sxaddpar,head,'TUCD9','obs.airMass'
sxaddpar,head,'TUCD10','meta.code'
sxaddpar,head,'TUCD11','pos.eq.ra;meta.main'
sxaddpar,head,'TUCD12','pos.eq.dec;meta.main'
sxaddpar,head,'TUCD13','stat.stdev'
sxaddpar,head,'TUCD14','instr.det.psf'
sxaddpar,head,'TUCD15',''
sxaddpar,head,'TUCD16','stat.stdev'
sxaddpar,head,'TUCD17','meta.number'
sxaddpar,head,'TUCD18',''
sxaddpar,head,'TUCD19','stat.fit.chi2'
sxaddpar,head,'TUCD20','meta.number'
sxaddpar,head,'TUCD21',''
sxaddpar,head,'TUCD22',''
sxaddpar,head,'TUCD23','phys.absorption'
sxaddpar,head,'TUCD24','meta.code'
sxaddpar,head,'TUCD25','meta.code'

; Add comments
sxaddhist,'',head,/comment
sxaddhist,'  ***  Column descriptions ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TCOMM1','Eight-digit exposure number'
sxaddpar,head,'TCOMM2','Number of chips with good data for this exposure, normally 60 or 61'
sxaddpar,head,'TCOMM3','SMASH Field ID'
sxaddpar,head,'TCOMM4','Filter, u, g, r, i, or z'
sxaddpar,head,'TCOMM5','Exposure time in seconds'
sxaddpar,head,'TCOMM6','Observation timestamp: YYYY-MM-DDTHH:MM:SS.SSSSSS'
sxaddpar,head,'TCOMM7','Observation Modified Julian Date'
sxaddpar,head,'TCOMM8','Night Modified Julian Date (integer)'
sxaddpar,head,'TCOMM9','Airmass'
sxaddpar,head,'TCOMM10','WCS type, always TPV'
sxaddpar,head,'TCOMM11','Right ascension at center of frame (in degrees)'
sxaddpar,head,'TCOMM12','Declination at center of frame (in degrees)'
sxaddpar,head,'TCOMM13','RMS in WCS fit with respect to reference catalog (normally USNO-B1) in arcsec, median across chips'
sxaddpar,head,'TCOMM14','FWHM of PSF (i.e. seeing) in pixels (multiply by 0.26 to get arcsec), median across chips'
sxaddpar,head,'TCOMM15','Median sky background level across chips'
sxaddpar,head,'TCOMM16','Median sigma in sky background across chips'
sxaddpar,head,'TCOMM17','Total number of DAOPHOT/ALLSTAR sources in frame'
sxaddpar,head,'TCOMM18','Median ALLSTAR "depth" (peak of the source histogram) across all chips, instrumental mags'
sxaddpar,head,'TCOMM19','Median DAOPHOT PSF "chi" value across all chips'
sxaddpar,head,'TCOMM20','Total number of DAOPHOT/ALLFRAMEs sources in frame. -1 if ALLFRAME was not run on this exposure'
sxaddpar,head,'TCOMM21','Median ALLFRAME depth (peak of the source histogram) across all chips, instrumental mags'
sxaddpar,head,'TCOMM22','Aperture correction'
sxaddpar,head,'TCOMM23','Median SFD E(B-V) across frame'
sxaddpar,head,'TCOMM24','Is this from a photometric night'
sxaddpar,head,'TCOMM25','Did we have a good transformation equation/solution for these data'
MODFITS,outdir+ifield+'_exposure.fits',0,head,exten_no=1

;COMMENT  *** Column names ***                                                   
;COMMENT                                                                         
;TTYPE1  = 'EXPNUM  '           /                                                
;TTYPE2  = 'NCHIPS  '           /                                                
;TTYPE3  = 'FIELDID '           /                                                
;TTYPE4  = 'FILTER  '           /                                                
;TTYPE5  = 'EXPTIME '           /                                                
;TTYPE6  = 'DATEOBS '           /                                                
;TTYPE7  = 'MJD     '           /                                                
;TTYPE8  = 'NIGHT_MJD'          /                                                
;TTYPE9  = 'AIRMASS '           /                                                
;TTYPE10 = 'WCSTYPE '           /                                                
;TTYPE11 = 'RA      '           /                                                
;TTYPE12 = 'DEC     '           /                                                
;TTYPE13 = 'WCSRMS  '           /                                                
;TTYPE14 = 'FWHM    '           /                                                
;TTYPE15 = 'SKYMODE '           /                                                
;TTYPE16 = 'SKYSIG  '           /                                                
;TTYPE17 = 'DAO_NSOURCES'       /                                                
;TTYPE18 = 'DAO_DEPTH'          /                                                
;TTYPE19 = 'DAO_PSFCHI'         /                                                
;TTYPE20 = 'ALF_NSOURCES'       /                                                
;TTYPE21 = 'ALF_DEPTH'          /                                                
;TTYPE22 = 'APCOR   '           /                                                
;TTYPE23 = 'EBV     '           /                                                
;TTYPE24 = 'PHOTOMETRIC'        /                                                
;TTYPE25 = 'BADSOLN '           /   

; ADD COMMENTS KEYWORDS!!!

; --- allsrc ---
; cmbindx, chipindx, fid, id, idref, x, y, xref, yref, mag, err,
; cmag, cerr, chip, sharp, flag, prob, ra, dec, raref, decref
;  -strip out cmbindx, chipindx
;  -add MJD, EXPNUM, and unique ID (maybe FIELDID but with 'Field' removed)
;  -add unique chip information
;  -change NANs and 99.99 to something else?

; NEED TO ADD:
; -raerr, decerr, forced

allsrc = mrdfits(dir+ifield+'_combined_allsrc.fits.gz',1)
nallsrc = n_elements(allsrc)
schema_allsrc = {id:'',origid:0L,refid:0L,fieldid:0,expnum:'',chip:0,mjd:0.0d0,filter:'',x:0.0,y:0.0,xref:0.0,yref:0.0,$
                 forced:0.0,mag:0.0,err:0.0,cmag:0.0,cerr:0.0,chi:0.0,sharp:0.0,$
                 flag:0,prob:0.0,ra:0.0d0,dec:0.0d0,raerr:0.0,decerr:0.0,raindiv:0.0d0,decindiv:0.0d0,raref:0.0d0,$
                 decref:0.0d0}
newsrc = replicate(schema_allsrc,nallsrc)
STRUCT_ASSIGN,allsrc,newsrc
newsrc.id = strmid(strtrim(allsrc.fid,2),5)  ; strip 'Field' portion
newsrc.fieldid = fieldid
newsrc.origid = allsrc.id         ; original allstar/allframe ID
newsrc.refid = allsrc.idref
newsrc.expnum = chips[allsrc.chipindx].expnum
newsrc.chip = chips[allsrc.chipindx].chip
newsrc.mjd = chips_mjd[allsrc.chipindx]
newsrc.filter = chips[allsrc.chipindx].filter
MWRFITS,newsrc,outdir+ifield+'_source.fits',/create
undefine,allsrc    ; free up memory
; Add the units
head = headfits(outdir+ifield+'_source.fits',exten=1)
sxaddhist,'',head,/comment
sxaddhist,'  ***  Column units ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TUNIT1','None'
sxaddpar,head,'TUNIT2','None'
sxaddpar,head,'TUNIT3','None'
sxaddpar,head,'TUNIT4','None'
sxaddpar,head,'TUNIT5','None'
sxaddpar,head,'TUNIT6','None'
sxaddpar,head,'TUNIT7','Days'
sxaddpar,head,'TUNIT8','None'
sxaddpar,head,'TUNIT9','Pixels'
sxaddpar,head,'TUNIT10','Pixels'
sxaddpar,head,'TUNIT11','Pixels'
sxaddpar,head,'TUNIT12','Pixels'
sxaddpar,head,'TUNIT13','None'
sxaddpar,head,'TUNIT14','Magnitude'
sxaddpar,head,'TUNIT15','Magnitude'
sxaddpar,head,'TUNIT16','Magnitude'
sxaddpar,head,'TUNIT17','Magnitude'
sxaddpar,head,'TUNIT18','None'
sxaddpar,head,'TUNIT19','None'
sxaddpar,head,'TUNIT20','None'
sxaddpar,head,'TUNIT21','None'
sxaddpar,head,'TUNIT22','Degrees'
sxaddpar,head,'TUNIT23','Degrees'
sxaddpar,head,'TUNIT24','Degrees'
sxaddpar,head,'TUNIT25','Degrees'
sxaddpar,head,'TUNIT26','Degrees'
sxaddpar,head,'TUNIT27','Degrees'
sxaddpar,head,'TUNIT28','Degrees'
sxaddpar,head,'TUNIT29','Degrees'

; Add UCD
sxaddhist,'',head,/comment
sxaddhist,'  ***  Column UCD ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TUCD1',''
sxaddpar,head,'TUCD2',''
sxaddpar,head,'TUCD3',''
sxaddpar,head,'TUCD4',''
sxaddpar,head,'TUCD5',''
sxaddpar,head,'TUCD6',''
sxaddpar,head,'TUCD7',''
sxaddpar,head,'TUCD8',''
sxaddpar,head,'TUCD9',''
sxaddpar,head,'TUCD10',''
sxaddpar,head,'TUCD11',''
sxaddpar,head,'TUCD12',''
sxaddpar,head,'TUCD13',''
sxaddpar,head,'TUCD14',''
sxaddpar,head,'TUCD15',''
sxaddpar,head,'TUCD16',''
sxaddpar,head,'TUCD17',''
sxaddpar,head,'TUCD18',''
sxaddpar,head,'TUCD19',''
sxaddpar,head,'TUCD20',''
sxaddpar,head,'TUCD21',''
sxaddpar,head,'TUCD22',''
sxaddpar,head,'TUCD23',''
sxaddpar,head,'TUCD24',''
sxaddpar,head,'TUCD25',''
sxaddpar,head,'TUCD26',''
sxaddpar,head,'TUCD27',''
sxaddpar,head,'TUCD28',''
sxaddpar,head,'TUCD29',''

; Add descriptions
sxaddhist,'',head,/comment
sxaddhist,'  ***  Column descriptions ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TCOMM1','Unique ID for this source, the field name plus a running number'
sxaddpar,head,'TCOMM2','ID used for this detection in the original als/alf chip file'
sxaddpar,head,'TCOMM3','ID used for this detection in the reference frame file'
sxaddpar,head,'TCOMM4','SMASH Field ID'
sxaddpar,head,'TCOMM5','Eight-digit exposure number'
sxaddpar,head,'TCOMM6','Chip number (1-62)'
sxaddpar,head,'TCOMM7','Exposure Modified Julian Date'
sxaddpar,head,'TCOMM8','Filter, u, g, r, i, or z'
sxaddpar,head,'TCOMM9','X-coordinate for this source in the original chip (1-indexed)'
sxaddpar,head,'TCOMM10','Y-coordinate for this source in the original chip (1-indexed)'
sxaddpar,head,'TCOMM11','X-coordinate for this source in the reference chip (1-indexed)'
sxaddpar,head,'TCOMM12','Y-coordinate for this source in the reference chip (1-indexed)'
sxaddpar,head,'TCOMM13','Boolean flag indicating if this is forced photometry (ALLFRAME) or single-frame (ALLSTAR)'
sxaddpar,head,'TCOMM14','Instrumental magnitude from ALLFRAME or ALLSTAR (not both)'
sxaddpar,head,'TCOMM15','Uncertainty of MAG'
sxaddpar,head,'TCOMM16','Calibrated magnitude version of MAG'
sxaddpar,head,'TCOMM17','Uncertainty of CMAG including errors in calibration'
sxaddpar,head,'TCOMM18','DAOPHOT chi value, i.e. how well the PSF fit this source'
sxaddpar,head,'TCOMM19','DAOPHOT sharp value, measurement of peakiness'
sxaddpar,head,'TCOMM20','Source Extractor FLAG value from coadd image (only if ALLFRAME was run otherwise -1)'
sxaddpar,head,'TCOMM21','Source Extractor stellaricity probability value (0~galaxy, 1~star, -1 if ALLFRAME not run)'
sxaddpar,head,'TCOMM22','Right Ascension (J2000.0) of source, in degrees'
sxaddpar,head,'TCOMM23','Declination (J2000.0) of source, in degrees'
sxaddpar,head,'TCOMM24',''
sxaddpar,head,'TCOMM25',''
sxaddpar,head,'TCOMM26',''
sxaddpar,head,'TCOMM27',''
sxaddpar,head,'TCOMM28',''
sxaddpar,head,'TCOMM29',''

MODFITS,outdir+ifield+'_source.fits',0,head,exten_no=1

;TTYPE1  = 'ID      '           /                                                
;TTYPE2  = 'ORIGID  '           /                                                
;TTYPE3  = 'REFID   '           /                                                
;TTYPE4  = 'FIELDID '           /                                                
;TTYPE5  = 'EXPNUM  '           /                                                
;TTYPE6  = 'CHIP    '           /                                                
;TTYPE7  = 'MJD     '           /                                                
;TTYPE8  = 'FILTER  '           /                                                
;TTYPE9  = 'X       '           /                                                
;TTYPE10 = 'Y       '           /                                                
;TTYPE11 = 'XREF    '           /                                                
;TTYPE12 = 'YREF    '           /                                                
;TTYPE13 = 'FORCED  '           /                                                
;TTYPE14 = 'MAG     '           /                                                
;TTYPE15 = 'ERR     '           /                                                
;TTYPE16 = 'CMAG    '           /                                                
;TTYPE17 = 'CERR    '           /                                                
;TTYPE18 = 'CHI     '           /                                                
;TTYPE19 = 'SHARP   '           /                                                
;TTYPE20 = 'FLAG    '           /                                                
;TTYPE21 = 'PROB    '           /                                                
;TTYPE22 = 'RA      '           /                                                
;TTYPE23 = 'DEC     '           /                                                
;TTYPE24 = 'RAERR  '            /                                                
;TTYPE25 = 'DECERR '            /       
;TTYPE26 = 'RAINDIV  '          /                                                
;TTYPE27 = 'DECINDIV '          /       
;TTYPE28 = 'RAREF   '           /                                                
;TTYPE29 = 'DECREF  '           /       

; ADD COMMENTS KEYWORDS!!!

; --- allobj ---
; id, ra, dec, rascatter, decscatter, ndet, depthflag, srcindx,
; srcfindx, u, uerr, uscatter, g, gerr, gscatter, r, rerr,
; rscatter, i, ierr, iscatter, z, zerr, zscatter, chi, sharp,
; flag, prob, ebv
;  -strip out srcindx, srcfindx
;  -remove 'Field' from ID
;  -change NANs and 99.99 to something else?

; NEED TO ADD:
; -raerr, decerr

allobj = mrdfits(dir+ifield+'_combined_allobj.fits.gz',1)
nallobj = n_elements(allobj)
schema_allobj = {id:'',fieldid:0,ra:0.0d0,dec:0.0d0,raerr:0.0,decerr:0.0,$
                ndet:0L,depthflag:0B,umag:0.0,uerr:0.0,uscatter:0.0,ndetu:0,gmag:0.0,$
                gerr:0.0,gscatter:0.0,ndetg:0,rmag:0.0,rerr:0.0,rscatter:0.0,ndetr:0L,$
                imag:0.0,ierr:0.0,iscatter:0.0,ndeti:0,zmag:0.0,zerr:0.0,zscatter:0.0,$
                ndetz:0,u_g:99.99,g_r:99.99,g_i:99.99,i_z:99.99,chi:0.0,sharp:0.0,$
                flag:0,prob:0.0,ebv:0.0}
newobj = replicate(schema_allobj,nallobj)
STRUCT_ASSIGN,allobj,newobj
newobj.fieldid = fieldid   ; field integer ID
newobj.id = strmid(strtrim(allobj.id,2),5)  ; strip 'Field' portion
newobj.umag = allobj.u
newobj.gmag = allobj.g
newobj.rmag = allobj.r
newobj.imag = allobj.i
newobj.zmag = allobj.z
; THESE ARE NOW IN THE ORIGINAL FITS CATALOGS
;; put in number of detections of each star per filter
;for i=0,nallobj-1 do begin
;  ;if i mod 50000 eq 0 then print,i
;  srcind = allobj[i].srcindx[0:allobj[i].ndet-1]
;  filter = chips[allsrc[srcind].chipindx].filter
;  MATCH,filter,'u',ind1,ind2,/sort,count=nu
;  MATCH,filter,'g',ind1,ind2,/sort,count=ng
;  MATCH,filter,'r',ind1,ind2,/sort,count=nr
;  MATCH,filter,'i',ind1,ind2,/sort,count=ni
;  MATCH,filter,'z',ind1,ind2,/sort,count=nz
;  newobj[i].ndetu = nu
;  newobj[i].ndetg = ng
;  newobj[i].ndetr = nr
;  newobj[i].ndeti = ni
;  newobj[i].ndetz = nz
;endfor
; colors
gdug = where(newobj.umag lt 50 and newobj.gmag lt 50,ngdug)
if ngdug gt 0 then newobj[gdug].u_g = newobj[gdug].umag - newobj[gdug].gmag
gdgr = where(newobj.gmag lt 50 and newobj.rmag lt 50,ngdgr)
if ngdgr gt 0 then newobj[gdgr].g_r = newobj[gdgr].gmag - newobj[gdgr].rmag
gdgi = where(newobj.gmag lt 50 and newobj.imag lt 50,ngdgi)
if ngdgi gt 0 then newobj[gdgi].g_i = newobj[gdgi].gmag - newobj[gdgi].imag
gdiz = where(newobj.imag lt 50 and newobj.zmag lt 50,ngdiz)
if ngdiz gt 0 then newobj[gdiz].i_z = newobj[gdiz].imag - newobj[gdiz].zmag
MWRFITS,newobj,outdir+ifield+'_object.fits',/create
; Add the units
head = headfits(outdir+ifield+'_object.fits',exten=1)
sxaddhist,'',head,/comment
sxaddhist,'  ***  Column units ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TUNIT1','None'
sxaddpar,head,'TUNIT2','None'
sxaddpar,head,'TUNIT3','Degrees'
sxaddpar,head,'TUNIT4','Degrees'
sxaddpar,head,'TUNIT5','Degrees'
sxaddpar,head,'TUNIT6','Degrees'
sxaddpar,head,'TUNIT7','None'
sxaddpar,head,'TUNIT8','None'
sxaddpar,head,'TUNIT9','Magnitude'
sxaddpar,head,'TUNIT10','Magnitude'
sxaddpar,head,'TUNIT11','Magnitude'
sxaddpar,head,'TUNIT12','None'
sxaddpar,head,'TUNIT13','Magnitude'
sxaddpar,head,'TUNIT14','Magnitude'
sxaddpar,head,'TUNIT15','Magnitude'
sxaddpar,head,'TUNIT16','None'
sxaddpar,head,'TUNIT17','Magnitude'
sxaddpar,head,'TUNIT18','Magnitude'
sxaddpar,head,'TUNIT19','Magnitude'
sxaddpar,head,'TUNIT20','None'
sxaddpar,head,'TUNIT21','Magnitude'
sxaddpar,head,'TUNIT22','Magnitude'
sxaddpar,head,'TUNIT23','Magnitude'
sxaddpar,head,'TUNIT24','None'
sxaddpar,head,'TUNIT25','Magnitude'
sxaddpar,head,'TUNIT26','Magnitude'
sxaddpar,head,'TUNIT27','Magnitude'
sxaddpar,head,'TUNIT28','None'
sxaddpar,head,'TUNIT29','Magnitude'
sxaddpar,head,'TUNIT30','Magnitude'
sxaddpar,head,'TUNIT31','Magnitude'
sxaddpar,head,'TUNIT32','Magnitude'
sxaddpar,head,'TUNIT33','None'
sxaddpar,head,'TUNIT34','None'
sxaddpar,head,'TUNIT35','None'
sxaddpar,head,'TUNIT36','None'
sxaddpar,head,'TUNIT37','Magnitude'
MODFITS,outdir+ifield+'_object.fits',0,head,exten_no=1

;TTYPE1  = 'ID      '           /                                                
;TTYPE2  = 'FIELDID '           /                                                
;TTYPE3  = 'RA      '           /                                                
;TTYPE4  = 'DEC     '           /                                                
;TTYPE5  = 'RASCATTER'          /                                                
;TTYPE6  = 'DECSCATTER'         /                                                
;TTYPE7  = 'NDET    '           /                                                
;TTYPE8  = 'DEPTHFLAG'          /                                                
;TTYPE9  = 'UMAG    '           /                                                
;TTYPE10 = 'UERR    '           /                                                
;TTYPE11 = 'USCATTER'           /                                                
;TTYPE12 = 'NDETU   '           /                                                
;TTYPE13 = 'GMAG    '           /                                                
;TTYPE14 = 'GERR    '           /                                                
;TTYPE15 = 'GSCATTER'           /                                                
;TTYPE16 = 'NDETG   '           /                                                
;TTYPE17 = 'RMAG    '           /                                                
;TTYPE18 = 'RERR    '           /                                                
;TTYPE19 = 'RSCATTER'           /                                                
;TTYPE20 = 'NDETR   '           /                                                
;TTYPE21 = 'IMAG    '           /                                                
;TTYPE22 = 'IERR    '           /                                                
;TTYPE23 = 'ISCATTER'           /                                                
;TTYPE24 = 'NDETI   '           /                                                
;TTYPE25 = 'ZMAG    '           /                                                
;TTYPE26 = 'ZERR    '           /                                                
;TTYPE27 = 'ZSCATTER'           /                                                
;TTYPE28 = 'NDETZ   '           /                                                
;TTYPE29 = 'U_G     '           /                                                
;TTYPE30 = 'G_R     '           /                                                
;TTYPE31 = 'G_I     '           /                                                
;TTYPE32 = 'I_Z     '           /                                                
;TTYPE33 = 'CHI     '           /                                                
;TTYPE34 = 'SHARP   '           /                                                
;TTYPE35 = 'FLAG    '           /                                                
;TTYPE36 = 'PROB    '           /                                                
;TTYPE37 = 'EBV     '           /               

; ADD COMMENTS KEYWORDS!!!!

; --- xmatch ---
; ADD XMATCH!!!!
xmatch = mrdfits(dir+ifield+'_combined_allobj_xmatch.fits.gz',1)
nxmatch = n_elements(xmatch)
;schema_xmatch = {id:'',gaia_match:0b,gaia_matchdist:0.0,gaia_source:0L,gaia_ra:0.0d0,gaia_dec:0.0d0,$
;         gaia_raerr:0.0,gaia_decerr:0.0,gaia_gmag:0.0,gaia_gerr:0.0,tmass_match:0b,tmass_matchdist:0.0,$
;         tmass_id:'',tmass_ra:0.0d0,tmass_dec:0.0d0,tmass_raerr:0.0,tmass_decerr:0.0,tmass_jmag:0.0,$
;         tmass_jerr:0.0,tmass_hmag:0.0,tmass_herr:0.0,tmass_kmag:0.0,tmass_kerr:0.0,tmass_qflg:'',$
;         wise_match:0b,wise_matchdist:0.0,wise_id:'',wise_ra:0.0d0,wise_dec:0.0d0,wise_raerr:0.0,$
;         wise_decerr:0.0,wise_w1mag:0.0,wise_w1err:0.0,wise_w2mag:0.0,wise_w2err:0.0,wise_w3mag:0.0,$
;         wise_w3err:0.0,wise_w4mag:0.0,wise_w4err:0.0,wise_qph:''}
;newxmatch = replicate(schema_xmatch,nxmatch)
;STRUCT_ASSIGN,allobj,new
MWRFITS,field,outdir+ifield+'_xmatch.fits',/create
; Add the units
head = headfits(outdir+ifield+'_xmatch.fits',exten=1)
sxaddhist,'',head,/comment
sxaddhist,'  ***  Column units ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TUNIT1','None'
sxaddpar,head,'TUNIT2','None'
sxaddpar,head,'TUNIT3','None'
sxaddpar,head,'TUNIT4','None'
sxaddpar,head,'TUNIT5','Degrees'
sxaddpar,head,'TUNIT6','Degrees'
sxaddpar,head,'TUNIT7','Degrees'
sxaddpar,head,'TUNIT8','Degrees'
sxaddpar,head,'TUNIT9','Magnitude'
sxaddpar,head,'TUNIT10','Magnitude'
sxaddpar,head,'TUNIT11','None'
sxaddpar,head,'TUNIT12','None'
sxaddpar,head,'TUNIT13','None'
sxaddpar,head,'TUNIT14','Degrees'
sxaddpar,head,'TUNIT15','Degrees'
sxaddpar,head,'TUNIT16','Degrees'
sxaddpar,head,'TUNIT17','Degrees'
sxaddpar,head,'TUNIT18','Magnitude'
sxaddpar,head,'TUNIT19','Magnitude'
sxaddpar,head,'TUNIT20','Magnitude'
sxaddpar,head,'TUNIT21','Magnitude'
sxaddpar,head,'TUNIT22','Magnitude'
sxaddpar,head,'TUNIT23','Magnitude'
sxaddpar,head,'TUNIT24','None'
sxaddpar,head,'TUNIT25','None'
sxaddpar,head,'TUNIT26','None'
sxaddpar,head,'TUNIT27','None'
sxaddpar,head,'TUNIT28','Degrees'
sxaddpar,head,'TUNIT29','Degrees'
sxaddpar,head,'TUNIT30','Degrees'
sxaddpar,head,'TUNIT31','Degrees'
sxaddpar,head,'TUNIT32','Magnitude'
sxaddpar,head,'TUNIT33','Magnitude'
sxaddpar,head,'TUNIT34','Magnitude'
sxaddpar,head,'TUNIT35','Magnitude'
sxaddpar,head,'TUNIT36','Magnitude'
sxaddpar,head,'TUNIT37','Magnitude'
sxaddpar,head,'TUNIT38','Magnitude'
sxaddpar,head,'TUNIT39','Magnitude'
sxaddpar,head,'TUNIT40','None'

; Add the UCD

sxaddhist,'',head,/comment
sxaddhist,'  ***  Column UCD ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TUCD1','meta.id;meta.main'
sxaddpar,head,'TUCD2',''
sxaddpar,head,'TUCD3',''
sxaddpar,head,'TUCD4','meta.id'
sxaddpar,head,'TUCD5','pos.eq.ra'
sxaddpar,head,'TUCD6','pos.eq.dec'
sxaddpar,head,'TUCD7','stat.error'
sxaddpar,head,'TUCD8','stat.error'
sxaddpar,head,'TUCD9','phot.mag'
sxaddpar,head,'TUCD10','stat.error;phot.mag'
sxaddpar,head,'TUCD11',''
sxaddpar,head,'TUCD12',''
sxaddpar,head,'TUCD13','meta.id'
sxaddpar,head,'TUCD14','pos.eq.ra'
sxaddpar,head,'TUCD15','pos.eq.dec'
sxaddpar,head,'TUCD16','stat.error'
sxaddpar,head,'TUCD17','stat.error'
sxaddpar,head,'TUCD18','phot.mag;em.IR.J'
sxaddpar,head,'TUCD19','stat.error;phot.mag;em.IR.J'
sxaddpar,head,'TUCD20','phot.mag;em.IR.H'
sxaddpar,head,'TUCD21','stat.error;phot.mag;em.IR.H'
sxaddpar,head,'TUCD22','phot.mag;em.IR.K'
sxaddpar,head,'TUCD23','stat.error;phot.mag;em.IR.K'
sxaddpar,head,'TUCD24','meta.code.qual;phot'
sxaddpar,head,'TUCD25',''
sxaddpar,head,'TUCD26',''
sxaddpar,head,'TUCD27','meta.id'
sxaddpar,head,'TUCD28','pos.eq.ra'
sxaddpar,head,'TUCD29','pos.eq.dec'
sxaddpar,head,'TUCD30','stat.error'
sxaddpar,head,'TUCD31','stat.error'
sxaddpar,head,'TUCD32','phot.mag;em.IR.3-4um'
sxaddpar,head,'TUCD33','stat.error;phot.mag'
sxaddpar,head,'TUCD34','phot.mag;em.IR.4-8um'
sxaddpar,head,'TUCD35','stat.error;phot.mag'
sxaddpar,head,'TUCD36','phot.mag;em.IR.8-15um'
sxaddpar,head,'TUCD37','stat.error;phot.mag'
sxaddpar,head,'TUCD38','phot.mag;em.IR.15-30um'
sxaddpar,head,'TUCD39','stat.error;phot.mag'
sxaddpar,head,'TUCD40','meta.code.qual;phot'

; Add the descriptions
sxaddhist,'',head,/comment
sxaddhist,'  ***  Column descriptions ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TCOMM1','SMASH OBJECT source ID.'
sxaddpar,head,'TCOMM2','Boolean flag indicating if a match with GAIA DR1 was found.'
sxaddpar,head,'TCOMM3','Distance to matched GAIA-DR1 source (arcseconds).'
sxaddpar,head,'TCOMM4','GAIA source ID.'
sxaddpar,head,'TCOMM5','GAIA right ascension (degrees).'
sxaddpar,head,'TCOMM6','GAIA declination (degrees)'
sxaddpar,head,'TCOMM7','GAIA RA uncertainty (arcseconds)'
sxaddpar,head,'TCOMM8','GAIA DEC uncertainty (arcseconds)'
sxaddpar,head,'TCOMM9','GAIA G-band mean magnitude'
sxaddpar,head,'TCOMM10','GAIA uncertainty in G-band magnitue (calculated from flux uncertainty)'
sxaddpar,head,'TCOMM11','Boolean flag indicating if a match with 2MASS-PSC was found'
sxaddpar,head,'TCOMM12','Distance to matched 2MASS-PSC source (arcseconds)'
sxaddpar,head,'TCOMM13','2MASS source ID'
sxaddpar,head,'TCOMM14','2MASS right ascension (degrees)'
sxaddpar,head,'TCOMM15','2MASS declination (degrees)'
sxaddpar,head,'TCOMM16','2MASS RA uncertainty (arcseconds)'
sxaddpar,head,'TCOMM17','2MASS DEC uncertainty (arcseconds)'
sxaddpar,head,'TCOMM18','2MASS J-band magnitude'
sxaddpar,head,'TCOMM19','2MASS uncertainty in J-band magnitude'
sxaddpar,head,'TCOMM20','2MASS H-band magnitude'
sxaddpar,head,'TCOMM21','2MASS uncertainty in H-band magnitude'
sxaddpar,head,'TCOMM22','2MASS Ks-band magnitude'
sxaddpar,head,'TCOMM23','2MASS uncertainty in Ks-band magnitude'
sxaddpar,head,'TCOMM24','2MASS JHKs photometry quality flag'
sxaddpar,head,'TCOMM25','Boolean flag indicating if a match with ALLWISE was found'
sxaddpar,head,'TCOMM26','Distance to matched ALLWISE source (arcseconds)'
sxaddpar,head,'TCOMM27','ALLWISE source ID'
sxaddpar,head,'TCOMM28','ALLWISE right ascension (degrees)'
sxaddpar,head,'TCOMM29','ALLWISE declination (degrees)'
sxaddpar,head,'TCOMM30','ALLWISE RA uncertainty (arcseconds)'
sxaddpar,head,'TCOMM31','ALLWISE DEC uncertainty (arcseconds)'
sxaddpar,head,'TCOMM32','ALLWISE W1 magnitude (3.35um)'
sxaddpar,head,'TCOMM33','ALLWISE uncertainty in W1 magnitude'
sxaddpar,head,'TCOMM34','ALLWISE W2 magnitude (4.6um)'
sxaddpar,head,'TCOMM35','ALLWISE uncertainty in W2 magnitude'
sxaddpar,head,'TCOMM36','ALLWISE W3 magnitude (11.6um)'
sxaddpar,head,'TCOMM37','ALLWISE uncertainty in W3 magnitude'
sxaddpar,head,'TCOMM38','ALLWISE W4 magnitude (22.1um)'
sxaddpar,head,'TCOMM39','ALLWISE uncertainty in W4 magnitude'
sxaddpar,head,'TCOMM40','ALLWISE photometric quality flag'

;MODFITS,outdir+ifield+'_xmatch.fits',0,head,exten_no=1

stop

; --- fields ---
;   which bands are calibrated
;   name, central RA/DEC, nexposures
field = {fieldid:0,name:'',ra:0.0d0,dec:0.0d0,glon:0.0d0,glat:0.0d0,$
         mslon:0.0d0,mslat:0.0d0,nexp:0L,nexp_u:0L,nexp_g:0L,nexp_r:0L,$
         nexp_i:0L,nexp_z:0L,nchips:0L,nsrc:0LL,nobj:0LL,ucalib:0B,gcalib:0B,$
         rcalib:0B,icalib:0B,zcalib:0B}
field.fieldid = fieldid
field.name = ifield
field.nexp = n_elements(exp)
field.nchips = n_elements(chips)
field.nsrc = nallsrc
field.nobj = nallobj
uexp = where(newexp.filter eq 'u',nuexp)
field.nexp_u = nuexp
gexp = where(newexp.filter eq 'g',ngexp)
field.nexp_g = ngexp
rexp = where(newexp.filter eq 'r',nrexp)
field.nexp_r = nrexp
iexp = where(newexp.filter eq 'i',niexp)
field.nexp_i = niexp
zexp = where(newexp.filter eq 'z',nzexp)
field.nexp_z = nzexp
; which bands are calibrated
uchip = where(chips.filter eq 'u',nuchips)
ucalib = where(chips.filter eq 'u' and chips.calibrated eq 1 and chips.photometric eq 1 and chips.badsoln eq 0,nucalib)
ucalibfrac = nucalib/float(nuchips) 
if nucalib gt 0 then field.ucalib=1 else field.ucalib=0
;if ucalibfrac gt 0.7 then field.ucalib=1 else field.ucalib=0
gchip = where(chips.filter eq 'g',ngchips)
gcalib = where(chips.filter eq 'g' and chips.calibrated eq 1 and chips.photometric eq 1 and chips.badsoln eq 0,ngcalib)
gcalibfrac = ngcalib/float(ngchips)
if ngcalib gt 0 then field.gcalib=1 else field.gcalib=0
;if gcalibfrac gt 0.70 then field.gcalib=1 else field.gcalib=0
rchip = where(chips.filter eq 'r',nrchips)
rcalib = where(chips.filter eq 'r' and chips.calibrated eq 1 and chips.photometric eq 1 and chips.badsoln eq 0,nrcalib)
rcalibfrac = nrcalib/float(nrchips)
if nrcalib gt 0 then field.rcalib=1 else field.rcalib=0
;if rcalibfrac gt 0.70 then field.rcalib=1 else field.rcalib=0
ichip = where(chips.filter eq 'i',nichips)
icalib = where(chips.filter eq 'i' and chips.calibrated eq 1 and chips.photometric eq 1 and chips.badsoln eq 0,nicalib)
icalibfrac = nicalib/float(nichips)
if nicalib gt 0 then field.icalib=1 else field.icalib=0
;if icalibfrac gt 0.70 then field.icalib=1 else field.icalib=0
zchip = where(chips.filter eq 'z',nzchips)
zcalib = where(chips.filter eq 'z' and chips.calibrated eq 1 and chips.photometric eq 1 and chips.badsoln eq 0,nzcalib)
zcalibfrac = nzcalib/float(nzchips)
if nzcalib gt 0 then field.zcalib=1 else field.zcalib=0
;if zcalibfrac gt 0.70 then field.zcalib=1 else field.zcalib=0
; get coordinates from smash_fields_final.txt file
MATCH,'Field'+strtrim(fieldstr.num,2),ifield,ind1,ind2,/sort,count=nmatch
field.ra = fieldstr[ind1[0]].radeg
field.dec = fieldstr[ind1[0]].dedeg
glactc,field.ra,field.dec,2000.0,glon,glat,1,/deg
field.glon = glon
field.glat = glat
field.mslon = fieldstr[ind1[0]].mslon
field.mslat = fieldstr[ind1[0]].mslat
MWRFITS,field,outdir+ifield+'_field.fits',/create
; Add the units
head = headfits(outdir+ifield+'_field.fits',exten=1)
sxaddhist,'',head,/comment
sxaddhist,'  ***  Column units ***',head,/comment
sxaddhist,'',head,/comment
sxaddpar,head,'TUNIT1','None'
sxaddpar,head,'TUNIT2','None'
sxaddpar,head,'TUNIT3','Degrees'
sxaddpar,head,'TUNIT4','Degrees'
sxaddpar,head,'TUNIT5','Degrees'
sxaddpar,head,'TUNIT6','Degrees'
sxaddpar,head,'TUNIT7','Degrees'
sxaddpar,head,'TUNIT8','Degrees'
sxaddpar,head,'TUNIT9','None'
sxaddpar,head,'TUNIT10','None'
sxaddpar,head,'TUNIT11','None'
sxaddpar,head,'TUNIT12','None'
sxaddpar,head,'TUNIT13','None'
sxaddpar,head,'TUNIT14','None'
sxaddpar,head,'TUNIT15','None'
sxaddpar,head,'TUNIT16','None'
sxaddpar,head,'TUNIT17','None'
sxaddpar,head,'TUNIT18','None'
sxaddpar,head,'TUNIT19','None'
sxaddpar,head,'TUNIT20','None'
sxaddpar,head,'TUNIT21','None'
sxaddpar,head,'TUNIT22','None'
MODFITS,outdir+ifield+'_field.fits',0,head,exten_no=1

;TTYPE1  = 'FIELDID '           /                                                
;TTYPE2  = 'NAME    '           /                                                
;TTYPE3  = 'RA      '           /                                                
;TTYPE4  = 'DEC     '           /                                                
;TTYPE5  = 'GLON    '           /                                                
;TTYPE6  = 'GLAT    '           /                                                
;TTYPE7  = 'MSLON   '           /                                                
;TTYPE8  = 'MSLAT   '           /                                                
;TTYPE9  = 'NEXP    '           /                                                
;TTYPE10 = 'NEXP_U  '           /                                                
;TTYPE11 = 'NEXP_G  '           /                                                
;TTYPE12 = 'NEXP_R  '           /                                                
;TTYPE13 = 'NEXP_I  '           /                                                
;TTYPE14 = 'NEXP_Z  '           /                                                
;TTYPE15 = 'NCHIPS  '           /                                                
;TTYPE16 = 'NSRC    '           /                                                
;TTYPE17 = 'NOBJ    '           /                                                
;TTYPE18 = 'UCALIB  '           /                                                
;TTYPE19 = 'GCALIB  '           /                                                
;TTYPE20 = 'RCALIB  '           /                                                
;TTYPE21 = 'ICALIB  '           /                                                
;TTYPE22 = 'ZCALIB  '           / 

;stop

end
