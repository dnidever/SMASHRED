pro fits2db_conv_expchip_dr2

;; Add some columns to the exposure files
version = 'v6'
dir = '/dl1/users/dnidever/smash/cp/red/photred/catalogs/final/'+version+'/'
outdir = dir+'db/'
if file_test(outdir,/directory) eq 0 then file_mkdir,outdir

fieldstr = importascii('/d1/users/dnidever/smash/cp/red/photred/check_calibrated_v6.fits',/header)
fieldstr.field = strtrim(fieldstr.field,2)
nfields = n_elements(fieldstr)

fieldid = fix(strmid(strtrim(ifield,2),5))
print,ifield

;; Loop over fields
for i=0,nfields-1 do begin
  ifield = fieldstr[i].field
  fieldid = fix(strmid(strtrim(ifield,2),5))
  expfile = dir+ifield+'_combined_exposures.fits.gz'
  expstr = mrdfits(expfile,1,/silent)
  nexp = n_elements(expstr)
  
  ;CHIP Table
  ;add smash fieldid column?
  ;gaiarms -> gaiawcsrms
  ;gaianmatch -> gaiawcsnmatch
  ;vertices_ra_# -> vertex_ra# ??
  ;vertices_dec_# -> vertex_dec# ??
  ;DROP calib_depth
  ;DROP calib_color
  ;DROP calib_zpterm
  ;DROP calib_amterm
  ;DROP calib_colorterm
  ;DROP calib_magname
  ;band -> calib_band
  ;colband -> calib_colband
  ;colsign -> calib_colsign
  ;zpterm -> calib_zpterm
  ;zptermsig -> calib_zptermsig
  ;amterm -> calib_amterm
  ;amtermsig -> calib_amtermsig
  ;colterm -> calib_colterm
  ;coltermsig -> calib_coltermsig
  ;DROP amcolterm
  ;DROP amcoltermsig
  ;DROP colsqterm
  ;DROP colsqterm  

  chips = mrdfits(dir+ifield+'_combined_chips.fits.gz',1)
  nchips = n_elements(chips)
  schema_chips = {chipid:'',fieldid:0,photred_field:'',file:'',expnum:'',chip:0,base:'',$
                  exptime:0.0,filter:'',dateobs:'',mjd:0.0d0,$
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
  newchips.chipid = strtrim(chips.expnum,2)+'_'+strtrim(chips.chip,2)
  newchips.fieldid = fieldid
  newchips.photred_field = chips.field
  newchips.dateobs = chips.utdate+'T'+chips.uttime
  for j=0,nchips-1 do newchips[j].mjd=date2jd(chips[j].utdate+'T'+chips[j].uttime,/mjd)  
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
  
  ;EXPOSURE table
  ;add smash fieldid column?
  ;add MJD column?
  ;add night_mjd column?
  ;add gaiawcsrms column?
  ;DROP magname

  exp = mrdfits(dir+ifield+'_combined_exposures.fits.gz',1)
  nexp = n_elements(exp)
  schema_exp = {expnum:'',nchips:0L,fieldid:0,filter:'',exptime:0.0,$
                dateobs:'',mjd:0.0d0,night_mjd:0L,airmass:0.0,wcstype:'',ra:0.0d0,$
                dec:0.0d0,wcsrms:0.0,gaiawcsrms:0.0,fwhm:0.0,skymode:0.0,skysig:0.0,$
                dao_nsources:0L,dao_depth:0.0,dao_psfchi:0.0,alf_nsources:0L,$
                alf_depth:0.0,apcor:0.0,ebv:0.0,photometric:0B,badsoln:0B}
  newexp = replicate(schema_exp,nexp)
  STRUCT_ASSIGN,exp,newexp
  newexp.fieldid = fieldid
  newexp.dateobs = exp.utdate+'T'+exp.uttime
  for j=0,nexp-1 do newexp[j].mjd=date2jd(newexp[j].dateobs,/mjd)
  MATCH,newexp.expnum,chips.expnum,ind1,ind2,/sort
  newexp[ind1].night_mjd = chips[ind2].mjd
  for j=0,nexp-1 do begin
    MATCH,newexp[j].expnum,chips.expnum,ind1,ind2,/sort
    newexp[j].gaiawcsrms = median(chips[ind2].gaiarms)
  endfor

  push,allchips,newchips
  push,allexp,newexp
  
  stop
  
endfor

;MWRFITS,allchips,'smash_dr2_chip.fits',/create
;MWRFITS,allexp,'smash_dr2_exposure.fits',/create

stop

end
