pro fits2db_conv,ifield

; Convert the SMASH binary fits files to the fits
; binary format needed to load the database

dir = '/data/smash/cp/red/photred/catalogs/final/'
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
chips = mrdfits(dir+ifield+'_combined_chips.fits.gz',1)
nchips = n_elements(chips)
schema_chips = {fieldid:0,photred_field:'',file:'',expnum:'',chip:0,base:'',$
                airmass:0.0,gain:0.0,rdnoise:0.0,nx:0L,ny:0L,wcstype:'',$
                pixscale:0.0,ra:0.0d0,dec:0.0d0,wcsrms:0.0,fwhm:0.0,$
                skymode:0.0,skysig:0.0,dao_nsources:0L,dao_depth:0.0,$
                dao_npsfstars:0L,dao_psftype:'',dao_psfboxsize:0L,$
                dao_psfvarorder:0L,dao_psfchi:0.0,alf_nosurces:0L,$
                alf_depth:0.0,apcor:0.0,ebv:0.0,refexpnum:'',$
                vertex_ra1:0.0d0,vertex_ra2:0.0d0,vertex_ra3:0.0d0,$
                vertex_ra4:0.0d0,vertex_dec1:0.0d0,vertex_dec2:0.0d0,$
                vertex_dec3:0.0d0,vertex_dec4:0.0d0,$
                nsrc:0L,calibrated:0B,ubercal_magoffset:0.0,$
                ubercal_flag:0,photometric:0B,badsoln:0B,calib_band:'',$
                calib_colband:'',calib_colsign:0,calib_zpterm:0.0,$
                calib_zptermsig:0.0,calib_amterm:0.0,calib_amtermsig:0.0,$
                calib_colterm:0.0,calib_coltermsig:0.0}
newchips = replicate(schema_chips,nchips)
STRUCT_ASSIGN,chips,newchips
newchips.fieldid = fieldid
newchips.photred_field = chips.field
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
sxaddpar,head,'TUNIT17','Pixels'
sxaddpar,head,'TUNIT18','Counts'
sxaddpar,head,'TUNIT19','Counts'
sxaddpar,head,'TUNIT20','None'
sxaddpar,head,'TUNIT21','Magnitude'
sxaddpar,head,'TUNIT22','None'
sxaddpar,head,'TUNIT23','None'
sxaddpar,head,'TUNIT24','Pixels'
sxaddpar,head,'TUNIT25','None'
sxaddpar,head,'TUNIT26','None'
sxaddpar,head,'TUNIT27','None'
sxaddpar,head,'TUNIT28','Magnitude'
sxaddpar,head,'TUNIT29','Magnitude'
sxaddpar,head,'TUNIT30','Magnitude'
sxaddpar,head,'TUNIT31','None'
sxaddpar,head,'TUNIT32','Degrees'
sxaddpar,head,'TUNIT33','Degrees'
sxaddpar,head,'TUNIT34','Degrees'
sxaddpar,head,'TUNIT35','Degrees'
sxaddpar,head,'TUNIT36','Degrees'
sxaddpar,head,'TUNIT37','Degrees'
sxaddpar,head,'TUNIT38','Degrees'
sxaddpar,head,'TUNIT39','Degrees'
sxaddpar,head,'TUNIT40','None'
sxaddpar,head,'TUNIT41','None'
sxaddpar,head,'TUNIT42','Magnitude'
sxaddpar,head,'TUNIT43','None'
sxaddpar,head,'TUNIT44','None'
sxaddpar,head,'TUNIT45','None'
sxaddpar,head,'TUNIT46','None'
sxaddpar,head,'TUNIT47','None'
sxaddpar,head,'TUNIT48','None'
sxaddpar,head,'TUNIT49','Magnitude'
sxaddpar,head,'TUNIT50','Magnitude'
sxaddpar,head,'TUNIT51','Magnitude'
sxaddpar,head,'TUNIT52','Magnitude'
sxaddpar,head,'TUNIT53','None'
sxaddpar,head,'TUNIT54','None'
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
;TTYPE17 = 'FWHM    '           /                                                
;TTYPE18 = 'SKYMODE '           /                                                
;TTYPE19 = 'SKYSIG  '           /                                                
;TTYPE20 = 'DAO_NSOURCES'       /                                                
;TTYPE21 = 'DAO_DEPTH'          /                                                
;TTYPE22 = 'DAO_NPSFSTARS'      /                                                
;TTYPE23 = 'DAO_PSFTYPE'        /                                                
;TTYPE24 = 'DAO_PSFBOXSIZE'     /                                                
;TTYPE25 = 'DAO_PSFVARORDER'    /                                                
;TTYPE26 = 'DAO_PSFCHI'         /                                                
;TTYPE27 = 'ALF_NOSURCES'       /                                                
;TTYPE28 = 'ALF_DEPTH'          /                                                
;TTYPE29 = 'APCOR   '           /                                                
;TTYPE30 = 'EBV     '           /                                                
;TTYPE31 = 'REFEXPNUM'          /                                                
;TTYPE32 = 'VERTEX_RA1'         /                                                
;TTYPE33 = 'VERTEX_RA2'         /                                                
;TTYPE34 = 'VERTEX_RA3'         /                                                
;TTYPE35 = 'VERTEX_RA4'         /                                                
;TTYPE36 = 'VERTEX_DEC1'        /                                                
;TTYPE37 = 'VERTEX_DEC2'        /                                                
;TTYPE38 = 'VERTEX_DEC3'        /                                                
;TTYPE39 = 'VERTEX_DEC4'        /                                                
;TTYPE40 = 'NSRC    '           /                                                
;TTYPE41 = 'CALIBRATED'         /                                                
;TTYPE42 = 'UBERCAL_MAGOFFSET'  /                                                
;TTYPE43 = 'UBERCAL_FLAG'       /                                                
;TTYPE44 = 'PHOTOMETRIC'        /                                                
;TTYPE45 = 'BADSOLN '           /                                                
;TTYPE46 = 'CALIB_BAND'         /                                                
;TTYPE47 = 'CALIB_COLBAND'      /                                                
;TTYPE48 = 'CALIB_COLSIGN'      /                                                
;TTYPE49 = 'CALIB_ZPTERM'       /                                                
;TTYPE50 = 'CALIB_ZPTERMSIG'    /                                                
;TTYPE51 = 'CALIB_AMTERM'       /                                                
;TTYPE52 = 'CALIB_AMTERMSIG'    /                                                
;TTYPE53 = 'CALIB_COLTERM'      /                                                
;TTYPE54 = 'CALIB_COLTERMSIG'   /    


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

; --- allsrc ---
; cmbindx, chipindx, fid, id, idref, x, y, xref, yref, mag, err,
; cmag, cerr, chip, sharp, flag, prob, ra, dec, raref, decref
;  -strip out cmbindx, chipindx
;  -add MJD, EXPNUM, and unique ID (maybe FIELDID but with 'Field' removed)
;  -add unique chip information
;  -change NANs and 99.99 to something else?
allsrc = mrdfits(dir+ifield+'_combined_allsrc.fits.gz',1)
nallsrc = n_elements(allsrc)
schema_allsrc = {id:'',origid:0L,refid:0L,fieldid:0,expnum:'',chip:0,mjd:0.0d0,filter:'',x:0.0,y:0.0,xref:0.0,yref:0.0,$
                 mag:0.0,err:0.0,cmag:0.0,cerr:0.0,chi:0.0,sharp:0.0,$
                 flag:0,prob:0.0,ra:0.0d0,dec:0.0d0,raref:0.0d0,$
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
sxaddpar,head,'TUNIT13','Magnitude'
sxaddpar,head,'TUNIT14','Magnitude'
sxaddpar,head,'TUNIT15','Magnitude'
sxaddpar,head,'TUNIT16','Magnitude'
sxaddpar,head,'TUNIT17','None'
sxaddpar,head,'TUNIT18','None'
sxaddpar,head,'TUNIT19','None'
sxaddpar,head,'TUNIT20','None'
sxaddpar,head,'TUNIT21','Degrees'
sxaddpar,head,'TUNIT22','Degrees'
sxaddpar,head,'TUNIT23','Degrees'
sxaddpar,head,'TUNIT24','Degrees'
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
;TTYPE13 = 'MAG     '           /                                                
;TTYPE14 = 'ERR     '           /                                                
;TTYPE15 = 'CMAG    '           /                                                
;TTYPE16 = 'CERR    '           /                                                
;TTYPE17 = 'CHI     '           /                                                
;TTYPE18 = 'SHARP   '           /                                                
;TTYPE19 = 'FLAG    '           /                                                
;TTYPE20 = 'PROB    '           /                                                
;TTYPE21 = 'RA      '           /                                                
;TTYPE22 = 'DEC     '           /                                                
;TTYPE23 = 'RAREF   '           /                                                
;TTYPE24 = 'DECREF  '           /       

; --- allobj ---
; id, ra, dec, rascatter, decscatter, ndet, depthflag, srcindx,
; srcfindx, u, uerr, uscatter, g, gerr, gscatter, r, rerr,
; rscatter, i, ierr, iscatter, z, zerr, zscatter, chi, sharp,
; flag, prob, ebv
;  -strip out srcindx, srcfindx
;  -remove 'Field' from ID
;  -change NANs and 99.99 to something else?
allobj = mrdfits(dir+ifield+'_combined_allobj.fits.gz',1)
nallobj = n_elements(allobj)
schema_allobj = {id:'',fieldid:0,ra:0.0d0,dec:0.0d0,rascatter:0.0,decscatter:0.0,$
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
; put in number of detections of each star per filter
for i=0,nallobj-1 do begin
  ;if i mod 50000 eq 0 then print,i
  srcind = allobj[i].srcindx[0:allobj[i].ndet-1]
  filter = chips[allsrc[srcind].chipindx].filter
  MATCH,filter,'u',ind1,ind2,/sort,count=nu
  MATCH,filter,'g',ind1,ind2,/sort,count=ng
  MATCH,filter,'r',ind1,ind2,/sort,count=nr
  MATCH,filter,'i',ind1,ind2,/sort,count=ni
  MATCH,filter,'z',ind1,ind2,/sort,count=nz
  newobj[i].ndetu = nu
  newobj[i].ndetg = ng
  newobj[i].ndetr = nr
  newobj[i].ndeti = ni
  newobj[i].ndetz = nz
endfor
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
ucalib = where(chips.filter eq 'u' and chips.calibrated eq 1,nucalib)
ucalibfrac = nucalib/float(nuchips) 
if ucalibfrac gt 0.7 then field.ucalib=1 else field.ucalib=0
gchip = where(chips.filter eq 'g',ngchips)
gcalib = where(chips.filter eq 'g' and chips.calibrated eq 1,ngcalib)
gcalibfrac = ngcalib/float(ngchips)
if gcalibfrac gt 0.70 then field.gcalib=1 else field.gcalib=0
rchip = where(chips.filter eq 'r',nrchips)
rcalib = where(chips.filter eq 'r' and chips.calibrated eq 1,nrcalib)
rcalibfrac = nrcalib/float(nrchips)
if rcalibfrac gt 0.70 then field.rcalib=1 else field.rcalib=0
ichip = where(chips.filter eq 'i',nichips)
icalib = where(chips.filter eq 'i' and chips.calibrated eq 1,nicalib)
icalibfrac = nicalib/float(nichips)
if icalibfrac gt 0.70 then field.icalib=1 else field.icalib=0
zchip = where(chips.filter eq 'z',nzchips)
zcalib = where(chips.filter eq 'z' and chips.calibrated eq 1,nzcalib)
zcalibfrac = nzcalib/float(nzchips)
if zcalibfrac gt 0.70 then field.zcalib=1 else field.zcalib=0
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
