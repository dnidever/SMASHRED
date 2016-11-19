;+
;
; SMASHRED_LOAD_CATPHOT
;
; This program loads all of the photometry for one PHOTRED catalog
;
; INPUTS:
;  info        The structure with the relevant input information needed
;                to load the data.
;  /useast     Use the .ast files, by default the .phot files are used.
;  /useorig    Use the original als/alf files.
;  =reduxdir   The reduction directory, the default is "/data/smash/cp/red/photred/"
;  =outputdir  The output directory, the default is reduxdir+"catalogs/final/"
;  /redo       Reget the photometry even if the temporary output file
;                already exists.
;
; OUTPUTS:
;  chstr       The structure with information on each chip file and
;                indices to the photometry information in ALLSRC.
;  allsrc      All source data concatenated together.
;  =error      The error message if one occurred.
;
; USAGE:
;  IDL>smashred_load_catphot,info,chstr,/useast
;
; By D.Nidever April 2016
;-

pro smashred_load_catphot,info,chstr,allsrc,useast=useast,useorig=useorig,reduxdir=reduxdir,redo=redo,error=error,$
                          outputdir=outputdir

undefine,error
undefine,chstr
undefine,allsrc

; Checking the inputs
if n_elements(info) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_load_catphot,info,chstr,allsrc,useast=useast,useorig=useorig,reduxdir=reduxdir,redo=redo,error=error'
  return
endif

; Defaults
if n_elements(reduxdir) eq 0 then reduxdir='/data/smash/cp/red/photred/'
if file_test(reduxdir,/directory) eq 0 then begin
  error = reduxdir+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif
if n_elements(outputdir) eq 0 then outputdir=reduxdir+'catalogs/final/'
if file_test(outputdir,/directory) eq 0 then begin
  if not keyword_set(silent) then print,outputdir+' does NOT exist.  Creating it.'
  FILE_MKDIR,outputdir
endif
; Temporary directory
tmpdir = outputdir+'/tmp/'
if file_test(tmpdir,/directory) eq 0 then FILE_MKDIR,tmpdir

; Load the list of bad exposures
badexp = IMPORTASCII('~/projects/SMASHRED/obslog/smash_badexposures.txt',/header,delim=string(9B),$
                     fieldtypes=[7,7,7,7],/silent)
nbadexp = n_elements(badexp)

; Construct the base name
fbase = file_basename(info.file,'_summary.fits')  ; the observed field name

; Output filename
outfile = tmpdir+fbase+'_'+info.night+'_photred.fits'


; --- Initalize ALLSRC structure ---
; ALLSRC structure schema
;  cmbindx       index into ALLOBJ, added later on in smashred_crossmatch.pro
;  fid           final ID, added later on in smashred_crossmatch.pro
;  id            original ID in ALS/ALF file
;  idref         ID in PHOT/AST file
;  x/y           original X/Y coordinate in the chip image
;  xref/yref     final X/Y coordinates in the reference frame
;                  NOT the original frame
;  mag/err       instrumental photometry from PHOT/AST
;  cmag/cerr     calibrated photometry that will be added later
;  chi/sharp     original chi/sharp from ALS/ALF file
;  ra/dec        ra/dec coordinates using X/Y and chip WCS with GAIA references
;  raindiv/decindiv  ra/dec coordinates using X/Y and chip WCS with
;                        USNO-B1/2MASS as reference
;  raref/decref  ra/dec coordinates using xref/yref and ref image WCS
;                        with USNO-B1/2MASS as reference
allsrc_schema = {cmbindx:-1L,chipindx:-1L,fid:'',id:-1L,idref:-1L,x:0.0,y:0.0,xref:0.0,yref:0.0,mag:0.0,err:0.0,$
                 cmag:-1.0,cerr:-1.0,chi:0.0,sharp:0.0,flag:-1,prob:-1.0,ra:0.0d0,dec:0.0d0,$
                 raindiv:0.0d0,decindiv:0.0d0,raref:0.0d0,decref:0.0d0}
;                 raerr:0.0,decerr:0.0,raindiv:0.0d0,decindiv:0.0d0,raref:0.0d0,decref:0.0d0}
allsrc = replicate(allsrc_schema,5000000L)
nallsrc = n_elements(allsrc)
cur_allsrc_indx = 0LL

; Load in the data
If file_test(outfile) eq 0 or keyword_set(redo) then begin

  ; Load the chip-level summary information
  chstr = MRDFITS(info.file,2,/silent)  ; load the chip structure
  nchstr = n_elements(chstr)
  night = info.night
  field = chstr[0].field

  ; --- Remove bad exposures ---
  undefine,badexpind,ubadexpnum
  for k=0,nbadexp-1 do begin
    MATCH,chstr.expnum,badexp[k].expnum,badind1,badind2,/sort,count=nbadexp
    if nbadexp gt 0 then begin
      push,badexpind,badind1
      push,ubadexpnum,badexp[k].expnum
    endif
  endfor
  nbadexpind = n_elements(badexpind)
  if nbadexpind gt 0 then begin
    print,'REMOVING ',strtrim(n_elements(ubadexpnum),2),' BAD exposure(s): ',strjoin(ubadexpnum,',  ')
    if nbadexpind eq nchstr then begin
      print,'NO EXPOSURES LEFT FOR THIS PHOTRED CATALOG'
      undefine,chstr,allsrc,fstr
      return
    endif else begin
      REMOVE,badexpind,chstr
    endelse
    nchstr = n_elements(chstr)
    ; Also remove from exposures FSTR structure
    fstr = *info.fstr
    MATCH,fstr.expnum,ubadexpnum,badind1,badind2,/sort,count=nbadexp2
    REMOVE,badind1,fstr
    info.nexp = n_elements(fstr)
    uibands = uniq(fstr.filter,sort(fstr.filter))
    expbands = strjoin(fstr[uibands].filter)
    info.bands = expbands
    info.fstr = ptr_new(fstr)
  endif

  ; KLUDGE! Removing DUPLICATE exposures in short and deep fields
  if night eq '20150318' and field eq 'F5' then begin  ; Field130sh
    MATCH,chstr.expnum,'00423440',bdind,dum,/sort,count=nbdind
    if nbdind gt 0 then begin
      print,'KLUDGE! Removing DUPLICATE exposure 00423440 (in short and deep field) from short field'
      REMOVE,bdind,chstr
      nchstr = n_elements(chstr)
    endif
  endif
  if night eq '20150330' and field eq 'F7' then begin  ; Field130sh
    MATCH,chstr.expnum,'00426607',bdind,dum,/sort,count=nbdind
    if nbdind gt 0 then begin
      print,'KLUDGE! Removing DUPLICATE exposure 00423440 (in short and deep field) from short field'
      REMOVE,bdind,chstr
      nchstr = n_elements(chstr)
    endif
  endif

  ; Add some tags to CHSTR
  add_tag,chstr,'night','',chstr
  chstr.night = night
  add_tag,chstr,'alftiletype','',chstr
  add_tag,chstr,'gaiarms',0.0,chstr
  add_tag,chstr,'gaianmatch',0L,chstr
  add_tag,chstr,'refexpnum','',chstr
  add_tag,chstr,'vertices_ra',dblarr(4),chstr
  add_tag,chstr,'vertices_dec',dblarr(4),chstr
  add_tag,chstr,'nsrc',-1L,chstr
  add_tag,chstr,'allsrcindx',-1LL,chstr

  ; --- Get data from the phot/ast files ----
  ; Use PHOT files
  if not keyword_set(useast) and not keyword_set(useorig) then begin
    ext = '.phot'
    photfiles = file_search(reduxdir+night+'/'+field+'/'+field+'-*_??'+ext,count=nphotfiles)
    print,strtrim(nphotfiles,2),' PHOT files for ',fbase,' in ',reduxdir+night+'/'+field+'/'
  ; Use AST files
  endif else begin
    ext = '.ast'
    photfiles = file_search(reduxdir+night+'/'+field+'/'+field+'-*_??'+ext,count=nphotfiles)
    print,strtrim(nphotfiles,2),' AST files for ',fbase,' in ',reduxdir+night+'/'+field+'/'
  endelse

  ; --- Loop through individual chip AST/PHOT files ---
  for i=0,nphotfiles-1 do begin
  
    print,strtrim(i+1,2),' ',photfiles[i]
    ; Load the PHOT/AST file
    phot = IMPORTASCII(photfiles[i],/header,/silent)
    tags = tag_names(phot)
    phbase = file_basename(photfiles[i],ext)
    ichip = long(first_el(strsplit(phbase,'_',/extract),/last))
    dum = first_el(strsplit(phbase,'-',/extract),/last)
    refexpnum = first_el(strsplit(dum,'_',/extract))  ; the reference exposure

    ; Load the MCH file to get the ordering of the exposures
    mchfile = file_dirname(photfiles[i])+'/'+file_basename(photfiles[i],ext)+'.mch'
    if file_test(mchfile) eq 0 then stop,mchfile,' NOT FOUND'
    LOADMCH,mchfile,mchfilelist
    mchfilelist = file_basename(mchfilelist,'.als')  ; remove .als ending

    chind = where(chstr.chip eq ichip,nchind)

    ; Check the combination tile type 
    if tag_exist(phot,'PROB') then begin
      combfile = reduxdir+night+'/'+field+'/'+phbase+'_comb.fits'
      combhead = HEADFITS(combfile)
      alftiletype = sxpar(combhead,'AFTILTYP',count=nalftiletype)
      if nalftiletype gt 0 then chstr[chind].alftiletype=alftiletype else chstr[chind].alftiletype='ORIG'
    endif

    ; Loop through all of the exposures for this chip
    for j=0,nchind-1 do begin
      ; Figure out which magnitude column we should use
      ;  PHOT: use CALIB_MAGNAME from chip structure
      if not keyword_set(useast) then begin
        mind = where(tags eq strtrim(chstr[chind[j]].calib_magname,2),nmind)
      ;  AST: use order in the MCH file
      endif else begin
        mchfilelist_index = where(mchfilelist eq chstr[chind[j]].base,nmchfilelist_index)
        if nmchfilelist_index eq 0 then stop,chstr[chind[j]].base+' NOT FOUND in '+mchfile
        mind = where(tags eq 'MAG'+strtrim(mchfilelist_index[0]+1,2),nmind)
      endelse

      ; Find sources in PHOT/AST file with good measurements for this exp
      gdphot = where(phot.(mind[0]) lt 50,nsrc)
      if nsrc eq 0 then begin
        stop,'No good sources in this chip'
        goto,nogoodsrc
      endif
      phot1 = phot[gdphot]

      ; Initialize and fill in the final structure
      src = replicate(allsrc_schema,nsrc)
      STRUCT_ASSIGN,phot1,src,/nozero
      src.idref = phot1.id
      src.xref = phot1.x
      src.yref = phot1.y
      src.raref = phot1.ra
      src.decref = phot1.dec
      src.mag = phot1.(mind[0])
      src.err = phot1.(mind[0]+1)   ; assume error is the next column
      src.chipindx = chind[j]
      ;  These values will be filled in below from orig file
      src.id = -1
      src.x = !values.f_nan
      src.y = !values.f_nan
      src.ra = !values.f_nan
      src.dec = !values.f_nan
      src.chi = !values.f_nan
      src.sharp = !values.f_nan

      ; Load the FITS header of the original chip file
      fitsfile = reduxdir+night+'/'+field+'/'+strtrim(chstr[chind[j]].base,2)+'.fits'
      if file_test(fitsfile) eq 0 then stop,fitsfile,' NOT FOUND'
      head = headfits(fitsfile)
      ; Get the MJD for this exposure
      mjd = PHOTRED_GETMJD(fitsfile,'ctio')


      ; LOAD the ORIGINAL ALF/ALS files to use the original X/Y/RA/DEC/CHI/SHARP
      ;--------------------------------------------------------------------------

      ; Are we using ALF or ALS files?
      ;   Determine if als/alf if prob/flag is there
      if tag_exist(phot,'PROB') then origphotexten='alf' else origphotexten='als'
      origphotfile = reduxdir+night+'/'+field+'/'+strtrim(chstr[chind[j]].base,2)+'.'+origphotexten
      if file_test(origphotfile) eq 0 then stop,origphotfile,' NOT FOUND'

      ; Load original photometry als/alf file
      LOADALS,origphotfile,origphot,origphothead,count=norigphot

      ; ALF: The star IDs are uniform across the alf files
      ;       makes the matching easy
      if origphotexten eq 'alf' then begin
        MATCH,src.idref,origphot.id,ind1,ind2,/sort,count=nmatch
      ; ALS: Need to use the TFR file to get indices for the als file
      endif else begin
        tfrfile = reduxdir+night+'/'+field+'/'+phbase+'.tfr'
        ; Check if the there is a .tfr.orig file, this is the DAOMASTER
        ;  tfr file if ALLFRAME was run
        ;if file_test(tfrfile+'.orig') eq 1 then tfrfile+='.orig'
        if file_test(tfrfile+'.orig') eq 1 then stop,'TFR.ORIG FOUND!  There could be problems'
        if file_test(tfrfile) eq 0 then stop,tfrfile+' NOT FOUND'
        LOADTFR,tfrfile,tfrlist,tfrstr
        ; What number is it in the frame list
        tfrlist_index = where(tfrlist eq file_basename(origphotfile),ntfrlist_index)          
        if ntfrlist_index eq 0 then stop,file_basename(origphotfile)+' NOT in '+tfrfile
        ; Get indices for the original ALS file, 0-based
        ;   # of stars in TFR file should be same as RAW/AST
        ind2 = reform( tfrstr[gdphot].index[tfrlist_index]-1 )
        ind1 = lindgen(nsrc)
      endelse

      ; Put in final structure
      src[ind1].id = origphot[ind2].id
      src[ind1].x = origphot[ind2].x
      src[ind1].y = origphot[ind2].y
      src[ind1].chi = origphot[ind2].chi
      src[ind1].sharp = origphot[ind2].sharp

      ; Get RA/DEC coordinates for X/Y
      ;  use IDL X/Y convention, starting at (0,)
      HEAD_XYAD,head,src.x-1,src.y-1,ra,dec,/degree
      ;src.ra = ra
      ;src.dec = dec
      src.raindiv = ra
      src.decindiv = dec

      ; Get GAIA-calibrated coordinates
      ;--------------------------------
      if mjd lt 57691 then begin
        gaiawcsfile = reduxdir+night+'/'+field+'/'+strtrim(chstr[chind[j]].base,2)+'.gaiawcs.head'
        if file_test(gaiawcsfile) eq 0 then stop,gaiawcsfile+' NOT FOUND'
        READLINE,gaiawcsfile,gaiahead
      endif else begin
        ; Used GAIA as the default reference catalog from 57691 onward
        gaiahead = head
      endelse
      HEAD_XYAD,gaiahead,src.x-1,src.y-1,gra,gdec,/degree
      src.ra = gra
      src.dec = gdec
      ; Get GAIA RMS and NMATCH
      ;  HISTORY WCSFIT: RMS=0.027 arcsec on Mon Sep 26 03:09:50 2016                    
      rmsind = first_el(where(stregex(gaiahead,'WCSFIT: RMS',/boolean) eq 1,nrmsind),/last)  ; last one
      wcsline = gaiahead[rmsind[0]]
      lo = strpos(wcsline,'RMS=')
      tmp = strmid(wcsline,lo+4)
      gaiarms = float( first_el(strsplit(tmp,' ',/extract)) )
      ;  HISTORY WCSFIT: NMATCH=269 
      nmatchind = first_el(where(stregex(gaiahead,'WCSFIT: NMATCH',/boolean) eq 1,n_nmatchind),/last)
      wcsline = gaiahead[nmatchind[0]]
      lo = strpos(wcsline,'NMATCH=')
      tmp = strmid(wcsline,lo+7)
      gaianmatch = long( first_el(strsplit(tmp,' ',/extract)) )


      ; Compute uncertainty in RA/DEC based on FWHM and S/N
      ;----------------------------------------------------
      ;snr = 1.087/src.err
      ;coorderr = 0.644*chstr[chind[j]].fwhm*chstr[chind[j]].pixscale/snr
      ;;  Might need to modify the uncertainties for ALLFRAME
      ;if (origphotexten eq 'alf') and (n_elements(mchfilelist) ne nchind) then begin
      ;  ; For ALLFRAME the actual error is ~coorderr/sqrt(Nexp) which
      ;  ; will get properly reduced in the OBJECT stage when the
      ;  ; uncertainty in the mean gets calculated.
      ;  ; But if some exposures are bad and thrown out then we need to
      ;  ; correct the individual uncertainties
      ;  coorderr *= sqrt(nchind)/sqrt(n_elements(mchfilelist))
      ;endif
      ;; Add GAIARMS as a uncertainty floor in quadrature
      ;coorderr = sqrt(coorderr^2 + gaiarms^2)
      ;src.raerr = coorderr  ; for now RA/DEC errors are the same
      ;src.decerr = coorderr

      ; Remove bad data for DECam chip 31
      ;----------------------------------
      if (ichip eq 31) and (mjd gt 56660) then begin
        print,'Removing bad data for DECam chip 31'

        ; Remove bad measurements
        ; X: 1-1024 okay
        ; X: 1025-2049 bad
        bdind = where(src.x gt 1024,nbdind,comp=gdind,ncomp=ngdind)
        if nbdind gt 0 then begin   ; some bad ones found
          if ngdind eq 0 then begin   ; all bad
            print,'NO useful measurements in ',fitsfile
            undefine,src
            nsrc = 0
          endif else begin
            print,'Removing '+strtrim(nbdind,2)+' bad measurements, '+strtrim(ngdind,2)+' left.'
            REMOVE,bdind,src
            nsrc = n_elements(src)
          endelse
       endif  ; some bad ones to remove
      endif  ; chip 31

      NOGOODSRC:

      ; Stuff some information in the chip structure
      chstr[chind[j]].gaiarms = gaiarms
      chstr[chind[j]].gaianmatch = gaianmatch
      chstr[chind[j]].refexpnum = refexpnum
      chstr[chind[j]].nsrc = nsrc
      chstr[chind[j]].allsrcindx = cur_allsrc_indx
      ; Get astrometric vertices from header
      nx = sxpar(head,'NAXIS1')
      ny = sxpar(head,'NAXIS2')
      HEAD_XYAD,gaiahead,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/degree
      chstr[chind[j]].vertices_ra = vra
      chstr[chind[j]].vertices_dec = vdec

      ; Add new elements to ALLSRC
      if nallsrc lt cur_allsrc_indx+nsrc then begin
        print,'Adding new elements to ALLSRC'
        new = replicate(allsrc_schema,nallsrc+5000000L)  ; add another 5 million elements
        new[0:nallsrc-1] = allsrc
        allsrc = new
        undefine,new
        nallsrc = n_elements(allsrc)  
      endif

      ; Add SRC to ALLSRC structure
      if nsrc gt 0 then allsrc[cur_allsrc_indx:cur_allsrc_indx+nsrc-1] = src

      ; Increment index
      cur_allsrc_indx += nsrc

      print,j+1,chstr[chind[j]].expnum,chstr[chind[j]].chip,nsrc,format='(I4,A12,I5,I10)'
    endfor
  endfor

  ; Pruning extra ALLSRC elements
  if nallsrc gt cur_allsrc_indx then allsrc = allsrc[0:cur_allsrc_indx-1]

  ; Save the output file
  print,'Saving info temporarily to ',outfile
  MWRFITS,chstr,outfile,/create
  MWRFITS,allsrc,outfile,/silent

; Load the existing temporary file
Endif else begin
  print,'Restoring previously saved ',outfile
  chstr = MRDFITS(outfile,1,/silent)
  allsrc = MRDFITS(outfile,2,/silent)

  ; --- Remove bad exposures from fstr ---
  fstr = *info.fstr
  MATCH,fstr.expnum,badexp.expnum,badind1,badind2,/sort,count=nbadexp
  if nbadexp gt 0 then begin
    REMOVE,badind1,fstr
    info.nexp = n_elements(fstr)
    uibands = uniq(fstr.filter,sort(fstr.filter))
    expbands = strjoin(fstr[uibands].filter)
    info.bands = expbands
    info.fstr = ptr_new(fstr)
  endif

Endelse

;stop

end
