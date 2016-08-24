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

pro smashred_load_catphot,info,chstr,allsrc,useast=useast,useorig=useorig,reduxdir=reduxdir,redo=redo,error=error

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


; Construct the base name
fbase = file_basename(info.file,'_summary.fits')  ; the observed field name

; Output filename
outfile = tmpdir+fbase+'_'+info.night+'_photred.fits'


; Initalize ALLSRC structure
allsrc_schema = {cmbindx:-1L,chipindx:-1L,fid:'',id:-1L,x:0.0,y:0.0,mag:0.0,err:0.0,$
                    cmag:-1.0,cerr:-1.0,chi:0.0,sharp:0.0,flag:-1,prob:-1.0,ra:0.0d0,dec:0.0d0}
allsrc = replicate(allsrc_schema,5000000L)
nallsrc = n_elements(allsrc)
cur_allsrc_indx = 0LL

; Load in the data
If file_test(outfile) eq 0 or keyword_set(redo) then begin

  ; Load the chip-level summary information
  chstr = MRDFITS(info.file,2,/silent)  ; load the chip structure
  nchstr = n_elements(chstr)
  night = info.night

  ; Get data from the phot/ast files
  add_tag,chstr,'refexpnum','',chstr
  add_tag,chstr,'vertices_ra',dblarr(4),chstr
  add_tag,chstr,'vertices_dec',dblarr(4),chstr
  add_tag,chstr,'nsrc',-1L,chstr
  add_tag,chstr,'allsrcindx',-1LL,chstr

  ; Use PHOT files
  if not keyword_set(useast) and not keyword_set(useorig) then begin
    ext = '.phot'
    photfiles = file_search(reduxdir+night+'/'+chstr[0].field+'/'+chstr[0].field+'-*_??'+ext,count=nphotfiles)
    print,strtrim(nphotfiles,2),' PHOT files for ',fbase,' in ',reduxdir+night+'/'+chstr[0].field+'/'
  ; Use AST files
  endif else begin
    ext = '.ast'
    photfiles = file_search(reduxdir+night+'/'+chstr[0].field+'/'+chstr[0].field+'-*_??'+ext,count=nphotfiles)
    print,strtrim(nphotfiles,2),' AST files for ',fbase,' in ',reduxdir+night+'/'+chstr[0].field+'/'
  endelse

  ; Loop through individual chip AST/PHOT files
  for i=0,nphotfiles-1 do begin
  
    print,strtrim(i+1,2),' ',photfiles[i]
    phot = IMPORTASCII(photfiles[i],/header,/silent)
    tags = tag_names(phot)
    phbase = file_basename(photfiles[i],ext)
    ichip = long(first_el(strsplit(phbase,'_',/extract),/last))
    dum = first_el(strsplit(phbase,'-',/extract),/last)
    refexpnum = first_el(strsplit(dum,'_',/extract))

    ; Load the MCH file, that gives the ordering of the exposures
    mchfile = file_dirname(photfiles[i])+'/'+file_basename(photfiles[i],ext)+'.mch'
    if file_test(mchfile) eq 0 then stop,mchfile,' NOT FOUND'
    LOADMCH,mchfile,mchfilelist
    mchfilelist = file_basename(mchfilelist,'.als')  ; remove .als ending

    chind = where(chstr.chip eq ichip,nchind)
    ; Loop through all of the exposures
    for j=0,nchind-1 do begin
      ; phot files
      if not keyword_set(useast) then begin
        mind = where(tags eq strtrim(chstr[chind[j]].calib_magname,2),nmind)
      ; ast files
      endif else begin
        ; figure out it's order in the MCH file
        mchfilelist_index = where(mchfilelist eq chstr[chind[j]].base,nmchfilelist_index)
        if nmchfilelist_index eq 0 then stop,chstr[chind[j]].base+' NOT FOUND in '+mchfile
        mind = where(tags eq 'MAG'+strtrim(mchfilelist_index[0]+1,2),nmind)
      endelse
      gd = where(phot.(mind[0]) lt 50,ngd)

      ; Load the FITS header
      fitsfile = reduxdir+night+'/'+chstr[0].field+'/'+strtrim(chstr[chind[j]].base,2)+'.fits'
      if file_test(fitsfile) eq 0 then stop,fitsfile,' NOT FOUND'
      head = headfits(fitsfile)

      mjd = PHOTRED_GETMJD(fitsfile,'ctio')

      ; Remove bad data for DECam chip 31
      if (ichip eq 31) and (mjd gt 56660) then begin
        print,'Removing bad data for DECam chip 31'

        ; Determine if als/alf if prob/flag is there
        ; Load ALS files
        if tag_exist(phot,'PROB') then origphotexten='alf' else origphotexten='als'
        origphotfile = reduxdir+night+'/'+chstr[0].field+'/'+strtrim(chstr[chind[j]].base,2)+'.'+origphotexten
        if file_test(origphotfile) eq 0 then stop,origphotfile,' NOT FOUND'

        ; Load original photometry als/alf file
        LOADALS,origphotfile,origphot,origphothead,count=norigphot

        ; ALF: The star IDs are unform across the alf files
        if origphotexten eq 'alf' then begin
          MATCH,phot[gd].id,origphot.id,ind1,ind2,/sort,count=nmatch
          xorig = fltarr(ngd)
          yorig = fltarr(ngd)
          xorig[ind1] = origphot[ind2].x
          yorig[ind1] = origphot[ind2].y
        ; ALS: Need to use the TFR file to get indices for the als file
        endif else begin
          tfrfile = reduxdir+night+'/'+chstr[0].field+'/'+phbase+'.tfr'
          if file_test(tfrfile) eq 0 then stop,tfrfile+' NOT FOUND'
          LOADTFR,tfrfile,tfrlist,tfrstr
          ; What number is it in the frame list
          tfrlist_index = where(tfrlist eq file_basename(origphotfile),ntfrlist_index)          
          if ntfrlist_index eq 0 then stop,file_basename(origphotfile)+' NOT in '+tfrfile
          ; Get indices for the original ALS file, 0-based
          ;   # of stars in TFR file should be same as RAW/AST
          alsindex = reform( tfrstr[gd].index[tfrlist_index]-1 )
          xorig = origphot[alsindex].x
          yorig = origphot[alsindex].y
        endelse

        ; Remove bad measurements
        ; X: 1-1024 okay
        ; X: 1025-2049 bad
        bdind = where(xorig gt 1024,nbdind,comp=gdind,ncomp=ngdind)
        if nbdind gt 0 then begin   ; some bad ones found
          if ngdind eq 0 then begin   ; all bad
            print,'NO useful measurements in ',fitsfile
            undefine,gd
            ngd = 0
          endif else begin
            print,'Removing '+strtrim(nbdind,2)+' bad measurements, '+strtrim(ngdind,2)+' left.'
            REMOVE,bdind,gd
            ngd = n_elements(gd)
          endelse
       endif  ; some bad ones to remove
      endif  ; chip 31



      ; Load the original values from ALS or ALF files
      ;  X, Y, CHI, SHARP, RA, DEC
      if keyword_set(useorig) then begin

      ; FOR NOW DON'T DO THIS!
      ; USING THE SAME X/Y/RA/DEC IS EASIER FOR MATCHING

        ; Determine if als/alf if prob/flag is there
        ; Load ALS files
        if tag_exist(phot,'PROB') then origphotexten='alf' else origphotexten='als'
        origphotfile = reduxdir+night+'/'+chstr[0].field+'/'+strtrim(chstr[chind[j]].base,2)+'.'+origphotexten
        if file_test(origphotfile) eq 0 then stop,origphotfile,' NOT FOUND'

        ; Load original photometry als/alf file
        LOADALS,origphotfile,origphot,origphothead,count=norigphot

        ; MATCH THEM UP

        ; ID, X, Y, MAG, ERR, SKY, ITER, CHI, SHARP

        ; Converting to IDL X/Y convention, starting at (0,0)
        ; DAOPHOT has X/Y start at (1,1)
        x = phot.x - 1.0
        y = phot.y - 1.0

        ; Get RA/DEC coordinates for X/Y
        HEAD_XYAD,head,x,y,ra,dec,/degree

        ; Get ID, MAG, ERR, PROB, FLAG from AST/PHOT file
        stop
      endif

      ; ALLSRC structure schema
      ;  cmbindx   index into ALLOBJ, added later on in smashred_crossmatch.pro
      ;  fix       final ID, added later on in smashred_crossmatch.pro
      ;  cmag/cerr are for calibrated photometry that will be added later
      if ngd gt 0 then begin
        temp = replicate({cmbindx:-1L,chipindx:-1L,fid:'',id:-1L,x:0.0,y:0.0,mag:0.0,err:0.0,$
                          cmag:-1.0,cerr:-1.0,chi:0.0,sharp:0.0,flag:-1,prob:-1.0,ra:0.0d0,dec:0.0d0},ngd)
        ;temp = replicate({id:-1L,x:0.0,y:0.0,mag:0.0,err:0.0,cmag:-1.0,cerr:-1.0,chi:0.0,sharp:0.0,flag:-1,prob:-1.0,ra:0.0d0,dec:0.0d0},ngd)
        struct_assign,phot[gd],temp,/nozero
        temp.mag = phot[gd].(mind[0])
        temp.err = phot[gd].(mind[0]+1)   ; assume error is the next column
        temp.chipindx = chind[j]
      endif

      chstr[chind[j]].refexpnum = refexpnum
      chstr[chind[j]].nsrc = ngd
      chstr[chind[j]].allsrcindx = cur_allsrc_indx

      ; Get astrometric vertices from header
      nx = sxpar(head,'NAXIS1')
      ny = sxpar(head,'NAXIS2')
      head_xyad,head,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/degree
      chstr[chind[j]].vertices_ra = vra
      chstr[chind[j]].vertices_dec = vdec

      ; Add new elements to ALLSRC
      if nallsrc lt cur_allsrc_indx+ngd then begin
        print,'Adding new elements to ALLSRC'
        new = replicate(allsrc_schema,nallsrc+5000000L)  ; add another 5 million elements
        new[0:nallsrc-1] = allsrc
        allsrc = new
        undefine,new
        nallsrc = n_elements(allsrc)  
      endif

      ; Add to ALLSRC structure
      if ngd gt 0 then allsrc[cur_allsrc_indx:cur_allsrc_indx+ngd-1] = temp

      ; Increment index
      cur_allsrc_indx += ngd

      print,j+1,chstr[chind[j]].expnum,chstr[chind[j]].chip,ngd,format='(I4,A12,I5,I10)'
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
Endelse

;stop

end
