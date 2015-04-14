pro imstack,input,noweight=noweight,nozero=nozero,imscale=imscale,clobber=clobber

;+
;
; Combine the deep images for each field/band
;
; INPUTS:
;  input     List of FITS files to combine.
;  /noweight Don't weight the images, i.e. equal weights.
;  /nozero   Don't correct for zero-point offsets.
;  /imscale  Scale the images.  Not set by default.
;  /clobber  Overwritten existing file.
;
; OUTPUTS:
;  Combined FITS file.
;
; USAGE:
;  IDL>imstack,input
;
; By D. Nidever  Nov. 2014

forward_function wcsfit_find
resolve_routine,'wcsfit',/compile_full_file

t0 = systime(1)

if n_elements(input) eq 0 then begin
  print,'Syntax - imstack,input'
  return
endif

; Default parameters
if n_elements(noweight) eq 0 then noweight=0
if n_elements(nozero) eq 0 then nozero=0
if n_elements(imscale) eq 0 then imscale=0


; Loading the input
LOADINPUT,input,infiles
nfiles = n_elements(infiles)
if nfiles eq 0 then begin
  print,'No files to process'
  return
endif

cd,current=curdir

; Check that we have the necessary files
; daofindphot.sh, photo.opt
if file_test('daofindphot.sh') eq 0 then begin
  print,'Script daofindphot.sh NOT FOUND'
  return
endif
if file_test('photo.opt') eq 0 then begin
  print,'Script photo.opt NOT FOUND'
  return
endif


nfiles = n_elements(infiles)
print,'Combining ',strtrim(nfiles,2),' files: ',infiles


; Check if the output file exists already
combfile = file_basename(infiles[0],'.fits')+'_comb.fits'
if file_test(combfile) eq 1 and not keyword_set(clobber) then begin
  print,combfile,' exists ALREADY and /clobber not set'
  return
endif

; Get source lists with coordinates
undefine,allcat
for l=0,nfiles-1 do begin

   tempbase = file_basename(infiles[l],'.fits')+'.temp'
   tempfile = tempbase+'.fits'
   file_copy,infiles[l],tempfile,/allow,/over

   ; Run daofindphot.sh
   MKOPT,tempfile
   print,'Running DAOPHOT find/phot on ',tempfile
   file_delete,tempbase+['.coo','.ap','.log'],/allow  ; delete output files
   spawn,['./daofindphot.sh',file_basename(tempfile,'.fits')],out,outerr,/noshell
   file_delete,tempfile,/allow

   LOADCOO,tempbase+'.coo',cat,cathead
   LOADAPER,tempbase+'.ap',phot,phothead

   ;WCSFIT_FIND,infiles[l],cat
   ;gdcat = where(cat.mag gt 0,ngdcat)
   naper = n_elements(phot[0].mag)
   useaper = (naper-3) > 0
   gdcat = where(phot.mag[useaper] lt 50 and phot.mag[useaper] gt 0 and $
                 cat.sharp gt 0.2 and cat.sharp lt 1.0 and abs(cat.round) lt 1 and abs(cat.round2) lt 1,ngdcat)
   print,strtrim(ngdcat,2),' sources with decent photometry'
   phot1 = phot[gdcat]
   cat1 = cat[gdcat]
   newcat = replicate({id:0L,x:0.0d0,y:0.0d0,mag:0.0,err:0.0,sky:0.0,niter:0L,chi:0.0,sharp:0.0},n_elements(phot1))
   struct_assign,phot1,newcat
   newcat.mag = phot1.mag[useaper]
   newcat.err = phot1.err[useaper]
   newcat.niter = 1
   newcat.chi = 1.0
   newcat.sharp = cat1.sharp
   ; construct the header
   ;fhead = headfits(infiles[l])
   ;nx = sxpar(fhead,'NAXIS1')
   ;ny = sxpar(fhead,'NAXIS2')
   ;rdnoise = sxpar(fhead,'RDNOISE')
   ;gain = sxpar(fhead,'GAIN')
   ;saturate = sxpar(fhead,'SATURATE',count=nsaturate)
   ;if nsaturate eq 0 then saturate=50000.
   ;head=' NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD'
   ;head=[head,'  1  '+strtrim(nx,2)+'  '+strtrim(ny,2)+'   0.0  '+string(saturate,format='(F8.1)')+'    5.00    0.00    '+$
   ;      stringize(gain,ndec=2)+'    '+stringize(rdnoise/gain,ndec=2)+'    5.00']
   alsfile = file_basename(infiles[l],'.fits')+'.temp.als'
   writeals,alsfile,newcat,cathead   ; head
   add_tag,cat,'file',infiles[l],cat
   push,allcat,cat
endfor
alsfiles = file_basename(infiles,'.fits')+'.temp.als'
DAOMATCH,alsfiles,/verbose

mchbase = file_basename(infiles[0],'.fits')
mchfile = mchbase+'.temp.mch'


;###########################################
; STEP 1: IMALIGN PREP

print,'STEP 1: GETTING WEIGHTS'

;-----------------------------------
; Computs Weights
ALLFRAME_GETWEIGHTS,mchfile,weights,scales,sky,raw2=raw2
invscales = 1.0/scales
bdscale = where(scales lt 1e-5 or invscales gt 900,nbdscale)
if nbdscale gt 0 then begin
  scales[bdscale] = 1.0
  invscales[bdscale] = 1.0
  weights[bdscale] = 0.0
endif
weightfile = mchbase+'.weights'
WRITECOL,weightfile,weights,fmt='(F10.6)'
scalefile = mchbase+'.scale'
WRITECOL,scalefile,invscales,fmt='(F10.5)'  ; want to scale it UP
zerofile = mchbase+'.zero'
WRITECOL,zerofile,-sky,fmt='(F10.2)'  ; want to remove the background, set to 1st frame


;---------------------------------------
; Get X/Y translations using DAOMASTER
;  NO ROTATION ONLY SHIFTS
;  Need these shifts for IMSHIFT
print,'Measuring X/Y shifts'
shiftmch = mchbase+'_shift'
FILE_COPY,mchbase+'.temp.mch',shiftmch+'.mch',/overwrite,/allow
;FILE_COPY,mchbase+'.mch',shiftmch+'.mch',/overwrite,/allow
; Make the DAOMASTER script
undefine,cmdlines
PUSH,cmdlines,'#!/bin/csh'
PUSH,cmdlines,'set input=${1}'
PUSH,cmdlines,'daomaster << DONE'
PUSH,cmdlines,'${input}.mch'
PUSH,cmdlines,'1,1,1'
PUSH,cmdlines,'99.'
PUSH,cmdlines,'2'
PUSH,cmdlines,'10'
;for i=0,nfiles-1 do PUSH,cmdlines,''  ; for MONGO daomaster version
PUSH,cmdlines,'5'
PUSH,cmdlines,'4'
PUSH,cmdlines,'3'
PUSH,cmdlines,'2'
PUSH,cmdlines,'1'
PUSH,cmdlines,'0'
PUSH,cmdlines,'n'
PUSH,cmdlines,'n'
PUSH,cmdlines,'n'
PUSH,cmdlines,'n'
PUSH,cmdlines,'y'
PUSH,cmdlines,''
PUSH,cmdlines,''
PUSH,cmdlines,'n'
PUSH,cmdlines,'n'
PUSH,cmdlines,'n'
PUSH,cmdlines,'DONE'
tempscript = MKTEMP('daomaster')   ; absolute filename
WRITELINE,tempscript,cmdlines
FILE_CHMOD,tempscript,'755'o
; Run DAOMASTER                                                                                                                           

              
cmd2 = tempscript+' '+shiftmch
SPAWN,cmd2,out2,errout2
; Remove temporary DAOMASTER script
FILE_DELETE,tempscript,/allow_non
LOADMCH,shiftmch+'.mch',files2,trans2

xshift = reform(trans2[*,0])
yshift = reform(trans2[*,1])
xyshifts = [[xshift],[yshift]]
print,'Image shifts'
for i=0,nfiles-1 do print,infiles[i],xshift[i],yshift[i]


;-----------------------------------
; Create imalign prep files
; This is done by preimalign_k.sh
; Need an input list of fits files
; Need an output list of fits files
; Shift file
;base = file_basename(files,'.als')
base = file_basename(infiles,'.fits')
fitsfiles = base+'.fits'
outfiles = base+'.shft.fits'
infile = mchbase+'.inlist'
outfile = mchbase+'.outlist'
;WRITELINE,infile,fitsfiles   ; this is done below now with the temp
;files
WRITELINE,outfile,outfiles
; Remove outfiles
FILE_DELETE,outfiles,/allow
; shift list
shiftfile = mchbase+'.shift'
WRITECOL,shiftfile,xshift,yshift,fmt='(2F15.4)'


; Make temporary files for bad pixel fixing and combining
;  FIXIMAGES doesn't work properly on the shifted images
;  because the interpolation can bring the bad pixel values down
;  below the saturation threshold, and we don't want to touch
;  the original images.
tempfits = base+'.temp.fits'
FILE_COPY,fitsfiles,tempfits,/overwrite,/allow
WRITELINE,infile,tempfits



;###########################################
; STEP 2: FIX BAD PIXELS
;printlog,logf,'---------------------------'
print,'STEP 2: FIXING BAD PIXELS'
;printlog,logf,'---------------------------'
FIXIMAGES,'@'+infile,satlevel=satlevel  ;6e4
; This also makes the FILE.mask.fits files for each image

; Find the maximum saturation level
satlevelarr = fltarr(nfiles)
for i=0,nfiles-1 do begin
  ;head = headfits(base[i]+'.fits')
  FITS_READ,base[i]+'.fits',im,head,/no_abort
  saturate = sxpar(head,'SATURATE',count=nsaturate)
  if nsaturate eq 0 then saturate=max(im)-1000.
  satlevelarr[i] = saturate
endfor
maxsatlevel = max(satlevelarr)


;###########################################
; STEP 3: IMALIGN
;  This figures out the X/Y-shifts between the images
;  and creates the shifted images (".shft.fits")
;printlog,logf,'-----------------'
print,'STEP 3: IMALIGN'
;printlog,logf,'-----------------'
reffile = mchbase+'.fits'

; IMALIGN basically is a script that runs:
;  IMCENTROID - to compute the shifts
;  IMSHIFT - shifts the images
;  IMCOPY - trims the images


; First, shift the images
print,'Shifting the images'
IRAF_IMSHIFT,'@'+infile,'@'+outfile,shifts_file=shiftfile,interp_type='linear',$
             boundary_type='constant',constant=0,irafdir=irafdir,error=imshifterror,/verbose
if n_elements(imshifterror) ne 0 then begin
  print,'ERROR in IRAF_IMSHIFT'
  print,imshifterror
  error = imshifterror
  return
endif

; Check that the output files are there
if total(file_test(outfiles)) ne n_elements(outfiles) then begin
  bd = where(file_test(outfiles) eq 0,nb)
  print,outfiles[bd],' NOT FOUND.'
  stop
  return
endif

; Trim the images
if keyword_set(trimcomb) then begin

  ; Calculate the trim section
  hd = headfits(reffile)
  xsize = lonarr(nfiles)+sxpar(hd,'NAXIS1')
  ysize = lonarr(nfiles)+sxpar(hd,'NAXIS2')
  IA_TRIM,xshift,yshift,xsize,ysize,trimsection
  xoff = trimsection[0]-1
  yoff = trimsection[2]-1

  ; Trim the shifted images
  print,'Trimming the shifted images'
  xstart = trimsection[0]-1
  xstop = trimsection[1]-1
  ystart = trimsection[2]-1
  ystop = trimsection[3]-1

  for i=0,nfiles-1 do begin
    FITS_READ,outfiles[i],im,head
    newim = im[xstart:xstop,ystart:ystop]
    MWRFITS,newim,outfiles[i],head,/create,/silent
 end
  ; could also use IRAF_IMCOPY here instead


; Don't trim the images
endif else begin
  hd = headfits(reffile)
  xsize = sxpar(hd,'NAXIS1')
  ysize = sxpar(hd,'NAXIS2')
  trimsection = [1,xsize,1,ysize]
  xoff = 0
  yoff = 0
endelse

; Delete the temporary FITS files
FILE_DELETE,tempfits,/allow



;###########################################
; STEP 4: MAKE BAD PIXEL/WEIGHT MAP
; in how many images does the pixel need to be bad??
; 1. shift the masks (created by FIXIMAGES.PRO)
;     using the shifts from IRAF_IMALIGN
; 2. trim the masks
; 3. combine the masks
;printlog,logf,'-------------------------------'
print,'STEP 4: Making BAD PIXEL MASK'
;printlog,logf,'-------------------------------'

; Make lists
maskfiles = FILE_DIRNAME(tempfits)+'/'+FILE_BASENAME(tempfits,'.fits')+'.mask.fits'
outmaskfiles = FILE_DIRNAME(tempfits)+'/'+FILE_BASENAME(tempfits,'.fits')+'.mask.shft.fits'
maskinfile = mchbase+'.maskinlist'
maskoutfile = mchbase+'.maskoutlist'
maskshiftsfile = mchbase+'.maskshifts'
WRITELINE,maskinfile,maskfiles
WRITELINE,maskoutfile,outmaskfiles
FILE_DELETE,outmaskfiles,/allow
strxyshifts = strarr(nfiles)
for i=0,nfiles-1 do strxyshifts[i] = strjoin(reform(xyshifts[i,*]),'  ')
WRITELINE,maskshiftsfile,strxyshifts

; Run IMSHIFT
;  set boundary to 0=BAD
print,'Shifting masks'
undefine,iraflines
push,iraflines,'cd '+curdir
push,iraflines,'images'
push,iraflines,'imgeom'
push,iraflines,'imshift("@'+maskinfile+'","@'+maskoutfile+'",shifts_file="'+maskshiftsfile+'",'+$
               'interp_type="linear",boundary_typ="constant",constant=0)'
push,iraflines,'logout'
imshiftscript = curdir+'/'+mchbase+'.imshift'
WRITELINE,imshiftscript,iraflines
IRAF_RUN,imshiftscript,irafdir,out=out,/silent,error=iraferror
if n_elements(iraferror) ne 0 then begin
  print,'ERROR in running IMSHIFT with IRAF_RUN'
  print,iraferror
  error = iraferror
  return
endif

; Trim
if keyword_set(trimcomb) then begin
  print,'Trimming masks'
  xstart = trimsection[0]-1     ; should be same as xoff
  xstop = trimsection[1]-1
  ystart = trimsection[2]-1     ; should be same as yoff
  ystop = trimsection[3]-1

  for i=0,nfiles-1 do begin
    FITS_READ,outmaskfiles[i],im,head
    sz = size(im)
    newim = im[xstart:xstop,ystart:ystop]
    ; Add LTV1/LTV2 to the header
    ;  these are IRAF keywords to convert from logical to physical coords
    ltv1 = sxpar(head,'LTV1')  ; 0 if not found
    ltv2 = sxpar(head,'LTV2')  ; 0 if not found
    sxaddpar,head,'LTV1',ltv1-xstart
    sxaddpar,head,'LTV2',ltv2-ystart
    MWRFITS,newim,outmaskfiles[i],head,/create,/silent
 endfor
endif

; Combining masks
print,'Combining masks'
undefine,bpmsum
for i=0,nfiles-1 do begin
  FITS_READ,outmaskfiles[i],im,head
  if i eq 0 then begin
    bpmsum = im
    whead = head
 endif else begin
    bpmsum += im
 endelse
end


;; masks have 0-bad, 1-good.
;; anything with less than 1.0 is considered bad
;; weight map, -1 is bad, +1 is good
;weightmap = -2.0*float(bpmsum lt nfiles) + 1.
;combweightfile = mchbase+'_comb.mask.fits'
;FITS_WRITE,combweightfile,weightmap,whead
;
; THIS IS NOW DONE BELOW AFTER THE IMAGE IS COMBINED
; DEPENDING ON IF THE IMAGES ARE SCALED OR NOT!!!

;stop

;###########################################
; STEP 5: COMBINE IMAGES
;print,'-------------------'
print,'STEP 5: IMCOMBINE'
;print,'-------------------'



; SCALE the images for combining
;-------------------------------
if keyword_set(imscale) then begin

  ; Put BPM mask names in file headers
  ;  these will be used by IMCOMBINE
  for i=0,nfiles-1 do begin
    head = headfits(outfiles[i])
    sxaddpar,head,'BPM',outmaskfiles[i]
    modfits,outfiles[i],0,head
  endfor

  ; Combine the frames WITH scaling/offset/masking, for the bright stars
  ;printlog,logf,'Creating SCALED image'
  combfile = mchbase+'_comb.fits'
  FILE_DELETE,combfile,/allow
  FILE_DELETE,mchbase+'_comb.bpm.pl',/allow
  IRAF_IMCOMBINE,'@'+outfile,combfile,combine='average',reject='avsigclip',$
                 weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',$
                 irafdir=irafdir,error=imcombineerror2,scale='@'+scalefile,zero='@'+zerofile,$
                 masktype='badvalue',maskvalue=0,bpmasks=mchbase+'_comb.bpm'

  if n_elements(imcombineerror2) ne 0 then begin
    printlog,logf,'ERROR in IRAF_IMCOMBINE'
    printlog,logf,imcombineerror2
    error = imcombineerror2
    return
  endif

  ; Convert BPM mask from PL to FITS
  FILE_DELETE,mchbase+'_comb.bpm.fits',/allow
  undefine,lines
  cd,current=curdir
  push,lines,'cd '+curdir
  push,lines,'imcopy '+mchbase+'_comb.bpm.pl '+mchbase+'_comb.bpm.fits'
  push,lines,'logout'
  tempfile = mktemp('tiraf')
  WRITELINE,tempfile,lines
  IRAF_RUN,tempfile,irafdir,silent=silent,out=out,error=error

  ; Delete temporary scripts and PL file
  FILE_DELETE,[tempfile,mchbase+'_comb.bpm.pl'],/allow


  ; Fix the rdnoise and background/sky level and saturate
  ;  the bad pixels for DAOPHOT
  ;------------------------------------------------------

  ; 10/02/12
  ; THE IMAGES ARE (1) ZERO-SUBTRACTED, (2) SCALED, AND (3) WEIGHT AVERAGED
  ; The algorithm is:
  ; 1.) add zero-level correction.  im = im+zero
  ; 2.) scale the images.  im = im*scale
  ; 3.) take weighted average.  combim=total(weight*im)
  ;      there is also clipping that takes place during the averaging
  ; The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2))
  ; A gain that changes from frame to frame could be problematic,
  ; but this shouldn't happen since it's the same chip from the same night.

  ; IMCOMBINE wants rdnoise in electrons and gain in electrons/DN.
  ; DAOPHOT expects rdnoise in DN.  That's why mkopt converts
  ;  it with the gain.  So we are fine.  The header should have
  ;  rdnoise in ELECTRONS.

  ; page 63-65 of daophot2.pdf shows how you need to modify rdnoise/gain
  ; when averaging/summing frames. in observing/mosaic/.

  ; Load the IMCOMBINE output combined file and BPM
  FITS_READ,combfile,combim,combhead
  FITS_READ,mchbase+'_comb.bpm.fits',badmask,maskhead  ; 0-good, 1-bad


  ; Fix the gain
  ; For N averaged frames gain(N)=N*gain(1)
  ; Leave the gain as is!  We are scaling everything to the reference
  ; and using its gain.  It's nearly impossible to figure out the real
  ; gain since we are scaling the images and then taking a weighted
  ; average with outlier rejection.  Find a gain that properly
  ; describes/follows Poisson errors for the final combined image is
  ; difficult/impossible.  But that's okay.  This is just for source
  ; detection and DAOPHOT FIND just cares about the noise in the
  ; background.  We just need to ensure that the sky and rdnoise
  ; are correct.

  ; Fix the rdnoise
  ; The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2))
  rdnoisearr = fltarr(nfiles)
  for i=0,nfiles-1 do rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits')
  ;  the "scales" array here is actually 1/scales used by IMCOMBINE.
  rdnoise = sqrt(total((weights*rdnoisearr/scales)^2))
  dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey) ; get keyword
  sxaddpar,combhead,rdnoisekey,rdnoise

  ; Fix the sky
  ; DAOPHOT FIND computes the random error per pixel in ADU as
  ; noise = sqrt( sky level/gain + rdnoise^2)
  ; So it assumes that the noise in the background is sqrt(sky/gain)
  ; in ADU.  We need to set the sky level so this is correct.
  ; The final noise should be 
  ; final sky noise = sqrt(total((weights*scale*sqrt(sky/gain))^2)) 
  ; So the final sky level should be
  ; final sky = total((weights*scale*sqrt(sky/gain))^2)*gain
  gain = PHOTRED_GETGAIN(combfile,keyword=gainkey)
  comb_sky = total((weights*sqrt((sky>0)/gain)/scales)^2)*gain
  ; the "scales" array here is actually 1/scale
  combim += float(comb_sky)  ; keep it float


  ; set the maximum to a "reasonable" level
  ; Rescale the image and increase the gain
  if max(combim) gt 50000 then begin
    rescale = 50000./max(combim)
    combim = combim*rescale
    sxaddpar,combhead,gainkey,gain/rescale
    ; rdnoise does NOT get modified since it's in electrons
    ; we just need to modify the gain which takes you from ADUs to electrons
  endif

  maskdatalevel = max(combim) + 10000       ; set "bad" data level above the highest "good" value
  combim2 = combim*(1-badmask) + maskdatalevel*badmask    ; set bad pixels to maskdatalevel
  MWRFITS,combim2,combfile,combhead,/create  ; fits_write can create an empty PDU

  ;; Create the weight map for Sextractor using the BPM output by IMCOMBINE
  ;;  bad only if bad in ALL images
  ;weightmap = -2.0*float(badmask eq 1) + 1.0
  ;combweightfile = mchbase+'_comb.mask.fits'
  ;MWRFITS,weightmap,combweightfile,whead,/create

  
; NO SCALING of the images for combining
;---------------------------------------
Endif else begin


  ; Put BPM mask names in file headers
  ;  these will be used by IMCOMBINE
  for i=0,nfiles-1 do begin
    head = headfits(outfiles[i])
    sxaddpar,head,'BPM',outmaskfiles[i]
    modfits,outfiles[i],0,head
  endfor

  ; no scaling of the images for combining

  combfile = mchbase+'_comb.fits'
  FILE_DELETE,combfile,/allow
  FILE_DELETE,mchbase+'_comb.bpm.pl',/allow
  if keyword_set(noweight) then weightinput='none' else weightinput='@'+weightfile
  if keyword_set(nozero) then zeroinput='none' else zeroinput='@'+zerofile
  IRAF_IMCOMBINE,'@'+outfile,combfile,combine='average',reject='avsigclip',$
                 weight=weightinput,scale='none',zero=zeroinput,$,rdnoise='!rdnoise',gain='!gain',$
                 ;weight='none',scale='none',zero='none',rdnoise='!rdnoise',gain='!gain',$
                 irafdir=irafdir,error=imcombineerror2,masktype='badvalue',maskvalue=0,bpmasks=mchbase+'_comb.bpm'
  ; "sum" doesn't work well, pixels bad in any image are bad in
  ; the sum image.  "average" works better, less "bad" pixels in final image.

  if n_elements(imcombineerror2) ne 0 then begin
    print,'ERROR in IRAF_IMCOMBINE'
    print,imcombineerror2
    error = imcombineerror2
    return
  endif


  ; Convert BPM mask from PL to FITS
  FILE_DELETE,mchbase+'_comb.bpm.fits',/allow
  undefine,lines
  cd,current=curdir
  push,lines,'cd '+curdir
  push,lines,'imcopy '+mchbase+'_comb.bpm.pl '+mchbase+'_comb.bpm.fits'
  push,lines,'logout'
  tempfile = mktemp('tiraf')
  WRITELINE,tempfile,lines
  IRAF_RUN,tempfile,irafdir,silent=silent,out=out,error=error

  ; Delete temporary scripts and PL file
  FILE_DELETE,[tempfile,mchbase+'_comb.bpm.pl'],/allow

  ; Fix the rdnoise and saturate
  ;  the bad pixels for DAOPHOT
  ;------------------------------------------------------


  ; Load the IMCOMBINE output combined file and BPM
  FITS_READ,combfile,combim,combhead,/no_abort
  FITS_READ,mchbase+'_comb.bpm.fits',badmask,maskhead,/no_abort  ; 0-good, 1-bad

  ;weights = fltarr(nfiles)+1.0
  ; do we want to weight by the exptime??

  ; Convert from average to sum
  combim *= nfiles

  ; Fix the rdnoise
  ; The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise)^2))
  rdnoisearr = fltarr(nfiles)
  for i=0,nfiles-1 do rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits')
  rdnoise = sqrt(total((weights*rdnoisearr)^2))
  dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey) ; get keyword
  sxaddpar,combhead,rdnoisekey,rdnoise

  ; Increase the Exposure Time
  exptimearr = fltarr(nfiles)
  for i=0,nfiles-1 do exptimearr[i] = PHOTRED_GETEXPTIME(base[i]+'.fits')
  totexptime = total(exptimearr)
  ;dummy = PHOTRED_GETEXPTIME(combfile,keyword=exptimekey) ; get keyword
  sxaddpar,combhead,'EXPTIME',totexptime

  ; Fix the sky
  ; So the final sky level should be
  ; final sky = total((weights*scale*sqrt(sky/gain))^2)*gain
  gain = PHOTRED_GETGAIN(combfile,keyword=gainkey)
  comb_sky = total((weights*sqrt(sky/gain))^2)*gain
  ; the "scales" array here is actually 1/scale
  combim += comb_sky


  ;; set the maximum to a "reasonable" level
  ;; Rescale the image and increase the gain
  ;if max(combim) gt 50000 then begin
  ;  rescale = 50000./max(combim)
  ;  combim = combim*rescale
  ;  sxaddpar,combhead,gainkey,gain/rescale
  ;  ; rdnoise does NOT get modified since it's in electrons
  ;  ; we just need to modify the gain which takes you from ADUs to electrons
  ;endif


  ;; Making Sextractor "weight" map file
  ;;------------------------------------
  ;; masks have 0-bad, 1-good.
  ;; anything with less than 1.0 is considered bad
  ;; weight map, -1 is bad, +1 is good
  ;; "bpmsum" is the SUM of the bad pixel masks
  ;; consider a pixel bad that is bad in ANY image
  ;weightmap = -2.0*float(bpm lt nfiles) + 1.
  ;combweightfile = mchbase+'_comb.mask.fits'
  ;FITS_WRITE,combweightfile,weightmap,whead

  ;---------------------------------------------
  ; SATURATE BAD pixels in the COMBINED IMAGE
  ; DAOPHOT needs to have the bad pixels "saturated",
  ; SExtractor will know which pixels are bad from the "weight" map.
  ;
  ; We could skip the fiximage.pro step but we still need the
  ; individual bpm masks and setting the bad pixels to the background
  ; probably helps in the IMALIGN/IMCOMBINE steps.
  ;print,''
  print,'"Saturating" bad pixels in the COMBINED image'
  ;print,''

  ; Use BOTH the bad mask from IMCOMBINE and the sum of the individual masks
  ; bpmsum (SUM of bpm masks) - 0-bad, 1-good
  ; badmask (IMCOMBINE) - 0-good, 1-bad
  fbadmask = float(bpmsum eq 0 or badmask eq 1)  ; bad in all or bad from IMCOMBINE (rejected)
  ;badmask = float(weightmap lt 0.5)
  maskdatalevel = max(combim) + 10000       ; set "bad" data level above the highest "good" value
  combim2 = combim*(1.0-fbadmask) + maskdatalevel*fbadmask    ; set bad pixels to 100,000

  ; Set SATURATE in the header
  sxaddpar,combhead,'SATURATE',maskdatalevel

  ; Add info to the header
  for i=0,nfiles-1 do sxaddhist,infiles[i]+' Xsh='+stringize(xshift[i],ndec=2)+' Ysh='+stringize(yshift[i],ndec=2),combhead

  ; Write final file
  ;FITS_WRITE,combfile,combim2,combhead
  MWRFITS,combim2,combfile,combhead,/create

  ; Set airmass to midpoint of exposure
  am = photred_getairmass(combfile,obs='CTIO',/update,/recalculate)

Endelse


; DELETE TEMPORARY FILES
;-----------------------

; Delete the temporary DAOPHOT files
FILE_DELETE,file_basename(infiles,'.fits')+'.temp.opt',/allow
FILE_DELETE,file_basename(infiles,'.fits')+'.temp.als.opt',/allow
FILE_DELETE,file_basename(infiles,'.fits')+'.temp.coo',/allow
FILE_DELETE,file_basename(infiles,'.fits')+'.temp.ap',/allow
FILE_DELETE,file_basename(infiles,'.fits')+'.temp.log',/allow

; Delete the shifted images
READLINE,outfile,shiftedfiles
FILE_DELETE,shiftedfiles,/allow,/quiet

; Delete mask files
FILE_DELETE,[maskfiles,outmaskfiles],/allow
FILE_DELETE,[maskinfile,maskoutfile,maskshiftsfile,imshiftscript],/allow
FILE_DELETE,[alsfiles,infile,outfile,shiftfile,shiftmch+'.mch'],/allow
FILE_DELETE,file_basename(alsfiles[0],'.als')+['.mch','.tfr','.raw'],/allow

; Final files should be: BASE_comb.fits, BASE_comb.bpm.fits, BASE_comb.mask.fits

dt = systime(1)-t0
print,'dt=',strtrim(dt,2),' sec.'

BOMB:

;stop

end
