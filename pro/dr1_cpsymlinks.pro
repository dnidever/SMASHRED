pro dr1_cpsymlinks

; Figure out the CP files we want to release and sym links to produce
; on zeus1

catdir = '/data/smash/cp/red/photred/'
;indir = '/data/smash/cp/red/photred/'
;outdir = '/data/smash/dr1/'
;indir = '/datalab/users/dnidever/smash/cp/red/photred/'
;outdir = '/datalab/users/dnidever/smash/dr1/'
outdir = '/data/smash/dr1/'

; Restore the DR1 fields
dr1 = mrdfits('~/projects/SMASHRED/data/smash_dr1_fields.fits',1)
dr1.field = strtrim(dr1.field,2)
ndr1 = n_elements(dr1)

; Get smash reduction info
smashred_getredinfo,allredinfo,/silent
bd = where(strtrim(allredinfo.field,2) eq '',nbd)
if nbd gt 0 then remove,bd,allredinfo

; Load the full file paths
;;fpath = mrdfits('alldec2.fits',1)
;fpath = mrdfits(catdir+'observations/mss_filenames.fits.gz',1)
;;for i=0,n_tags(fpath)-1 do fpath.(i)=strtrim(fpath.(i),2)
;fpath.base = strtrim(fpath.base,2)

; All RAW files
allrawstr = mrdfits(catdir+'observations/decam_all_raw.vot.fits',1)
add_tag,allrawstr,'expnum','',allrawstr
allrawstr.expnum = strmid(file_basename(allrawstr.dtacqnam,'.fits.fz'),6)

; Read in the VO tables
;read_votable,catdir+'observations/smash_raw_vot.xml',rawstr
;add_tag,rawstr,'expnum','',rawstr
;rawstr.expnum = strmid(file_basename(rawstr.dtacqnam,'.fits.fz'),6)
;add_tag,rawstr,'smashfield','',rawstr
;smashred_getfieldname,rawstr,rawinfo
;rawstr.smashfield = rawinfo.field
;add_tag,rawstr,'nsaname','',rawstr
;rawstr.nsaname = strtrim(strmid(file_basename(rawstr.reference),9),2)
;;mwrfits,rawstr,catdir+'observations/smash_raw_vot.fits',/create
rawstr = mrdfits(catdir+'observations/smash_raw_vot.fits',1)

;read_votable,catdir+'observations/smash_instcal_vot.xml',instcalstr
;add_tag,instcalstr,'expnum','',instcalstr
;instcalstr.expnum = strmid(file_basename(instcalstr.dtacqnam,'.fits.fz'),6)
;add_tag,instcalstr,'smashfield','',instcalstr
;smashred_getfieldname,instcalstr,instcalinfo
;instcalstr.smashfield = instcalinfo.field
;add_tag,instcalstr,'nsaname','',instcalstr
;instcalstr.nsaname = strtrim(strmid(file_basename(instcalstr.reference),9),2)
;add_tag,instcalstr,'fullpath','',instcalstr
;match,fpath.base,instcalstr.nsaname,ind1,ind2,/sort
;instcalstr[ind2].fullpath = strtrim(fpath[ind1].filename,2)
;;mwrfits,instcalstr,catdir+'observations/smash_instcal_vot.fits',/create
instcalstr = mrdfits(catdir+'observations/smash_instcal_vot.fits',1)

;read_votable,catdir+'observations/smash_resampled_vot.xml',resampstr
;add_tag,resampstr,'expnum','',resampstr
;resampstr.expnum = strmid(file_basename(resampstr.dtacqnam,'.fits.fz'),6)
;add_tag,resampstr,'smashfield','',resampstr
;smashred_getfieldname,resampstr,resampinfo
;resampstr.smashfield = resampinfo.field
;add_tag,resampstr,'nsaname','',resampstr
;resampstr.nsaname = strtrim(strmid(file_basename(resampstr.reference),9),2)
;add_tag,resampstr,'fullpath','',resampstr
;match,fpath.base,resampstr.nsaname,ind1,ind2,/sort
;resampstr[ind2].fullpath = strtrim(fpath[ind1].filename,2)
;;mwrfits,resampstr,catdir+'observations/smash_resampled_vot.fits',/create
resampstr = mrdfits(catdir+'observations/smash_resampled_vot.fits',1)

;read_votable,catdir+'observations/smash_stacked_vot.xml',stackstr
;add_tag,stackstr,'smashfield','',stackstr
;smashred_getfieldname,stackstr,stackinfo
;stackstr.smashfield = stackinfo.field
;add_tag,stackstr,'nsaname','',stackstr
;stackstr.nsaname = strtrim(strmid(file_basename(stackstr.reference),9),2)
;add_tag,stackstr,'fullpath','',stackstr
;match,fpath.base,stackstr.nsaname,ind1,ind2,/sort
;stackstr[ind2].fullpath = strtrim(fpath[ind1].filename,2)
;;mwrfits,stackstr,catdir+'observations/smash_stacked_vot.fits',/create
stackstr = mrdfits(catdir+'observations/smash_stacked_vot.fits',1)
stackstr.fullpath = strtrim(stackstr.fullpath,2)

; Load the header information from the mass store stack files
stackhdstr = mrdfits(catdir+'observations/smash_stacked_vot_scrape.fits',1)
stackhdstr.filename = strtrim(stackhdstr.filename,2)
add_tag,stackhdstr,'minexptime',0.0,stackhdstr
add_tag,stackhdstr,'maxexptime',0.0,stackhdstr
add_tag,stackhdstr,'totexptime',0.0,stackhdstr
add_tag,stackhdstr,'exptimebin','',stackhdstr

; Get the minimum and maximum exptime per stack
for i=0,n_elements(stackhdstr)-1 do begin
  names = [stackhdstr[i].imcmb001, stackhdstr[i].imcmb002, stackhdstr[i].imcmb003, stackhdstr[i].imcmb004,$
           stackhdstr[i].imcmb005, stackhdstr[i].imcmb006, stackhdstr[i].imcmb007, stackhdstr[i].imcmb008,$
           stackhdstr[i].imcmb009, stackhdstr[i].imcmb010]
  gdnames = where(names ne 'None',ngdnames)
  names = names[gdnames]
  expnum = strmid(names,6)
  exptime = fltarr(ngdnames)-1
  filter = strarr(ngdnames)
  match,allrawstr.expnum,expnum,ind1,ind2,/sort,count=nmatch
  exptime[ind2] = allrawstr[ind1].exposure
  filter[ind2] = allrawstr[ind1].filter
  if nmatch ne ngdnames then stop,'not all matched'
  ;for j=0,ngdnames-1 do begin
  ;  match,allrawstr.expnum,expnum[j],ind1,ind2,/sort,count=nmatch
  ;  if nmatch gt 0 then exptime[j]=allrawstr[ind1[0]].exposure
  ;endfor
  ; Some stacks are missing filter information from the NSA VO table
  if strtrim(stackhdstr[i].filter,2) eq 'None' then stackhdstr[i].filter=filter[0]
  gdexptime = where(exptime gt 0,ngdexptime)
  stackhdstr[i].minexptime = min(exptime[gdexptime])
  stackhdstr[i].maxexptime = max(exptime[gdexptime])
  stackhdstr[i].totexptime = total(exptime[gdexptime])
  
  ; Decide on exptime bin names
  ;stop

endfor

; Each group has five files which go together
;  use almost all of the PLFNAME for the group name
;DEC13A_20130317_8475ffb-g-20130318T002924_sd   dqmask
;DEC13A_20130317_8475ffb-g-20130318T002924_si   image
;DEC13A_20130317_8475ffb-g-20130318T002924_sw   wtmap
;DEC13A_20130317_8475ffb-g-20130318T002924_si1   image1
;DEC13A_20130317_8475ffb-g-20130318T002924_se   expmap

; Create stack group IDs that are unique and can be used to associate
;  the various files together (image, dqmask, wtmap, expmap)
add_tag,stackhdstr,'groupid','',stackhdstr
;Make stack group id, maybe something like: Field Night filter skysub expbin procdir plver 
;groupid = stackhdstr.dtcaldat+'.'+strmid(stackhdstr.filter,0,1)+'.'+stackhdstr.prodtype
;groupid = stackhdstr.night
groupid = strmid(stackhdstr.plfname,0,41)
stackhdstr.groupid = groupid
; Create stack id that can be duplicated
;Make id for stack that could be duplicated, maybe something like: Field Night filter skysub expbin 
add_tag,stackhdstr,'stackid','',stackhdstr
stackid = stackhdstr.dtcaldat+'.'+strmid(stackhdstr.filter,0,1)
stackhdstr.stackid = stackid

; Fix PLVER in stackhdstr
bd = where(strtrim(stackhdstr.plver,2) eq 'None',nbd)
if nbd gt 0 then stackhdstr[bd].plver='V1.0.0'

; All of the unique SMASH file names
;nsaname = [rawstr.nsaname,instcalstr.nsaname,resampstr.nsaname,stackstr.nsaname]

; We need more information for these files
;gstackim = where(strtrim(stackstr.prodtype,2) eq 'image' or strtrim(stackstr.prodtype,2) eq 'image1')
;stackimstr = stackstr[gstackim]
;mwrfits,stackimstr,'smash_stack_images.fits',/create

; Make main output directories
if file_test(outdir+'raw',/directory) eq 0 then file_mkdir,outdir+'raw'
if file_test(outdir+'instcal',/directory) eq 0 then file_mkdir,outdir+'instcal'
if file_test(outdir+'resampled',/directory) eq 0 then file_mkdir,outdir+'resampled'
if file_test(outdir+'stacked',/directory) eq 0 then file_mkdir,outdir+'stacked'

; Loop through the fields and construct the list of files
cinfiles = strarr(2e6)
coutfiles = strarr(2e6)
cnt = 0LL
push,reddirarr
for i=0,ndr1-1 do begin
  print,strtrim(i+1,2),'/'+strtrim(ndr1,2),' ',dr1[i].field
  chstr = mrdfits(catdir+'catalogs/final/v4/'+dr1[i].field+'_combined_chips.fits.gz',1)
  chstr.field = strtrim(chstr.field,2)
  chstr.base = strtrim(chstr.base,2)
  chstr.night = strtrim(chstr.night,2)
  chstr.expnum = strtrim(chstr.expnum,2)
  chstr.refexpnum = strtrim(chstr.refexpnum,2)
  fstr = mrdfits(catdir+'catalogs/final/v4/'+dr1[i].field+'_combined_exposures.fits.gz',1)

  ; Unique exposures
  uiexp = uniq(chstr.expnum,sort(chstr.expnum))
  uexpnum = chstr[uiexp].expnum
  nexpnum = n_elements(uexpnum)

  ; Loop through the exposures
  undefine,frawstr,finstcalstr,fresampstr,fstackstr
  For j=0,nexpnum-1 do begin

    iexpnum = uexpnum[j]

    ; Raw files
    gdraw = where(rawstr.expnum eq iexpnum,ngdraw)
    if ngdraw eq 0 then stop,'No raw file for ',dr1[i].field,' ',iexpnum
    push,frawstr,rawstr[gdraw]
    ;dec059577.fits.fz -> ../../../Volumes/archive//mtn/20130317/ct4m/2013A-0411/dec059577.fits.fz
    ;/net/mss1/archive/mtn/20130317/ct4m/2013A-0411/

    ; InstCal files
    gdinstcal = where(instcalstr.expnum eq iexpnum,ngdinstcal)
    if ngdinstcal eq 0 then stop,'No InstCal files for ',dr1[i].field,' ',iexpnum
    ; multiple reduction version
    if ngdinstcal gt 3 then begin
      maxplver = max(instcalstr[gdinstcal].plver)  ; will return "maximum" string
      gdinstcal = where(instcalstr.expnum eq iexpnum and instcalstr.plver eq maxplver,ngdinstcal)
      if ngdinstcal gt 3 then stop,'Too many files for ',dr1[i].field,' ',iexpnum
    endif
    push,finstcalstr,instcalstr[gdinstcal]

    ; Resampled files
    gdresamp = where(resampstr.expnum eq iexpnum,ngdresamp)
    if ngdresamp eq 0 then stop,'No Resampled files for ',dr1[i].field,' ',iexpnum
    ; multiple reduction version
    if ngdresamp gt 3 then begin
      maxplver = max(resampstr[gdresamp].plver)  ; will return "maximum" string
      gdresamp = where(resampstr.expnum eq iexpnum and resampstr.plver eq maxplver,ngdresamp)
      if ngdresamp gt 3 then stop,'Too many files for ',dr1[i].field,' ',iexpnum
    endif
    push,fresampstr,resampstr[gdresamp]
  Endfor

  ; Stacked files
  gdstack = where(stackstr.smashfield eq dr1[i].field,ngdstack)
  if ngdstack eq 0 then stop,'No stack files for ',dr1[i].field
  fstackstr = stackstr[gdstack]  
  ; Get more detailed FITS header information for these files
  match,fstackstr.fullpath,stackhdstr.filename,ind1,ind2,/sort,count=nmatch
  if nmatch lt ngdstack then stop,'not all matched to stackhdstr'
  fstackstr = fstackstr[ind1]   ; want them in the same order
  fstackhdstr = stackhdstr[ind2]

  ; Create unique "ID"s for the stacks
  gdim = where(strtrim(fstackhdstr.prodtype,2) eq 'image',ngdim)
  stackid = fstackhdstr[gdim].stackid
  uistackid = uniq(stackid,sort(stackid))
  ustackid = stackid[uistackid]
  nustackid = n_elements(ustackid)

  ; For each unique stack ID get the group ID with the largest PLVER
  undefine,useind
  for j=0,nustackid-1 do begin
    ind1 = where(fstackhdstr.stackid eq ustackid[j] and strtrim(fstackhdstr.prodtype,2) eq 'image',nind1)
    ; Duplicates, pick the most recent version
    if nind1 gt 1 then begin
      vers = strtrim(fstackhdstr[ind1].sb_dir1,2)+'-'+strtrim(fstackhdstr[ind1].plver,2)
      bestind = first_el(maxloc(vers))
      ; Find the groupID for the latest one
      ;   and add all the associated files to the final list
      groupid1 = fstackhdstr[ind1[bestind]].groupid
      grpind = where(fstackhdstr.groupid eq groupid1,ngrpind)
      push,useind,grpind
    endif
  endfor
  ; Only keep the unique ones
  fstackstr = fstackstr[useind]
  fstackhdstr = fstackhdstr[useind]

; MAYBE LOOP THROUGH THE RELEVANT NIGHTS!!!  

  ; Create empty versions of the files in the right directory
  ;   structure on bambam
  ; Also create summary/inventory files for each directory
  ; one "global" one for each "main" directory.

  ; Make raw directory and files
  if file_test(outdir+'raw/'+dr1[i].field,/directory) eq 0 then file_mkdir,outdir+'raw/'+dr1[i].field
  for j=0,n_elements(frawstr)-1 do spawn,['touch',outdir+'raw/'+dr1[i].field+'/'+frawstr[j].nsaname],/noshell
  ; inventory file
  rawinvfile = outdir+'raw/'+dr1[i].field+'/INVENTORY.raw.'+dr1[i].field+'.txt'
  fmt = '(I-5,A-36,A-10,A-11,A-8,A-12,A-15,A-4,A-25,F-9.2,F-13.6,F-13.6)'
  dateobs = repstr(strtrim(frawstr.date_obs,2),' ','T')   ; replace space with T
  writecol,rawinvfile,lindgen(n_elements(frawstr))+1,frawstr.nsaname,frawstr.expnum,frawstr.proctype,frawstr.prodtype,frawstr.smashfield,$
           frawstr.object,strmid(frawstr.filter,0,1),dateobs,frawstr.exposure,frawstr.ra,frawstr.dec,fmt=fmt
  head = '#NUM   FILENAME                          EXPNUM    PROCTYPE PRODTYPE  SMASHFIELD  OBJECT       FILTER        DATEOBS          EXPTIME      RA           DEC'
  writeline,rawinvfile,head,/prepend

  ; Make instcal directory and files
  if file_test(outdir+'instcal/'+dr1[i].field,/directory) eq 0 then file_mkdir,outdir+'instcal/'+dr1[i].field
  for j=0,n_elements(finstcalstr)-1 do spawn,['touch',outdir+'instcal/'+dr1[i].field+'/'+finstcalstr[j].nsaname],/noshell
  ; inventory file
  instcalinvfile = outdir+'instcal/'+dr1[i].field+'/INVENTORY.instcal.'+dr1[i].field+'.txt'
  fmt = '(I-5,A-36,A-10,A-11,A-8,A-12,A-15,A-4,A-25,F-9.2,F-13.6,F-13.6)'
  dateobs = repstr(strtrim(finstcalstr.date_obs,2),' ','T')   ; replace space with T
  writecol,instcalinvfile,lindgen(n_elements(finstcalstr))+1,finstcalstr.nsaname,finstcalstr.expnum,finstcalstr.proctype,finstcalstr.prodtype,finstcalstr.smashfield,$
           finstcalstr.object,strmid(finstcalstr.filter,0,1),dateobs,finstcalstr.exposure,finstcalstr.ra,finstcalstr.dec,fmt=fmt
  head = '#NUM   FILENAME                          EXPNUM    PROCTYPE PRODTYPE  SMASHFIELD  OBJECT       FILTER        DATEOBS          EXPTIME      RA           DEC'
  writeline,instcalinvfile,head,/prepend

  ; Make resampled directory and files
  if file_test(outdir+'resampled/'+dr1[i].field,/directory) eq 0 then file_mkdir,outdir+'resampled/'+dr1[i].field
  for j=0,n_elements(fresampstr)-1 do spawn,['touch',outdir+'resampled/'+dr1[i].field+'/'+fresampstr[j].nsaname],/noshell
  ; inventory file
  resampinvfile = outdir+'resampled/'+dr1[i].field+'/INVENTORY.resampled.'+dr1[i].field+'.txt'
  fmt = '(I-5,A-36,A-10,A-11,A-8,A-12,A-15,A-4,A-25,F-9.2,F-13.6,F-13.6)'
  dateobs = repstr(strtrim(fresampstr.date_obs,2),' ','T')   ; replace space with T
  writecol,resampinvfile,lindgen(n_elements(fresampstr))+1,fresampstr.nsaname,fresampstr.expnum,fresampstr.proctype,fresampstr.prodtype,fresampstr.smashfield,$
           fresampstr.object,strmid(fresampstr.filter,0,1),dateobs,fresampstr.exposure,fresampstr.ra,fresampstr.dec,fmt=fmt
  head = '#NUM   FILENAME                          EXPNUM    PROCTYPE PRODTYPE  SMASHFIELD  OBJECT       FILTER        DATEOBS          EXPTIME      RA           DEC'
  writeline,resampinvfile,head,/prepend

; GET THE NIGHT INFORMATION FROM THE RAW FILE PATH

  ; Make stacked directory and files
  if file_test(outdir+'stacked/'+dr1[i].field,/directory) eq 0 then file_mkdir,outdir+'stacked/'+dr1[i].field
  for j=0,n_elements(fstackstr)-1 do spawn,['touch',outdir+'stacked/'+dr1[i].field+'/'+fstackstr[j].nsaname],/noshell
  ; inventory file
  stackinvfile = outdir+'stacked/'+dr1[i].field+'/INVENTORY.stacked.'+dr1[i].field+'.txt'
  fmt = '(I-5,A-36,A-11,A-8,A-12,A-15,A-4,A-25,F-9.2,F-13.6,F-13.6)'
  dateobs = repstr(strtrim(fstackstr.date_obs,2),' ','T')   ; replace space with T
  writecol,stackinvfile,lindgen(n_elements(fstackstr))+1,fstackstr.nsaname,fstackstr.proctype,fstackstr.prodtype,fstackstr.smashfield,$
           fstackstr.object,strmid(fstackstr.filter,0,1),dateobs,fstackstr.exposure,fstackstr.ra,fstackstr.dec,fmt=fmt
  head = '#NUM   FILENAME                          PROCTYPE PRODTYPE  SMASHFIELD  OBJECT       FILTER        DATEOBS          EXPTIME      RA           DEC'
  writeline,stackinvfile,head,/prepend

  stop

endfor  ; field

stop

end
