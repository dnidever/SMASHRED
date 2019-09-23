pro dr2_cpsymlinks,data

; Figure out the CP files we want to release and sym links to produce


catdir = '/dl1/users/dnidever/smash/cp/red/photred/'
;indir = '/data/smash/cp/red/photred/'
;outdir = '/data/smash/dr1/'
;indir = '/datalab/users/dnidever/smash/cp/red/photred/'
;outdir = '/datalab/users/dnidever/smash/dr1/'
outdir = '/dl1/users/dnidever/smash/dr2/'

; Restore the DR2 fields
dr2 = importascii('~/projects/SMASHRED/data/smash_fields_final.txt',/header)
add_tag,dr2,'field','',dr2
dr2.field = 'Field'+strtrim(dr2.num,2)
ndr2 = n_elements(dr2)

; Get smash reduction info
;smashred_getredinfo,allredinfo,/silent
;bd = where(strtrim(allredinfo.field,2) eq '',nbd)
;if nbd gt 0 then remove,bd,allredinfo

; Load the full file paths
;;fpath = mrdfits('alldec2.fits',1)
;fpath = mrdfits(catdir+'observations/mss_filenames.fits.gz',1)
;;for i=0,n_tags(fpath)-1 do fpath.(i)=strtrim(fpath.(i),2)
;fpath.base = strtrim(fpath.base,2)

;; Load all the DECam information
if n_elements(data) eq 0 then begin
  print,'Loading DECam Archive information'
  data = mrdfits('/dl1/users/dnidever/nsc/instcal/v3/lists/decam_archive_info_070819.fits.gz',1)
  data.prodtype = strtrim(data.prodtype,2)
  data.proctype = strtrim(data.proctype,2)
  data.dtacqnam = strtrim(data.dtacqnam,2)
endif
data_expnum = strmid(file_basename(data.dtacqnam,'.fits.fz'),6)

;; Get SMASH exposures
expstr = mrdfits('/dl1/users/dnidever/smash/cp/red/photred/catalogs/dr2/smash_dr2_exposure.fits',1)

;; VOT of PropIDs from the NOAO Science Archive
nsa1 = mrdfits('/dl1/users/dnidever/smash/cp/red/photred/observations/nsa_2013A-0411.vot.fits.gz',1)  ; Nidever, LAF proposal
nsa2 = mrdfits('/dl1/users/dnidever/smash/cp/red/photred/observations/nsa_2013A-0719.vot.fits.gz',1)  ; Saha
;; only keep the few SMASH fields
g=where(nsa2.ra lt 50 and nsa2.dec lt -50,ng)
nsa2 = nsa2[g]
nsa3 = mrdfits('/dl1/users/dnidever/smash/cp/red/photred/observations/nsa_2013B-0440.vot.fits.gz',1)  ; Nidever, main SMASH proposal
nsa4 = mrdfits('/dl1/users/dnidever/smash/cp/red/photred/observations/nsa_2013B-0614.vot.fits.gz',1)  ; Munoz
nsa = [nsa1,nsa2,nsa3,nsa4]
nsa.telescope = strtrim(nsa.telescope,2)
nsa.instrument = strtrim(nsa.instrument,2)
g = where(nsa.telescope eq 'ct4m' and nsa.instrument eq 'decam',ng)
nsa = nsa[g]
add_tag,nsa,'expnum','',nsa
nsa.expnum = strmid(file_basename(nsa.dtacqnam,'.fits.fz'),6)
;; this includes all of the standard star fields as well

;; Get RAW filenames
graw = where(stregex(data.proctype,'raw',/boolean,/fold_case) eq 1 and data.prodtype eq 'image',ngraw)
MATCH,data_expnum[graw],expstr.expnum,ind1,ind2,/sort,count=nmatch
rawstr = data[graw[ind1]]
; missing 00627575

;; InstCal
;;=========
;; Get InstCal filenames
ginstcal_im = where(stregex(data.proctype,'InstCal',/boolean,/fold_case) eq 1 and data.prodtype eq 'image',ninstcal_im)
MATCH,data_expnum[ginstcal_im],expstr.expnum,ind1,ind2,/sort,count=nmatch
instcalstr_im = data[ginstcal_im[ind1]]
; 5982, all

ginstcal_dq = where(stregex(data.proctype,'InstCal',/boolean,/fold_case) eq 1 and data.prodtype eq 'dqmask',ninstcal_dq)
MATCH,data_expnum[ginstcal_dq],expstr.expnum,ind1,ind2,/sort,count=nmatch
instcalstr_dq = data[ginstcal_dq[ind1]]
; 5982, all

ginstcal_wt = where(stregex(data.proctype,'InstCal',/boolean,/fold_case) eq 1 and data.prodtype eq 'wtmap',ninstcal_wt)
MATCH,data_expnum[ginstcal_wt],expstr.expnum,ind1,ind2,/sort,count=nmatch
instcalstr_wt = data[ginstcal_wt[ind1]]
; 5977, 5 missing

instcalstr = [instcalstr_im,instcalstr_dq,instcalstr_wt]


;; Resampled
;;===========
;; Get Resampled filenames
gresamp_im = where(stregex(data.proctype,'Resampled',/boolean,/fold_case) eq 1 and data.prodtype eq 'image',nresamp_im)
MATCH,data_expnum[gresamp_im],expstr.expnum,ind1,ind2,/sort,count=nmatch
resampstr_im = data[gresamp_im[ind1]]
; 5980, 2 missing

gresamp_dq = where(stregex(data.proctype,'Resampled',/boolean,/fold_case) eq 1 and data.prodtype eq 'dqmask',nresamp_dq)
MATCH,data_expnum[gresamp_dq],expstr.expnum,ind1,ind2,/sort,count=nmatch
resampstr_dq = data[gresamp_dq[ind1]]
; 5979, 3 missing

gresamp_wt = where(stregex(data.proctype,'Resampled',/boolean,/fold_case) eq 1 and data.prodtype eq 'wtmap',nresamp_wt)
MATCH,data_expnum[gresamp_wt],expstr.expnum,ind1,ind2,/sort,count=nmatch
resampstr_wt = data[gresamp_wt[ind1]]
; 5976, 6 missing

resampstr = [resampstr_im,resampstr_dq,resampstr_wt]


;; Stacked
;;===========
;; Get Stacked filenames
gstacked_dq = where(stregex(data.proctype,'Stacked',/boolean,/fold_case) eq 1 and data.prodtype eq 'dqmask',ngstacked_dq)
MATCH,data_expnum[gstacked_dq],expstr.expnum,ind1,ind2,/sort,count=nmatch
stackstr_dq = data[gstacked_dq[ind1]]
; 1507

gstacked_exp = where(stregex(data.proctype,'Stacked',/boolean,/fold_case) eq 1 and data.prodtype eq 'expmap',ngstacked_exp)
MATCH,data_expnum[gstacked_exp],expstr.expnum,ind1,ind2,/sort,count=nmatch
stackstr_exp = data[gstacked_exp[ind1]]
; 1506 

gstacked_im = where(stregex(data.proctype,'Stacked',/boolean,/fold_case) eq 1 and data.prodtype eq 'image',ngstacked_im)
MATCH,data_expnum[gstacked_im],expstr.expnum,ind1,ind2,/sort,count=nmatch
stackstr_im = data[gstacked_im[ind1]]
; 1501

gstacked_im1 = where(stregex(data.proctype,'Stacked',/boolean,/fold_case) eq 1 and data.prodtype eq 'image1',ngstacked_im1)
MATCH,data_expnum[gstacked_im1],expstr.expnum,ind1,ind2,/sort,count=nmatch
stackstr_im1 = data[gstacked_im1[ind1]]
; 1454

gstacked_wt = where(stregex(data.proctype,'Stacked',/boolean,/fold_case) eq 1 and data.prodtype eq 'wtmap',ngstacked_wt)
MATCH,data_expnum[gstacked_wt],expstr.expnum,ind1,ind2,/sort,count=nmatch
stackstr_wt = data[gstacked_wt[ind1]]
; 1505

stackstr = [stackstr_dq,stackstr_exp,stackstr_im,stackstr_im1,stackstr_wt]



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
;stackhdstr = mrdfits(catdir+'observations/smash_stacked_vot_scrape.fits',1)
;stackhdstr.filename = strtrim(stackhdstr.filename,2)
;add_tag,stackhdstr,'minexptime',0.0,stackhdstr
;add_tag,stackhdstr,'maxexptime',0.0,stackhdstr
;add_tag,stackhdstr,'totexptime',0.0,stackhdstr
;add_tag,stackhdstr,'exptimebin','',stackhdstr
;
;; Get the minimum and maximum exptime per stack
;for i=0,n_elements(stackhdstr)-1 do begin
;  names = [stackhdstr[i].imcmb001, stackhdstr[i].imcmb002, stackhdstr[i].imcmb003, stackhdstr[i].imcmb004,$
;           stackhdstr[i].imcmb005, stackhdstr[i].imcmb006, stackhdstr[i].imcmb007, stackhdstr[i].imcmb008,$
;           stackhdstr[i].imcmb009, stackhdstr[i].imcmb010]
;  gdnames = where(names ne 'None',ngdnames)
;  names = names[gdnames]
;  expnum = strmid(names,6)
;  exptime = fltarr(ngdnames)-1
;  filter = strarr(ngdnames)
;  match,allrawstr.expnum,expnum,ind1,ind2,/sort,count=nmatch
;  exptime[ind2] = allrawstr[ind1].exposure
;  filter[ind2] = allrawstr[ind1].filter
;  if nmatch ne ngdnames then stop,'not all matched'
;  ;for j=0,ngdnames-1 do begin
;  ;  match,allrawstr.expnum,expnum[j],ind1,ind2,/sort,count=nmatch
;  ;  if nmatch gt 0 then exptime[j]=allrawstr[ind1[0]].exposure
;  ;endfor
;  ; Some stacks are missing filter information from the NSA VO table
;  if strtrim(stackhdstr[i].filter,2) eq 'None' then stackhdstr[i].filter=filter[0]
;  gdexptime = where(exptime gt 0,ngdexptime)
;  stackhdstr[i].minexptime = min(exptime[gdexptime])
;  stackhdstr[i].maxexptime = max(exptime[gdexptime])
;  stackhdstr[i].totexptime = total(exptime[gdexptime]) 
;  ; Decide on exptime bin names
;  ;stop
;
;endfor
;mwrfits,stackhdstr,catdir+'observations/smash_stacked_vot_scrape.fits',/create
stackhdstr = mrdfits(catdir+'observations/smash_stacked_vot_scrape.fits',1)
stackhdstr.filename = strtrim(stackhdstr.filename,2)

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
  gdstack = where(strtrim(stackstr.smashfield,2) eq dr1[i].field,ngdstack)
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

  ; Create empty versions of the files in the right directory
  ;   structure on bambam
  ; Also create summary/inventory files for each directory
  ; one "global" one for each "main" directory.

  ; Make raw directory and files
  if file_test(outdir+'raw/'+dr1[i].field,/directory) eq 0 then file_mkdir,outdir+'raw/'+dr1[i].field
  for j=0,n_elements(frawstr)-1 do spawn,['touch',outdir+'raw/'+dr1[i].field+'/'+frawstr[j].nsaname],/noshell
  ; inventory file
  rawinvfile = outdir+'raw/'+dr1[i].field+'/INVENTORY.raw.'+dr1[i].field+'.txt'
  fmt = '(I-5,A-36,A-10,A-11,A-8,A-12,A-15,A-4,A-10,A-25,F-9.2,F-13.6,F-13.6)'
  dateobs = repstr(strtrim(frawstr.date_obs,2),' ','T')   ; replace space with T
  object = repstr(frawstr.object,' ','')
  writecol,rawinvfile,lindgen(n_elements(frawstr))+1,frawstr.nsaname,frawstr.expnum,frawstr.proctype,frawstr.prodtype,frawstr.smashfield,$
           object,strmid(frawstr.filter,0,1),frawstr.night,dateobs,frawstr.exposure,frawstr.ra,frawstr.dec,fmt=fmt
  head = '#NUM   FILENAME                          EXPNUM    PROCTYPE PRODTYPE  SMASHFIELD  OBJECT       FILTER  NIGHT          DATEOBS          EXPTIME      RA           DEC'
  writeline,rawinvfile,head,/prepend
  ftypes = [3,7,7,7,7,7,7,7,7,7,4,5,5]
  rawinv = importascii(rawinvfile,/header,fieldtypes=ftypes)
  mwrfits,rawinv,repstr(rawinvfile,'.txt','.fits'),/create
  push,allrawinv,rawinv

  ; Make instcal directory and files
  if file_test(outdir+'instcal/'+dr1[i].field,/directory) eq 0 then file_mkdir,outdir+'instcal/'+dr1[i].field
  for j=0,n_elements(finstcalstr)-1 do spawn,['touch',outdir+'instcal/'+dr1[i].field+'/'+finstcalstr[j].nsaname],/noshell
  ; inventory file
  instcalinvfile = outdir+'instcal/'+dr1[i].field+'/INVENTORY.instcal.'+dr1[i].field+'.txt'
  fmt = '(I-5,A-36,A-10,A-11,A-8,A-12,A-15,A-4,A-10,A-25,F-9.2,F-13.6,F-13.6)'
  dateobs = repstr(strtrim(finstcalstr.date_obs,2),' ','T')   ; replace space with T
  object = repstr(finstcalstr.object,' ','')
  writecol,instcalinvfile,lindgen(n_elements(finstcalstr))+1,finstcalstr.nsaname,finstcalstr.expnum,finstcalstr.proctype,finstcalstr.prodtype,finstcalstr.smashfield,$
           object,strmid(finstcalstr.filter,0,1),finstcalstr.night,dateobs,finstcalstr.exposure,finstcalstr.ra,finstcalstr.dec,fmt=fmt
  head = '#NUM   FILENAME                          EXPNUM    PROCTYPE PRODTYPE  SMASHFIELD  OBJECT       FILTER  NIGHT          DATEOBS          EXPTIME      RA           DEC'
  writeline,instcalinvfile,head,/prepend
  ftypes = [3,7,7,7,7,7,7,7,7,7,4,5,5]
  instcalinv = importascii(instcalinvfile,/header,fieldtypes=ftypes)
  mwrfits,instcalinv,repstr(instcalinvfile,'.txt','.fits'),/create
  push,allinstcalinv,instcalinv

  ; Make resampled directory and files
  if file_test(outdir+'resampled/'+dr1[i].field,/directory) eq 0 then file_mkdir,outdir+'resampled/'+dr1[i].field
  for j=0,n_elements(fresampstr)-1 do spawn,['touch',outdir+'resampled/'+dr1[i].field+'/'+fresampstr[j].nsaname],/noshell
  ; inventory file
  resampinvfile = outdir+'resampled/'+dr1[i].field+'/INVENTORY.resampled.'+dr1[i].field+'.txt'
  fmt = '(I-5,A-36,A-10,A-11,A-8,A-12,A-15,A-4,A-10,A-25,F-9.2,F-13.6,F-13.6)'
  dateobs = repstr(strtrim(fresampstr.date_obs,2),' ','T')   ; replace space with T
  object = repstr(fresampstr.object,' ','')
  writecol,resampinvfile,lindgen(n_elements(fresampstr))+1,fresampstr.nsaname,fresampstr.expnum,fresampstr.proctype,fresampstr.prodtype,fresampstr.smashfield,$
           object,strmid(fresampstr.filter,0,1),fresampstr.night,dateobs,fresampstr.exposure,fresampstr.ra,fresampstr.dec,fmt=fmt
  head = '#NUM   FILENAME                          EXPNUM    PROCTYPE PRODTYPE  SMASHFIELD  OBJECT       FILTER  NIGHT          DATEOBS          EXPTIME      RA           DEC'
  writeline,resampinvfile,head,/prepend
  ftypes = [3,7,7,7,7,7,7,7,7,7,4,5,5]
  resampinv = importascii(resampinvfile,/header,fieldtypes=ftypes)
  mwrfits,resampinv,repstr(resampinvfile,'.txt','.fits'),/create
  push,allresampinv,resampinv

  ; Make stacked directory and files
  if file_test(outdir+'stacked/'+dr1[i].field,/directory) eq 0 then file_mkdir,outdir+'stacked/'+dr1[i].field
  for j=0,n_elements(fstackstr)-1 do spawn,['touch',outdir+'stacked/'+dr1[i].field+'/'+fstackstr[j].nsaname],/noshell
  ; inventory file
  stackinvfile = outdir+'stacked/'+dr1[i].field+'/INVENTORY.stacked.'+dr1[i].field+'.txt'
  fmt = '(I-5,A-36,A-11,A-8,A-12,A-15,A-4,A-25,F-9.2,F-13.6,F-13.6)'
  dateobs = repstr(strtrim(fstackstr.date_obs,2),' ','T')   ; replace space with T
  object = repstr(fstackstr.object,' ','')
  writecol,stackinvfile,lindgen(n_elements(fstackstr))+1,fstackstr.nsaname,fstackstr.proctype,fstackstr.prodtype,fstackstr.smashfield,$
           object,strmid(fstackstr.filter,0,1),dateobs,fstackstr.exposure,fstackstr.ra,fstackstr.dec,fmt=fmt
  head = '#NUM   FILENAME                          PROCTYPE PRODTYPE  SMASHFIELD  OBJECT       FILTER        DATEOBS          EXPTIME      RA           DEC'
  writeline,stackinvfile,head,/prepend
  ftypes = [3,7,7,7,7,7,7,7,4,5,5]
  stackinv = importascii(stackinvfile,/header,fieldtypes=ftypes)
  mwrfits,stackinv,repstr(stackinvfile,'.txt','.fits'),/create
  push,allstackinv,stackinv

endfor  ; field

; Make combined inventory files, ASCII and FITS versions

; RAW
rawinvfile = outdir+'raw/INVENTORY.raw.txt'
fmt = '(I-7,A-36,A-10,A-11,A-8,A-12,A-15,A-4,A-10,A-25,F-9.2,F-13.6,F-13.6)'
allrawinv.num = lindgen(n_elements(allrawinv))+1
writecol,rawinvfile,allrawinv.num,allrawinv.filename,allrawinv.expnum,allrawinv.proctype,allrawinv.prodtype,allrawinv.smashfield,$
         allrawinv.object,allrawinv.filter,allrawinv.night,allrawinv.dateobs,allrawinv.exptime,allrawinv.ra,allrawinv.dec,fmt=fmt
head = '# NUM    FILENAME                          EXPNUM    PROCTYPE PRODTYPE  SMASHFIELD  OBJECT       FILTER  NIGHT          DATEOBS          EXPTIME      RA           DEC'
writeline,rawinvfile,head,/prepend
mwrfits,allrawinv,repstr(rawinvfile,'.txt','.fits'),/create

; InstCal
instcalinvfile = outdir+'instcal/INVENTORY.instcal.txt'
fmt = '(I-7,A-36,A-10,A-11,A-8,A-12,A-15,A-4,A-10,A-25,F-9.2,F-13.6,F-13.6)'
allinstcalinv.num = lindgen(n_elements(allinstcalinv))+1
writecol,instcalinvfile,allinstcalinv.num,allinstcalinv.filename,allinstcalinv.expnum,allinstcalinv.proctype,allinstcalinv.prodtype,allinstcalinv.smashfield,$
         allinstcalinv.object,allinstcalinv.filter,allinstcalinv.night,allinstcalinv.dateobs,allinstcalinv.exptime,allinstcalinv.ra,allinstcalinv.dec,fmt=fmt
head = '# NUM    FILENAME                          EXPNUM    PROCTYPE PRODTYPE  SMASHFIELD  OBJECT       FILTER  NIGHT          DATEOBS          EXPTIME      RA           DEC'
writeline,instcalinvfile,head,/prepend
mwrfits,allinstcalinv,repstr(instcalinvfile,'.txt','.fits'),/create

; Resampled
resampinvfile = outdir+'resampled/INVENTORY.resampled.txt'
fmt = '(I-7,A-36,A-10,A-11,A-8,A-12,A-15,A-4,A-10,A-25,F-9.2,F-13.6,F-13.6)'
allresampinv.num = lindgen(n_elements(allresampinv))+1
writecol,resampinvfile,allresampinv.num,allresampinv.filename,allresampinv.expnum,allresampinv.proctype,allresampinv.prodtype,allresampinv.smashfield,$
         allresampinv.object,allresampinv.filter,allresampinv.night,allresampinv.dateobs,allresampinv.exptime,allresampinv.ra,allresampinv.dec,fmt=fmt
head = '# NUM    FILENAME                          EXPNUM    PROCTYPE PRODTYPE  SMASHFIELD  OBJECT       FILTER  NIGHT          DATEOBS          EXPTIME      RA           DEC'
writeline,resampinvfile,head,/prepend
mwrfits,allresampinv,repstr(resampinvfile,'.txt','.fits'),/create

; Stacked
rawinvfile = outdir+'stacked/INVENTORY.stacked.txt'
fmt = '(I-7,A-36,A-11,A-8,A-12,A-15,A-4,A-25,F-9.2,F-13.6,F-13.6)'
allstackinv.num = lindgen(n_elements(allstackinv))+1
writecol,stackinvfile,allstackinv.num,allstackinv.filename,allstackinv.proctype,allstackinv.prodtype,allstackinv.smashfield,$
         allstackinv.object,allstackinv.filter,allstackinv.dateobs,allstackinv.exptime,allstackinv.ra,allstackinv.dec,fmt=fmt
head = '# NUM    FILENAME                          PROCTYPE PRODTYPE  SMASHFIELD  OBJECT       FILTER        DATEOBS          EXPTIME      RA           DEC'
writeline,stackinvfile,head,/prepend
mwrfits,stackinv,repstr(stackinvfile,'.txt','.fits'),/create


stop

end
