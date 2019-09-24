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

; Calibration information
calstr = mrdfits(catdir+'check_calibrated_v6.fits',1)
calstr.field = strtrim(calstr.field,2)
MATCH,dr2.field,calstr.field,ind1,ind2,/sort,count=nmatch
dr2 = dr2[ind1]  ; 197 fields
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
  data.date_obs = strtrim(data.date_obs,2)
  data.filter = strtrim(data.filter,2)
  data.uri = strtrim(data.uri,2)
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
nsa.date_obs = strtrim(nsa.date_obs,2)
nsa.proctype = strtrim(nsa.proctype,2)
nsa.prodtype = strtrim(nsa.prodtype,2)
nsa.filter = strtrim(nsa.filter,2)
nsa.archive_file = strtrim(nsa.archive_file,2)
g = where(nsa.telescope eq 'ct4m' and nsa.instrument eq 'decam',ng)
nsa = nsa[g]
add_tag,nsa,'expnum','',nsa
nsa.expnum = strmid(file_basename(nsa.dtacqnam,'.fits.fz'),6)
;; this includes all of the standard star fields as well

;; RAW
;;======
print,'Getting RAW'
;; Get RAW filenames
graw = where(stregex(data.proctype,'raw',/boolean,/fold_case) eq 1 and data.prodtype eq 'image',ngraw)
MATCH,data_expnum[graw],expstr.expnum,ind1,ind2,/sort,count=nmatch
rawstr = data[graw[ind1]]
; missing 00627575

;; InstCal
;;=========
print,'Getting InstCal'
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
print,'Getting Resampled'
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
print,'Getting Stacked'
;; Get stacked images near the SMASH fields
gstack = where(stregex(nsa.proctype,'stack',/boolean,/fold_case) eq 1,ngstack)
nsastack = nsa[gstack]
undefine,stackind
for i=0,ndr2-1 do begin
  dist = sphdist(dr2[i].radeg,dr2[i].dedeg,nsastack.ra,nsastack.dec,/deg)
  gd = where(dist lt 2.0,ngd)
  if ngd gt 0 then push,stackind,gd
  print,strtrim(i+1,2),' ',dr2[i].field,' ',strtrim(ngd,2)
endfor
; 5 files per exposure bin and filter
ui = uniq(stackind,sort(stackind))
stackind = stackind[ui]
nsastackstr = nsastack[stackind]

;dataid = data.date_obs+'-'+data.proctype+'-'+data.prodtype+'-'+strmid(data.filter,0,1)+'-'+strtrim(string(data.exposure,format='(F10.2)'),2)
;nsaid = nsastackstr.date_obs+'-'+nsastackstr.proctype+'-'+nsastackstr.prodtype+'-'+strmid(nsastackstr.filter,0,1)+'-'+strtrim(string(nsastackstr.exposure,format='(F10.2)'),2)
datafile = file_basename(data.uri)
match,datafile,nsastackstr.archive_file,ind1,ind2,/sort,count=nmatch
stackstr = data[ind1]

;; Combine all lists
filestr = [rawstr, instcalstr, resampstr, stackstr]
add_tag,filestr,'fullpath','',filestr
filestr.fullpath = repstr(strtrim(filestr.uri,2),'irods:///noao-tuc-z1/','/net/mss1/archive/')

;mwrfits,filestr,catdir+'observations/smash_dr2_images.fits',/create

stop

end
