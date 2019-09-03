pro download_gaiadr2,field,redo=redo

outdir = '/dl1/users/dnidever/smash/cp/red/photred/gaiadr2/'
outfile = outdir+field+'_gaiadr2.fits'
if (file_test(outfile) eq 1 or file_test(outfile+'.gz') eq 1) and not keyword_set(redo) then begin
  print,outfile,' EXISTS and /redo not set'
  return
endif

;; Normal SMASH field
if valid_num(field) eq 0 then begin
  smashred_getredinfo,info,/silent
  ind = where(info.field eq field,nind)
  info1 = info[ind]

  ; Load the summary chip files
  undefine,chips
  for i=0,n_elements(info1)-1 do begin
    chips1 = mrdfits(info1[i].file,2,/silent)
    push,chips,chips1
  endfor

  cendec = mean(minmax(chips.dec))
  decr = range(chips.dec)*1.1 > 2.3
  if range(chips.ra) gt 180 then begin
    ra = chips.ra
    over = where(ra gt 180,nover,comp=under,ncomp=nunder)
    if nover gt 0 then ra[over]-=360
    cenra = mean(minmax(ra))
    rar = range(ra)*cos(cendec/!radeg)*1.1 > 2.3
  endif else begin
    cenra = mean(minmax(chips.ra))
    rar = range(chips.ra)*cos(cendec/!radeg)*1.1 > 2.3
  endelse
  print,field,' ',cenra,cendec,rar,decr
  rad = 1.8

;; HEALpix
endif else begin
  nside= 64
  radeg= 180.0d0 / !dpi
  PIX2ANG_RING,nside,long(field),centheta,cenphi
  cenra = first_el(cenphi*radeg)
  cendec = first_el(90-centheta*radeg)
  ;; resolution is
  ;; ~55'                                                                                                                                                              
  rad = 1.5
endelse

;pylines = 'python -c "from dl import authClient as ac, queryClient as qc, helpers;'+$
;          "token=ac.login('dnidever');"+$
;          "qc.set_timeout_request(10000);"+$
;          "res=qc.query(token,sql='select * from gaia_dr2.gaia_source where q3c_radial_query(ra,dec,"+$
;           stringize(cenra,ndec=5)+","+stringize(cendec,ndec=5)+","+stringize(rad,ndec=3)+")');"+$
;          "df = helpers.utils.convert(res,'table'); df.write('"+outfile+"')"+'"'
;spawn,pylines,out,errout
;if file_test(outfile) eq 0 then begin
;  print,outfile,' NOT FOUND. Something went wrong.'
;  printline,out
;  printine,errout
;  return
;endif
;print,outfile


tablename = 'gaia_dr2.gaia_source'
cols = 'source_id,ra,ra_error,dec,dec_error,pmra,pmra_error,pmdec,pmdec_error,phot_g_mean_flux,phot_g_mean_flux_error,'+$
       'phot_g_mean_mag,phot_bp_mean_mag,phot_bp_mean_flux,phot_bp_mean_flux_error,'+$
       'phot_rp_mean_mag,phot_rp_mean_flux,phot_rp_mean_flux_error'
server = 'gp04.datalab.noao.edu'
; Use Postgres command with q3c cone search
tmpfile = MKTEMP('gaiadr2')
cmd = "psql -h "+server+" -U datalab -d tapdb -w --pset footer -c 'SELECT "+cols+" FROM "+tablename+$
      " WHERE q3c_radial_query(ra,dec,"+stringize(cenra,ndec=4,/nocomma)+","+stringize(cendec,ndec=4,/nocomma)+$
      ","+stringize(rad,ndec=3)+")' > "+tmpfile
file_delete,tmpfile,/allow
spawn,cmd,out,outerr
; Check for empty query
READLINE,tmpfile,tlines,nlineread=4
if n_elements(tlines) lt 4 then begin
  if not keyword_set(silent) then print,'No Results'
endif
;;  Load ASCII file and create the FITS file
ref = importascii(tmpfile,/header,delim='|',skipline=2,/silent)
print,'Saving catalog to file '+outfile
MWRFITS,ref,outfile,/create
file_delete,tmpfile,/allow

print,strtrim(n_elements(ref),2),' GAIADR2 sources'

;; Compress the file
spawn,['gzip',outfile],/noshell

end
