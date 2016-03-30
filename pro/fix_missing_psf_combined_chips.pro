pro fix_missing_psf_combined_chips

; Fix the missing PSF information in the combined_chips.fits files

dir = '/data/smash/cp/red/photred/catalogs/'

smashred_getredinfo,redinfo,/silent

; Loop through the combined fields
cmbfiles = file_search(dir+'inst/comb/Field*_combined_chips.fits',count=ncmbfiles)
for i=0,ncmbfiles-1 do begin

  print,strtrim(i+1,2),' ',cmbfiles[i]
  field = file_basename(cmbfiles[i],'_combined_chips.fits')
  cmbdir = file_dirname(cmbfiles[i])
  cmbbase = file_basename(cmbfiles[i])
  cmbstr = mrdfits(cmbfiles[i],1)

  bdpsf = where(strtrim(cmbstr.dao_psftype,2) eq 'NAN' or cmbstr.dao_psfboxsize le 0,nbdpsf)
  if nbdpsf gt 0 then begin
    print,strtrim(nbdpsf,2),' bad psfs'

    ; Load ALL of the individual catalog files
    sumind = where(redinfo.field eq field,nsumind)
    undefine,sumstr
    for j=0,nsumind-1 do begin
      sumstr1 = mrdfits(redinfo[sumind[j]].file,2,/silent)
      push,sumstr,sumstr1
    endfor
    undefine,sumstr1

    ; Match up the bad chips
    for j=0,nbdpsf-1 do begin
      ind = where(sumstr.base eq cmbstr[bdpsf[j]].base,nind)
      ind = ind[0]
      cmbstr[bdpsf[j]].dao_psftype = sumstr[ind].dao_psftype
      cmbstr[bdpsf[j]].dao_psfboxsize = sumstr[ind].dao_psfboxsize 
      cmbstr[bdpsf[j]].dao_psfvarorder = sumstr[ind].dao_psfvarorder
      print,sumstr[ind].base,' ',sumstr[ind].dao_psftype,' ',sumstr[ind].dao_psfboxsize,$
            sumstr[ind].dao_psfvarorder
      ;stop
    endfor ; bad psfs

    ; Copy original to backup
    bakfile = cmbdir+'/bak/'+cmbbase
    if file_test(bakfile) eq 0 then file_copy,cmbfiles[i],bakfile,/verbose
    print,'Rewriting ',cmbbase
    MWRFITS,cmbstr,cmbfiles[i],/create

    ;stop

  endif  ; there are some bad psfs

  ;stop

endfor  ; combined fields/files

stop

end
