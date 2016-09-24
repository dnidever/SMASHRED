pro smashred_wcsprep_gaia,fstr,redo=redo

;if n_elements(fstr) eq 0 then begin
;  print,'Syntax - smashred_wcsprep_gaia,fstr,redo=redo'
;  return
;endif

smashred_getredinfo,fstr,/silent

bd = where(strtrim(fstr.field,2) eq '',nbd)
if nbd gt 0 then remove,bd,fstr
ui = uniq(fstr.field,sort(fstr.field))
fields = fstr[ui].field
nfields = n_elements(fields)

refcatname='GAIA/GAIA'

for i=0,nfields-1 do begin
  print,strtrim(i+1,2),'/',strtrim(nfields,2),'  ',fields[i]
  ind = where(fstr.field eq fields[i],nind)

  ;ra = double(sexig2ten(fieldstr[i].ra))*15
  ;dec = double(sexig2ten(fieldstr[i].dec))
  outfile = '/data/smash/cp/red/photred/gaia/'+fields[i]+'_gaia.fits'

  if file_test(outfile) eq 0 then begin
    ;print,'Querying VizieR for ',refcatname,' sources'
    ;canada = 0 ;1
    ;refcatall = queryvizier(refcatname,[ra,dec],1.5*60.,canada=canada,/allcolumns)
    ;print,strtrim(n_elements(refcatall),2),' 2MASS sources'
    ;print,'Saving to ',outfile
    ;save,refcatall,file=outfile
    stop
  endif else begin
    print,'Loading previously-saved ',outfile
    ;restore,outfile
    refcatall = mrdfits(outfile,1)
    add_tag,refcatall,'raj2000',0.0d0,refcatall
    add_tag,refcatall,'dej2000',0.0d0,refcatall
    refcatall.raj2000 = refcatall.ra_icrs
    refcatall.dej2000 = refcatall.de_icrs
  endelse

  ; Load the chip summary files
  undefine,chips
  for j=0,nind-1 do begin
    chips1 = mrdfits(fstr[ind[j]].file,2,/silent)
    add_tag,chips1,'dir','',chips1
    chips1.dir = file_dirname(fstr[ind[j]].file)
    push,chips,chips1
  endfor
  nchips = n_elements(chips)

  ; Loop through all chips for this field
  print,strtrim(nchips,2),' chips for this field'
  for j=0,nchips-1 do begin
    chipfile = chips[j].dir+'/'+chips[j].field+'/'+chips[j].file
    outfile1 = file_dirname(chipfile)+'/'+file_basename(chipfile,'.fits')+'_refcat_gaia.dat'
    if file_test(outfile1) eq 0 or keyword_set(redo) then begin
      head = headfits(chipfile)
      nx = sxpar(head,'naxis1')
      ny = sxpar(head,'naxis2')
      ;ra1 = double(sexig2ten(sxpar(head,'RA')))
      ;dec1 = double(sexig2ten(sxpar(head,'DEC')))
      head_xyad,head,nx/2,ny/2,ra1,dec1,/deg,error=wcserror
      if n_elements(wcserror) ne 0 then begin
        print,'WCS error for ',chipfile,'.  Skipping.'
        goto,BOMB
      endif
      ; just get sources in a box
      ;  the chips are ~19' wide E/W and  ~11' high N/S, want catalog
      ;  a little bit wider
      gdrefcat = where(abs(refcatall.ra_icrs-ra1) lt 15./60./cos(dec1/!radeg) and $
                      abs(refcatall.de_icrs-dec1) lt 10./60.,ngdrefcat)
      refcat = refcatall[gdrefcat]

      print,'  ',strtrim(i+1,2),' ',strtrim(j+1,2),' ',strtrim(ngdrefcat,2),' ',refcatname,' sources'
      print,'  Saving to ',outfile1
      save,refcat,file=outfile1
    endif else print,'  ',strtrim(i+1,2),' ',strtrim(j+1,2),' ',outfile1,' already exists'
    BOMB:
  endfor ; chip loop
endfor ; field loop

;stop

end
