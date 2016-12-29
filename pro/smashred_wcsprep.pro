pro smashred_wcsprep,fstr,refcatname=refcatname,redo=redo

if n_elements(fstr) eq 0 then begin
  print,'Syntax - smashred_wcsprep,fstr,refcatname=refcatname,redo=redo'
  return
endif

ui = uniq(fstr.object,sort(fstr.object))
fieldstr = fstr[ui]  ; unique fields 

if n_elements(refcatname) eq 0 then refcatname='USNO-B1'
;if n_elements(refcatname) eq 0 then refcatname='2MASS-PSC'

nfields = n_elements(fieldstr)
for i=0,nfields-1 do begin
  print,strtrim(i+1,2),'/',strtrim(nfields,2),'  ',fieldstr[i].object

  ra = double(sexig2ten(fieldstr[i].ra))*15
  dec = double(sexig2ten(fieldstr[i].dec))
  outfile = fieldstr[i].object+'_refcat.dat'
  if file_test(outfile) eq 0 or keyword_set(redo) then begin
    print,'Querying VizieR for ',refcatname,' sources'
    cfa = 1
    userefcatname = refcatname
    if refcatname eq '2MASS-PSC' and cfa eq 1 then userefcatname='II/246'  ; cfa issue
    refcatall = queryvizier(userefcatname,[ra,dec],1.5*60.,cfa=cfa,/allcolumns)
    if size(refcatall,/type) ne 8 then begin
      print,'Problem with queryvizier.  Bombing'
      goto,bomb0
    endif
    print,strtrim(n_elements(refcatall),2),' 2MASS sources'
    print,'Saving to ',outfile
    save,refcatall,file=outfile
  endif else begin
    print,'Loading previously-saved ',outfile
    restore,outfile
  endelse

  if (strupcase(refcatname) eq 'GAIA' or strupcase(refcatname) eq 'GAIA/GAIA') and $
     tag_exist(refcatall,'RAJ2000') eq 0 then begin
    add_tag,refcatall,'raj2000',0.0d0,refcatall
    add_tag,refcatall,'dej2000',0.0d0,refcatall
    refcatall.raj2000 = refcatall.ra_icrs
    refcatall.dej2000 = refcatall.de_icrs
  endif

  ; Loop through all exposures for this field
  expind = where(fstr.object eq fieldstr[i].object,nexpind)
  print,strtrim(nexpind,2),' exposures for this field'
  for j=0,nexpind-1 do begin
    fstr1 = fstr[expind[j]]
    chipfiles = file_search(fstr1.newfile+'_[0-9][0-9].fits',count=nchipfiles)
    ; Loop through the chip files for this exposure/field
    for k=0,nchipfiles-1 do begin
      outfile1 = file_basename(chipfiles[k],'.fits')+'_refcat.dat'
      if file_test(outfile1) eq 0 or keyword_set(redo) then begin
        head = headfits(chipfiles[k])
        nx = sxpar(head,'naxis1')
        ny = sxpar(head,'naxis2')
        ;ra1 = double(sexig2ten(sxpar(head,'RA')))
        ;dec1 = double(sexig2ten(sxpar(head,'DEC')))
        head_xyad,head,nx/2,ny/2,ra1,dec1,/deg,error=wcserror
        if n_elements(wcserror) ne 0 then begin
          print,'WCS error for ',chipfiles[k],'.  Skipping.'
          goto,BOMB
        endif
        ; just get sources in a box
        ;  the chips are ~19' wide E/W and  ~11' high N/S, want catalog
        ;  a little bit wider
        gdrefcat = where(abs(refcatall.raj2000-ra1) lt 15./60./cos(dec1/!radeg) and $
                        abs(refcatall.dej2000-dec1) lt 10./60.,ngdrefcat)
        refcat = refcatall[gdrefcat]

        ; Add Rmag, Bmag, Rerr, Berr
        if (refcatname eq 'USNO-B1') then begin

          ; Getting average Rmag, Bmag 
          rmag = median([[refcat.r1mag],[refcat.r2mag]],dim=2,/even)
          bmag = median([[refcat.b1mag],[refcat.b2mag]],dim=2,/even)
          rmagbd = where(finite(rmag) eq 0,nrmagbd)      ; replace NaNs with 99.9999       
          if nrmagbd gt 0 then rmag[rmagbd] = 99.9999
          bmagbd = where(finite(bmag) eq 0,nbmagbd)
          if nbmagbd gt 0 then bmag[bmagbd] = 99.9999
    
          ; Getting Err, Berr
          rerr = abs(refcat.r1mag-refcat.r2mag)
          berr = abs(refcat.b1mag-refcat.b2mag)
          rerrbd = where(finite(rerr) eq 0,nrerrbd)      ; replace NaNs with 0.2
          if nrerrbd gt 0 then rerr[rerrbd] = 0.2
          if nrmagbd gt 0 then rerr[rmagbd] = 9.9999     ; bad mag -> bad error
          berrbd = where(finite(berr) eq 0,nberrbd)
          if nberrbd gt 0 then berr[berrbd] = 0.2
          if nbmagbd gt 0 then berr[bmagbd] = 9.9999     ; bad mag -> bad error

          ; Adding the BMAG, RMAG, BERR, RERR tags
          ADD_TAG,refcat,'BMAG',0.0,refcat
          ADD_TAG,refcat,'RMAG',0.0,refcat
          ADD_TAG,refcat,'BERR',0.0,refcat
          ADD_TAG,refcat,'RERR',0.0,refcat
          refcat.bmag = bmag
          refcat.rmag = rmag
          refcat.berr = berr
          refcat.rerr = rerr

          ; Remove stars with bad photometry
          gd = where(refcat.bmag lt 50. and refcat.rmag lt 50. and $
                     refcat.berr lt 2. and refcat.rerr lt 2.,ngd)
          refcat = refcat[gd]
        endif

        print,'  ',strtrim(i+1,2),' ',strtrim(j+1,2),' ',strtrim(k+1,2),' ',strtrim(ngdrefcat,2),' ',refcatname,' sources'
        print,'  Saving to ',outfile1
        save,refcat,file=outfile1
      endif else print,'  ',strtrim(i+1,2),' ',strtrim(j+1,2),' ',strtrim(k+1,2),' ',outfile1,' already exists'
      BOMB:
    endfor ; chip loop
  endfor ; exposures loop
  BOMB0:
endfor ; field loop

;stop

end
