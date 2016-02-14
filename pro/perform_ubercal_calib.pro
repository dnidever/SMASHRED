pro perform_ubercal_calib,field,reduxdir=reduxdir,redo=redo

; Perform ubercal calibration on catalog with overlaps

if not keyword_set(reduxdir) then reduxdir='/data/smash/cp/red/photred/'
outdir = reduxdir+'catalogs/inst/comb/'

; Get reduction info
smashred_getredinfo,allinfo

gdinfo = where(allinfo.field eq field,ngdinfo)
if ngdinfo eq 0 then begin
  print,'No reduced exposures for ',field
  return
endif
info = allinfo[gdinfo]

;goto,starthere

; Loop through the catalogs
ninfo = n_elements(info)
undefine,allfstr,allchstr
For c=0,ninfo-1 do begin

  info1 = info[c]
  fstr = *info1.fstr
  fbase = file_basename(info1.file,'_summary.fits')
  print,strtrim(c+1,2),'/',strtrim(ninfo,2),'  ',fbase
  chstr = mrdfits(info1.file,2,/silent)  ; load the chip structure
  nchstr = n_elements(chstr)
  night = info1.night

  print,'Performing UBERCAL calibration on ',fbase

  ; Load in the data
  ;   don't have the "reference" frame information, just search
  ;   for PHOT files
  print,'---Loading the data for ',fbase,'---'
  outfile = reduxdir+'catalogs/inst/comb/'+fbase+'_'+night+'_ubercal.dat'
  if file_test(outfile) eq 0 or keyword_set(redo) then begin
    add_tag,chstr,'refexpnum','',chstr
    add_tag,chstr,'vertices_ra',dblarr(4),chstr
    add_tag,chstr,'vertices_dec',dblarr(4),chstr
    add_tag,chstr,'ndata',-1L,chstr
    add_tag,chstr,'data',ptr_new(),chstr
    photfiles = file_search(reduxdir+night+'/'+chstr[0].field+'/'+chstr[0].field+'-*_??.phot',count=nphotfiles)
    print,strtrim(nphotfiles,2),' PHOT files for ',fbase,' in ',reduxdir+night+'/'+chstr[0].field+'/'
    for i=0,nphotfiles-1 do begin
  
      print,strtrim(i+1,2),' ',photfiles[i]
      phot = importascii(photfiles[i],/header,/silent)
      tags = tag_names(phot)
      phbase = file_basename(photfiles[i],'.phot')
      ichip = long(first_el(strsplit(phbase,'_',/extract),/last))
      dum = first_el(strsplit(phbase,'-',/extract),/last)
      refexpnum = first_el(strsplit(dum,'_',/extract))

      chind = where(chstr.chip eq ichip,nchind)
      for j=0,nchind-1 do begin
        mind = where(tags eq chstr[chind[j]].calib_magname,nmind)
        gd = where(phot.(mind[0]) lt 50,ngd)
        temp = replicate({id:-1L,x:0.0,y:0.0,mag:0.0,err:0.0,chi:0.0,sharp:0.0,flag:-1,prob:-1.0,ra:0.0d0,dec:0.0d0},ngd)
        struct_assign,phot[gd],temp,/nozero
        temp.mag = phot[gd].(mind[0])
        temp.err = phot[gd].(mind[0]+1)   ; assume error is the next column

        chstr[chind[j]].refexpnum = refexpnum
        chstr[chind[j]].ndata = ngd
        chstr[chind[j]].data = ptr_new(temp)  ; save the data

        ; get astrometric vertices from header
        fitsfile = reduxdir+night+'/'+chstr[0].field+'/'+chstr[chind[j]].base+'.fits'
        head = headfits(fitsfile)
        nx = sxpar(head,'NAXIS1')
        ny = sxpar(head,'NAXIS2')
        head_xyad,head,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/degree
        chstr[chind[j]].vertices_ra = vra
        chstr[chind[j]].vertices_dec = vdec

        print,j+1,chstr[chind[j]].expnum,chstr[chind[j]].chip,ngd,format='(I4,A12,I5,I10)'
      endfor
    endfor

    print,'Saving info temporarily to ',outfile
    save,info1,fstr,chstr,file=outfile
  ; save file already exists
  endif else begin
    print,'Restoring previously saved ',outfile
    restore,outfile
  endelse

  ; combine
  push,allfstr,fstr
  push,allchstr,chstr

  ;stop

Endfor ; catalog loop

fstr = allfstr
chstr = allchstr
undefine,allfstr,allchstr

;return

;stop

; photmatch matches up catalogs and finds photometric offsets

; do we want to use the combined CAT as well and 
;  figure out what the indices are for each chip??
;  only do astrometric matching for stars that have
;  detections in that exposure

; but since we are going to combine the short and long fields
; chip-by-chip, maybe we need to rematch/recombine everything
; anyway????

; Do a WHILE loop (for each filter) and figure out which chips overlap
; each other, when there's an overlap figure out the magnitude offset

;stop


; How many unique filters do we have
ui = uniq(fstr.filter,sort(fstr.filter))
ufilter = fstr[ui].filter
nufilter = n_elements(ufilter)

print,strtrim(nufilter,2),' unique filters: ',fstr[ui].filter

;tags = tag_names(cat)
;outcat = cat
; set all averaged mags to 99.99, will recalculate those
; set all individual "calib" mags to 99.99 as well, will
;  copy those in as we go along
;
;for i=0,nufilter-1 do begin
;  magind = where(stregex(tags,'^'+strupcase(ufilter[i])+'MAG',/boolean) eq 1,nmagind)
;  if nmagind gt 0 then for j=0,nmagind-1 do outcat.(magind[j])=99.99
;  errind = where(stregex(tags,'^'+strupcase(ufilter[i])+'[0-9]*ERR',/boolean) eq 1,nerrind)
;  if nerrind gt 0 then for j=0,nerrind-1 do outcat.(errind[j])=9.99
;endfor

; add medoff tag
add_tag,chstr,'magoffset',0.0,chstr

; Loop through the unique filters
FOR i=0,nufilter-1 do begin

  ifilter = ufilter[i]
  ;filtind = where(fstr.filter eq ifilter,nfiltind)
  ;filt_fstr = fstr[filtind]
  chfiltind = where(chstr.filter eq ifilter,nchfiltind)
  chfiltstr = chstr[chfiltind]

  print,'--- FILTER = ',ifilter,' ',strtrim(nchfiltind,2),' nchips'

  ; Calculate all overlaps and magnitude offsets
  overlapstr = replicate({expnum1:'',chip1:0L,expnum2:'',chip2:0L,overlap:-1,magoff:99.99,magoffsig:9.99,nmatch:-1L},nchfiltind,nchfiltind)
  for j=0,nchfiltind-1 do begin
    overlapstr[j,*].expnum1 = chfiltstr[j].expnum
    overlapstr[j,*].chip1 = chfiltstr[j].chip
    overlapstr[*,j].expnum2 = chfiltstr[j].expnum
    overlapstr[*,j].chip2 = chfiltstr[j].chip
  endfor
  for j=0,nchfiltind-1 do begin
    for k=j+1,nchfiltind-1 do begin
      ; must be different exposures to overlap
      if chfiltstr[j].expnum ne chfiltstr[k].expnum then begin
        ; use vertices to check for overlap
        ;   use code from printVisitOverlap.py
        overlap = dopolygonsoverlap(chfiltstr[j].vertices_ra, chfiltstr[j].vertices_dec, chfiltstr[k].vertices_ra, chfiltstr[k].vertices_dec)
        overlapstr[j,k].overlap = overlap
        overlapstr[k,j].overlap = overlap
        ; measure mag offsets
        if overlap eq 1 then begin
          str1 = *chfiltstr[j].data
          str2 = *chfiltstr[k].data
          dcr = 0.5
          ;  Match the two catalogs and get photometric offsets
          photmatch,str1,str2,ind1,ind2,dcr=dcr,magoffset=magoffset,magoffsig=magoffsig,astoffset=astoffset,count=nmatch,/silent
          ; too few matches, increase matching radius
          ;if nmatch lt 5 then $
          ;  photmatch,str1,str2,ind1,ind2,dcr=1.0,magoffset=magoffset,magoffsig=magoffsig,astoffset=astoffset,count=nmatch,/silent
          ;if nmatch lt 5 then $
          ;  photmatch,str1,str2,ind1,ind2,dcr=1.5,magoffset=magoffset,magoffsig=magoffsig,astoffset=astoffset,count=nmatch,/silent
          if nmatch gt 1 and n_elements(magoffset) gt 0 then begin
            magoffsig = magoffsig > 1e-6  ; let lower limit, can sometimes be zero
            overlapstr[j,k].magoff = magoffset
            overlapstr[j,k].magoffsig = magoffsig
            overlapstr[j,k].nmatch = nmatch
            print,overlapstr[j,k].expnum1,overlapstr[j,k].chip1,overlapstr[j,k].expnum2,overlapstr[j,k].chip2,overlapstr[j,k].nmatch,$
                  overlapstr[j,k].magoff,overlapstr[j,k].magoffsig,format='(A10,I5,A10,I5,I7,F10.4,F10.5)'
            ; fill in the reverse situation
            overlapstr[k,j].magoff = -magoffset
            overlapstr[k,j].magoffsig = magoffsig
            overlapstr[k,j].nmatch = nmatch
          endif ; some matches
        endif ; overlap
      endif
    endfor
  endfor

  ;stop

  ; Now iterate to put them all on the same system
  count = 0
  flag = 0
  davg0 = 100
  dmax0 = 100
  fmagoff = fltarr(nchfiltind)
  fmagflag = lonarr(nchfiltind)  ; 1-good, 0-no good
  magoff = overlapstr.magoff  ; initialize
  bd = where(magoff gt 50,nbd,comp=gd)
  if nbd gt 0 then magoff[bd]=!values.f_nan
  WHILE (flag eq 0) do begin

    ; maybe calculate for each chip how (on average) it compares to
    ; its neighbors, then remove that and start again
    medoff = median(magoff,dim=1,/even)
    ; we actually use the weighted mean now, but just in case

    ; update the catalogs with these offsets
    ; the trick is only to go halfway, otherwise we facilate a lot
    ;for k=0,nchfiltind-1 do (*chfiltstr[k].data).mag += medoff[k]*0.5

    ; Compute weighted mean offset and remove
    ;  from the magoffset values
    mnoff = fltarr(nchfiltind)
    for k=0,nchfiltind-1 do begin
      ; get good matches with low sigma
      gd = where(finite(reform(magoff[*,k])) eq 1 and overlapstr[*,k].nmatch gt 3 and overlapstr[*,k].magoffsig lt 0.05,ngd)
;if ngd eq 0 then stop,'no good matches'
      if ngd gt 0 then begin
        ; calculate the weighted mean
        robust_mean,magoff[gd,k],robmean,robsig,sig=overlapstr[gd,k].magoffsig
        if finite(robmean) eq 0 then robmean=medoff[k] ; just in case
        magoff[k,gd] += robmean*0.5
        magoff[gd,k] -= robmean*0.5
        fmagoff[k] += robmean*0.5
        mnoff[k] = robmean
        fmagflag[k] = 1
      endif
    endfor

    ; How much have things changed
    davg = total(abs(mnoff))/nchfiltind
    dmax = max(abs(mnoff))
    print,count,davg,dmax

    ; Do we need to stop
    if count ge 300 or (abs(davg-davg0)/davg0 lt 1e-2 and abs(dmax-dmax0)/dmax0 lt 1e-2) then flag=1

    davg0 = davg
    dmax0 = dmax

    count++

  ENDWHILE

  ; Calculate the relative offset for each chip
  ;  first calculate the relative mag offset from exposure median
  rel_fmagoff = fmagoff
  uiexp = uniq(chfiltstr.expnum,sort(chfiltstr.expnum))
  uexp = chfiltstr[uiexp].expnum
  nuexp = n_elements(uexp)
  for k=0,nuexp-1 do begin
    expind = where(chfiltstr.expnum eq uexp[k],nexpind)
    if nexpind gt 0 then rel_fmagoff[expind] -= median(fmagoff[expind])
  endfor
  uichips = uniq(chfiltstr.chip,sort(chfiltstr.chip))
  uchips = chfiltstr[uichips].chip
  nuchips = n_elements(uchips)
  deltamagoff_chip = fltarr(nuchips)
  for k=0,nuchips-1 do begin
    chind = where(chfiltstr.chip eq uchips[k],nchind)
    if nchind gt 0 then deltamagoff_chip[k] = median(rel_fmagoff[chind])
  endfor

  ; chips with no overlap
  ;   use the median offset for all chips of that exposure
  totoverlap = total( (overlapstr.overlap eq 1 and overlapstr.magoff lt 50),1)
  bdoverlap = where(totoverlap eq 0,nbdoverlap)
  for k=0,nbdoverlap-1 do begin
    expind = where(chfiltstr.expnum eq chfiltstr[bdoverlap[k]].expnum and totoverlap gt 0 and fmagflag eq 1,nexpind)
    chind = where(uchips eq chfiltstr[bdoverlap[k]].chip,nchind)
    fmagoff[bdoverlap[k]] = median(fmagoff[expind]) + deltamagoff_chip[chind]
  endfor

  ; Apply these magnitude offsets to the data
  for k=0,nchfiltind-1 do begin
    (*chfiltstr[k].data).mag += fmagoff[k]
    chfiltstr[k].magoffset = fmagoff[k]
  endfor

  ; Stick the results back into CHSTR
  chstr[chfiltind] = chfiltstr

  ;stop

endfor

; the chip number information isn't preseved on an exposure level
; only the chip number of the first detection is saved
; I might need to load the individual calibrated chip PHOT files
; to get the information that I need.
; I might not even need the "final" catalog.

stop

; --- Now combine all of the data ---

; Unique exposures
uiexp = uniq(chstr.expnum,sort(chstr.expnum))
uexp = chstr[uiexp].expnum
nuexp = n_elements(uexp)

; Make the final structure
;fdum = {id:'',field:'',ra:0.0d0,dec:0.0d0,ndet:0L,sepindx:lonarr(nind)-1,sepfindx:lonarr(nind)-1,$
fdum = {id:'',ra:0.0d0,dec:0.0d0,ndet:0L,sepindx:lonarr(nuexp)-1,sepfindx:lonarr(nuexp)-1,$
        u:99.99,uerr:9.99,g:99.99,gerr:9.99,r:99.99,rerr:9.99,i:99.99,ierr:9.99,z:99.99,zerr:9.99,ebv:99.99}
; SEPALL schema
sepalldum = {cmbindx:-1L,chipindx:-1L,id:-1L,x:0.0,y:0.0,mag:0.0,err:0.0,$
             chi:0.0,sharp:0.0,flag:-1,prob:-1.0,ra:0.0d0,dec:0.0d0}
; SEPINDX has NDET indices at the front of the array
; SEPFINDX has them in the element that matches the frame they were
; detected in
finaltags = tag_names(fdum)
lfinaltags = strlowcase(finaltags)

; This is the CHIPSTR schema
nan = !values.f_nan
dnan = !values.d_nan
chipstrdum = {field:'NAN',file:'NAN',expnum:'NAN',chip:-1L,base:'NAN',filter:'NAN',exptime:nan,utdate:'NAN',uttime:'NAN',$
              airmass:nan,gain:nan,rdnoise:nan,nx:-1L,ny:-1L,wcstype:'NAN',pixscale:nan,ra:dnan,dec:dnan,wcsrms:nan,fwhm:nan,$
              skymode:nan,skysig:nan,dao_nsources:-1L,dao_depth:nan,dao_npsfstars:-1L,dao_psftype:'NAN',dao_psfboxsize:-1L,$
              dao_psfvarorder:-1L,dao_psfchi:nan,alf_nsources:-1L,alf_depth:nan,calib_depth:nan,calib_color:'NAN',calib_zpterm:nan,$
              calib_amterm:nan,calib_colorterm:nan,calib_magname:'NAN',apcor:nan,ebv:nan}

; Loop through the exposures
fstr0 = fstr
undefine,fstr,fchstr,final,sepall
for i=0,nuexp-1 do begin
  expind = where(chstr.expnum eq uexp[i],nexpind)
  fstr0ind = where(fstr0.expnum eq uexp[i],nfstr0ind)
  push,fstr,fstr0[fstr0ind]

  ; Get all chips for this exposure
  undefine,expnew,expnewsepallindx
  for j=0,nexpind-1 do begin
    ; Add chip-level information
    ;  remove DATA, add sepallindx, to first element
    chtemp = chipstrdum
    add_tag,chtemp,'magoffset',nan,chtemp
    add_tag,chtemp,'ndata',-1L,chtemp
    add_tag,chtemp,'sepallindx',-1L,chtemp
    struct_assign,chstr[exind[j]],chtemp,/nozero
    chtemp.sepallindx = n_elements(sepall)

    ; Add catalog data for this chip
    temp = replicate(sepalldum,chstr[expind[j]].ndata)
    struct_assign,*chstr[expind[j]].data,temp,/nozero
    temp.chipindx = n_elements(fchstr)

    push,expnew,temp         ; add to exposure "new" structure
    push,fchstr,chtemp    ; add to FINAL chstr
    push,expnewsepallindx,lindgen(n_elements(temp))+n_elements(sepall)  ; sepall index for this exposure

    ; Concatenate TEMP structure to SEPALL
    ; concatenating two large structures causes lots to be zerod out
    if n_elements(sepall) gt 0 then begin
      nold = n_elements(sepall)
      nnew = n_elements(temp)
      new = replicate(sepalldum,nold+nnew)
      new[0:nold-1] = sepall
      new[nold:*] = temp
      sepall = new
      undefine,new
    endif else sepall=temp
  endfor ; chips in exposure loop


  ;------------------------------
  ; PUT IN FINAL MERGED CATALOG
  ;------------------------------

  ; Copy to new structure type
  If i eq 0 then begin
    final = replicate(fdum,n_elements(expnew))
    ;final.id = field+'_'+strtrim(all.chip,2)+'.'+strtrim(all.number,2)
    final.id = field+'.'+strtrim(lindgen(n_elements(expnew))+1,2)
    final.field = field
    final.ra = expnew.ra
    final.dec = expnew.dec

    ; Put SEPALL index in FINAL
    final.sepindx[0] = lindgen(n_elements(expnew))
    final.sepfindx[i] = lindgen(n_elements(expnew))
    final.ndet = 1
    ; Put ID, CMBINDX in SEPALL
    sepall.id = final.id
    sepall.cmbindx = lindgen(n_elements(final))

  ; 2nd and later exposures, check for repeats/overlap
  Endif else begin

    ; Match sources
    ;print,'Matching'
    dcr = 0.5  ; 0.7
    srcmatch,final.ra,final.dec,expnew.ra,expnew.dec,dcr,ind1,ind2,count=nmatch,/sph,domains=4000
    ; some simple tests suggest 4000 domains works best
    print,' ',strtrim(nmatch,2),' matched sources'
    ; Some matches, add data to existing record for these sources
    if nmatch gt 0 then begin
      ;sepallind = lindgen(n_elements(expnew))+n_elements(sepall)-n_elements(expnew)  ; index of all in sepall
      ;final[ind1].sepindx[final[ind1].ndet]= sepallind[ind2]     ; put SEPINDX
      for k=0LL,nmatch-1 do final[ind1[k]].sepindx[final[ind1[k]].ndet] = expnewsepallindx[ind2[k]]     ; put SEPINDX in FINAL
      final[ind1].sepfindx[j] = expnewsepallindx[ind2]    ; put SEPFINDX in FINAL
      final[ind1].ndet++
      sepall[expnewsepallindx[ind2]].id = final[ind1].id  ; put ID, CMBINDX in SEPALL
      sepall[expnewsepallindx[ind2]].cmbindx = ind1
      ; Remove stars
      if nmatch lt n_elements(expnew) then remove,ind2,expnew,expnewsepallindx else undefine,expnew,expnewsepallindx
    endif

    ; Some left, add records for these sources
    if n_elements(expnew) gt 0 then begin
      print,' ',strtrim(n_elements(expnew),2),' sources left to add'
      newfinal = replicate(fdum,n_elements(expnew))
      newfinal.id = field+'.'+strtrim(lindgen(n_elements(expnew))+1+n_elements(final),2)
      ;newfinal.id = field+'_'+strtrim(all.chip,2)+'.'+strtrim(all.number,2)  ;
      newfinal.field = field
      newfinal.ra = expnew.ra
      newfinal.dec = expnew.dec
      newfinal.sepindx[0] = expnewsepallindx      ; put SEPINDX in FINAL
      newfinal.sepfindx[j] = expnewsepallindx     ; put SEPFINDX in FINAL
      newfinal.ndet = 1
      sepall[expnewsepallindx].id = newfinal.id  ; put ID, CMBINDX in SEPALL
      sepall[expnewsepallindx].cmbindx = lindgen(n_elements(expnew))+n_elements(final)
      ; concatenating two large structures causes lots to be zerod out
      nold = n_elements(final)
      nnew = n_elements(newfinal)
      new = replicate(fdum,nold+nnew)
      new[0:nold-1] = final
      new[nold:*] = newfinal
      final = new
      undefine,new
      ;PUSH,final,newfinal         ; DOESN'T WORK THIS WAY
    endif

  Endelse  ; 2nd or later

bd = where(final.ra eq 0.0,nbd)
if nbd gt 0 then stop,'FINAL Zerod out elements problem!!!'

    ;stop

Endfor  ; exposure loop


; Combine photometry from same filter
for i=0,n_elements(ufilter)-1 do begin
  filtind = where(fstr.filter eq ufilter[i],nfiltind)
  ; only one exposure for this filter, copy
  if nfiltind eq 1 then begin
    magind = where(lfinaltags eq ufilter[i])
    ; get stars that have detections in this frame
    gd = where(final.sepfindx[filtind] ge 0,ngd)
    final[gd].(magind) = sepall[final[gd].sepfindx[filtind[0]]].mag
    errind = where(lfinaltags eq ufilter[i]+'err')
    final[gd].(errind) = sepall[final[gd].sepfindx[filtind[0]]].err
    bd = where(final.(magind) gt 50,nbd)
    if nbd gt 0 then begin
      final[bd].(magind) = 99.99
      final[bd].(errind) = 9.99
    endif

  ; multiple exposures for this filter to combine
  endif else begin

    mag = fltarr(n_elements(final),nfiltind)+99.99
    err = mag*0+9.99
    for k=0,nfiltind-1 do begin
      gd = where(final.sepfindx[filtind[k]] ge 0,ngd)
      mag[gd,k] = sepall[final[gd].sepfindx[filtind[k]]].mag
      err[gd,k] = sepall[final[gd].sepfindx[filtind[k]]].err
    endfor
    bd = where(mag eq 0.0 or mag gt 50,nbd)
    if nbd gt 0 then mag[bd]=!values.f_nan
    if nbd gt 0 then err[bd]=!values.f_nan

    ; copied from phot_overlap.pro
    flux = 2.511864d^mag
    wt = 1.0d0/err^2
    totalwt = total(wt,2,/nan)
    totalflux = total(flux*wt,2,/nan)
    totalerr = total((err^2)*wt,2,/nan) 
    newflux = totalflux/totalwt
    newmag = 2.50*alog10(newflux)
    newerr = sqrt(1.0/totalwt)

    bd = where(finite(newmag) eq 0,nbd)
    if nbd gt 0 then newmag[bd]=99.99
    if nbd gt 0 then newerr[bd]=9.99

    magind = where(lfinaltags eq ufilter[i])
    errind = where(lfinaltags eq ufilter[i]+'err')
    final.(magind) = newmag
    final.(errind) = newerr
  endelse  ; combine multiple exposures for this filter
endfor ; unique filter loop

print,'Getting SFD E(B-V)'
glactc,final.ra,final.dec,2000.0,lon,lat,1,/deg
ebv = dust_getval(lon,lat,/interp,/noloop)
; add EBV, G0 and I0 to the catalogs
final.ebv = ebv

; write final file out
outfile = outdir+field+'_combined'
print,'Writing combined file to ',outfile
mwrfits,fstr,outfile+'_exposures.fits',/create
mwrfits,fchstr,outfile+'_chips.fits',/create
mwrfits,sepall,outfile+'_sepall.fits',/create
mwrfits,final,outfile+'_final.fits',/create

;stop

end
