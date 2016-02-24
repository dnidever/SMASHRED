pro smashred_combine_ubercal,field,reduxdir=reduxdir,redo=redo

; Perform ubercal calibration on catalog with overlaps

if not keyword_set(reduxdir) then reduxdir='/data/smash/cp/red/photred/'
outdir = reduxdir+'catalogs/inst/comb/'

; Get reduction info
smashred_getredinfo,allinfo,/silent

gdinfo = where(allinfo.field eq field,ngdinfo)
if ngdinfo eq 0 then begin
  print,'No reduced exposures for ',field
  return
endif
info = allinfo[gdinfo]

;goto,starthere
;goto,combhere

; The various steps for this calibration:
; 1.) Load all of the chip data from .phot files
; 2.) Measure the photometric offsets + perform the Ubercal calibration
; 3.) Crossmatch all of the sources and build SEPALL and FINAL structures
; 4.) Combine photometry for all bands and get EBV

; Loop through the catalogs
print,'----------------------------------------------------------'
print,'--- STEP 1. Load all of the chip data from .phot files ---'
print,'=========================================================='
ninfo = n_elements(info)
undefine,allfstr,allchstr
For c=0,ninfo-1 do begin

  info1 = info[c]
  fstr = *info1.fstr
  fbase = file_basename(info1.file,'_summary.fits')  ; the observed field name
  print,strtrim(c+1,2),'/',strtrim(ninfo,2),'  ',field
  chstr = mrdfits(info1.file,2,/silent)  ; load the chip structure
  nchstr = n_elements(chstr)
  night = info1.night

  print,'Loading data for ',field

  ; Load in the data
  ;   don't have the "reference" frame information, just search
  ;   for PHOT files
  print,'---Loading the data for ',field,'---'
  outfile = reduxdir+'catalogs/inst/comb/'+field+'_'+night+'_ubercal.dat'  ; use REAL field name!!
  if file_test(outfile) eq 0 or keyword_set(redo) then begin
    add_tag,chstr,'refexpnum','',chstr
    add_tag,chstr,'vertices_ra',dblarr(4),chstr
    add_tag,chstr,'vertices_dec',dblarr(4),chstr
    add_tag,chstr,'ndata',-1L,chstr
    add_tag,chstr,'data',ptr_new(),chstr
    photfiles = file_search(reduxdir+night+'/'+chstr[0].field+'/'+chstr[0].field+'-*_??.phot',count=nphotfiles)
    print,strtrim(nphotfiles,2),' PHOT files for ',field,' in ',reduxdir+night+'/'+chstr[0].field+'/'
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
        mind = where(tags eq strtrim(chstr[chind[j]].calib_magname,2),nmind)
        gd = where(phot.(mind[0]) lt 50,ngd)
        temp = replicate({id:-1L,x:0.0,y:0.0,mag:0.0,err:0.0,chi:0.0,sharp:0.0,flag:-1,prob:-1.0,ra:0.0d0,dec:0.0d0},ngd)
        struct_assign,phot[gd],temp,/nozero
        temp.mag = phot[gd].(mind[0])
        temp.err = phot[gd].(mind[0]+1)   ; assume error is the next column

        chstr[chind[j]].refexpnum = refexpnum
        chstr[chind[j]].ndata = ngd
        chstr[chind[j]].data = ptr_new(temp)  ; save the data

        ; get astrometric vertices from header
        fitsfile = reduxdir+night+'/'+chstr[0].field+'/'+strtrim(chstr[chind[j]].base,2)+'.fits'
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

; Step 2. Measure the photometric offsets
print,'-----------------------------------------------------------------------------'
print,'--- STEP 2. Measure the photometric offsets + Perform Ubercal calibration ---'
print,'============================================================================='

; Add magnitude offset tag
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
  outfile_overlap = outdir+field+'_'+ifilter+'overlap.dat'  ; fits
  if file_test(outfile_overlap) eq 0 or keyword_set(redo) then begin

    SMASHRED_MEASURE_MAGOFFSET,chfiltstr,overlapstr

    ; Save the overlap structure
    print,'Writing overlaps to ',outfile_overlap
    save,overlapstr,file=outfile_overlap
    ;MWRFITS,overlapstr,outfile_overlap,/create
    ; saving it as FITS makes it 1D

  ; Using previously saved results
  Endif else begin
    print,'Using previously saved overlap file ',outfile_overlap
    restore,outfile_overlap
    ;overlapstr = mrdfits(outfile_overlap,1)
  Endelse

  ; Solve the UBERCAL problem
  PERFORM_UBERCAL_CALIB,overlapstr,fmagoff

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

starthere:
;restore,'/data/smash/cp/red/photred/catalogs/inst/comb/'+field+'_perform_ubercal_calib.dat'
;stop

; 3.) Crossmatch all of the sources and build SEPALL and FINAL structures
print,'-----------------------------------------------------------------------------------'
print,'--- STEP 3. Crossmatch all of the sources and build SEPALL and FINAL structures ---'
print,'==================================================================================='

; Unique exposures
uiexp = uniq(chstr.expnum,sort(chstr.expnum))
uexp = chstr[uiexp].expnum
nuexp = n_elements(uexp)

nan = !values.f_nan
dnan = !values.d_nan

; Make the final structure
;fdum = {id:'',field:'',ra:0.0d0,dec:0.0d0,ndet:0L,sepindx:lonarr(nind)-1,sepfindx:lonarr(nind)-1,$
fdum = {id:'',ra:0.0d0,dec:0.0d0,ndet:0L,sepindx:lonarr(nuexp)-1,sepfindx:lonarr(nuexp)-1,$
        u:99.99,uerr:9.99,g:99.99,gerr:9.99,r:99.99,rerr:9.99,i:99.99,ierr:9.99,z:99.99,zerr:9.99,chi:nan,sharp:nan,flag:-1,prob:nan,ebv:99.99}
; SEPALL schema
sepalldum = {cmbindx:-1L,chipindx:-1L,fid:'',id:-1L,x:0.0,y:0.0,mag:0.0,err:0.0,$
             chi:0.0,sharp:0.0,flag:-1,prob:-1.0,ra:0.0d0,dec:0.0d0}
sepall = replicate(sepalldum,5000000L)  ; initial sepall with 5 million elements
cur_sepall_indx = 0LL         ; next one starts from HERE
nsepall = n_elements(sepall)  ; current number of total SEPALl elements, NOT all filled
; SEPINDX has NDET indices at the front of the array
; SEPFINDX has them in the element that matches the frame they were
; detected in
finaltags = tag_names(fdum)
lfinaltags = strlowcase(finaltags)

; This is the CHIPSTR schema
chipstrdum = {field:'NAN',file:'NAN',expnum:'NAN',chip:-1L,base:'NAN',filter:'NAN',exptime:nan,utdate:'NAN',uttime:'NAN',$
              airmass:nan,gain:nan,rdnoise:nan,nx:-1L,ny:-1L,wcstype:'NAN',pixscale:nan,ra:dnan,dec:dnan,wcsrms:nan,fwhm:nan,$
              skymode:nan,skysig:nan,dao_nsources:-1L,dao_depth:nan,dao_npsfstars:-1L,dao_psftype:'NAN',dao_psfboxsize:-1L,$
              dao_psfvarorder:-1L,dao_psfchi:nan,alf_nsources:-1L,alf_depth:nan,calib_depth:nan,calib_color:'NAN',calib_zpterm:nan,$
              calib_amterm:nan,calib_colorterm:nan,calib_magname:'NAN',apcor:nan,ebv:nan}

; Loop through the exposures
fstr0 = fstr
undefine,fstr,fchstr,final
for i=0,nuexp-1 do begin
  expind = where(chstr.expnum eq uexp[i],nexpind)
  fstr0ind = where(fstr0.expnum eq uexp[i],nfstr0ind)
  push,fstr,fstr0[fstr0ind]

  print,strtrim(i+1,2),'/',strtrim(nuexp,2),' adding exposure ',uexp[i]

  ; Get all chips for this exposure
  undefine,expnew,expnewsepallindx
;dtall = 0
  for j=0,nexpind-1 do begin
    ; Add chip-level information
    ;  remove DATA, add sepallindx, to first element
    chtemp = chipstrdum
    add_tag,chtemp,'magoffset',nan,chtemp
    add_tag,chtemp,'ndata',-1L,chtemp
    add_tag,chtemp,'sepallindx',-1L,chtemp
    struct_assign,chstr[expind[j]],chtemp,/nozero
    chtemp.sepallindx = cur_sepall_indx

    ; Add catalog data for this chip
    temp = replicate(sepalldum,chstr[expind[j]].ndata)
    struct_assign,*chstr[expind[j]].data,temp,/nozero
    temp.chipindx = n_elements(fchstr)
    ntemp = n_elements(temp)
    temp_sepallindx = lindgen(ntemp)+cur_sepall_indx

    push,expnew,temp         ; add to exposure "new" structure
    push,fchstr,chtemp    ; add to FINAL chstr
    push,expnewsepallindx,temp_sepallindx  ; sepall index for this exposure

    ; New elements needed for SEPALL
    if max(temp_sepallindx) gt nsepall-1 then begin
t0 = systime(1)
      print,'Adding new elements to SEPALL'
      new = replicate(sepalldum,nsepall+5000000L)  ; add another 5 million elements
      new[0:nsepall-1] = sepall
      sepall = new
      undefine,new
      nsepall = n_elements(sepall)
print,'dt=',systime(1)-t0,' sec.  to add new SEPALL elements'
    endif

    ; Add TEMP structure to SEPALL
    sepall[temp_sepallindx] = temp
    cur_sepall_indx += ntemp
;    ; Concatenate TEMP structure to SEPALL
;    ; concatenating two large structures causes lots to be zerod out
;t0 = systime(1)
;    if n_elements(sepall) gt 0 then begin
;      nold = n_elements(sepall)
;      nnew = n_elements(temp)
;      new = replicate(sepalldum,nold+nnew)
;      new[0:nold-1] = sepall
;      new[nold:*] = temp
;      sepall = new
;      undefine,new
;    endif else sepall=temp
;dtall += systime(1)-t0
  endfor ; chips in exposure loo
;print,'dt=',dtall,' sec.  putting chips together'
  ;------------------------------
  ; PUT IN FINAL MERGED CATALOG
  ;------------------------------

  ; Copy to new structure type
  If i eq 0 then begin
    final = replicate(fdum,n_elements(expnew))
    ;final.id = field+'_'+strtrim(all.chip,2)+'.'+strtrim(all.number,2)
    final.id = field+'.'+strtrim(lindgen(n_elements(expnew))+1,2)
    ;final.field = field
    final.ra = expnew.ra
    final.dec = expnew.dec

    ; Put SEPALL index in FINAL
    final.sepindx[0] = lindgen(n_elements(expnew))
    final.sepfindx[i] = lindgen(n_elements(expnew))
    final.ndet = 1
    ; Put ID, CMBINDX in SEPALL
    sepall[0:cur_sepall_indx-1].fid = final.id
    sepall[0:cur_sepall_indx-1].cmbindx = lindgen(n_elements(final))

  ; 2nd and later exposures, check for repeats/overlap
  Endif else begin

    ; Match sources
    ;print,'Matching'
    dcr = 0.5  ; 0.7
    ;srcmatch,final.ra,final.dec,expnew.ra,expnew.dec,dcr,ind1,ind2,count=nmatch,/sph,domains=4000
t0 = systime(1)
    srcmatch,final.ra,final.dec,expnew.ra,expnew.dec,dcr,ind1,ind2,count=nmatch,/sph,/usehist  ; use faster histogram_nd method
print,'dt=',systime(1)-t0,' sec.  matching time'
    ; some simple tests suggest 4000 domains works best
    print,' ',strtrim(nmatch,2),' matched sources'
    ; Some matches, add data to existing record for these sources
    if nmatch gt 0 then begin
      ;sepallind = lindgen(n_elements(expnew))+n_elements(sepall)-n_elements(expnew)  ; index of all in sepall
      ;final[ind1].sepindx[final[ind1].ndet]= sepallind[ind2]     ; put SEPINDX
      for k=0LL,nmatch-1 do final[ind1[k]].sepindx[final[ind1[k]].ndet] = expnewsepallindx[ind2[k]]     ; put SEPINDX in FINAL
      final[ind1].sepfindx[i] = expnewsepallindx[ind2]    ; put SEPFINDX in FINAL
      final[ind1].ndet++
      sepall[expnewsepallindx[ind2]].fid = final[ind1].id  ; put ID, CMBINDX in SEPALL
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
      ;newfinal.field = field
      newfinal.ra = expnew.ra
      newfinal.dec = expnew.dec
      newfinal.sepindx[0] = expnewsepallindx      ; put SEPINDX in FINAL
      newfinal.sepfindx[i] = expnewsepallindx     ; put SEPFINDX in FINAL
      newfinal.ndet = 1
      sepall[expnewsepallindx].fid = newfinal.id  ; put ID, CMBINDX in SEPALL
      sepall[expnewsepallindx].cmbindx = lindgen(n_elements(expnew))+n_elements(final)
      ; concatenating two large structures causes lots to be zerod out
;t0 = systime(1)
      nold = n_elements(final)
      nnew = n_elements(newfinal)
      new = replicate(fdum,nold+nnew)
      new[0:nold-1] = final
      new[nold:*] = newfinal
      final = new
      undefine,new
;print,'dt=',systime(1)-t0
      ;PUSH,final,newfinal         ; DOESN'T WORK THIS WAY
    endif

  Endelse  ; 2nd or later

bd = where(final.ra eq 0.0,nbd)
if nbd gt 0 then stop,'FINAL Zerod out elements problem!!!'

    ;stop

Endfor  ; exposure loop

; Pruning extra SEPALL elements
if nsepall gt cur_sepall_indx then begin
  print,'Pruning extra SEPALl elements'
  sepall = sepall[0:cur_sepall_indx-1]
endif


; Step 4. Combine photometry for all bands and get EBV
print,'------------------------------------------------------------'
print,'--- STEP 4. Combine photometry for all bands and get EBV ---'
print,'============================================================'


combhere:
;outfile = outdir+field+'_combined'
;fstr = mrdfits(outfile+'_exposures.fits',1)
;chstr = mrdfits(outfile+'_chips.fits',1)
;sepall = mrdfits(outfile+'_sepall.fits',1)
;final = mrdfits(outfile+'_final.fits',1)
;add_tag,final,'chi',0.0,final
;add_tag,final,'sharp',0.0,final
;add_tag,final,'flag',0,final
;add_tag,final,'prob',0.0,final
;
;ui = uniq(fstr.filter,sort(fstr.filter))
;ufilter = fstr[ui].filter
;nufilter = n_elements(ufilter)
;
;finaltags = tag_names(final)
;lfinaltags = strlowcase(finaltags)

; Combine photometry from same filter
print,'Combing all of the photometry'
for i=0,n_elements(ufilter)-1 do begin
  filtind = where(fstr.filter eq ufilter[i],nfiltind)
  print,ufilter[i]

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
    chi = mag*0+!values.f_nan
    sharp = mag*0+!values.f_nan
    flag = long(mag)*0
    prob = mag*0+!values.f_nan
    for k=0,nfiltind-1 do begin
      gd = where(final.sepfindx[filtind[k]] ge 0,ngd)
      mag[gd,k] = sepall[final[gd].sepfindx[filtind[k]]].mag
      err[gd,k] = sepall[final[gd].sepfindx[filtind[k]]].err
      chi[gd,k] = sepall[final[gd].sepfindx[filtind[k]]].chi
      sharp[gd,k] = sepall[final[gd].sepfindx[filtind[k]]].sharp
      flag[gd,k] = sepall[final[gd].sepfindx[filtind[k]]].flag
      prob[gd,k] = sepall[final[gd].sepfindx[filtind[k]]].prob
    endfor
    ; non-ALLFRAME have prob=-1
    bd = where(finite(prob) eq 1 and prob lt -0.5,nbd)
    if nbd gt 0 then begin
      flag[bd] = 0
      prob[bd] = !values.f_nan
    endif
    ; Ignore non-detections
    bd = where(mag eq 0.0 or mag gt 50,nbd)
    if nbd gt 0 then begin
      mag[bd] = !values.f_nan
      err[bd] = !values.f_nan
      chi[bd] = !values.f_nan
      sharp[bd] = !values.f_nan
      flag[bd] = 0
      prob[bd] = !values.f_nan
    endif

    ; copied from phot_overlap.pro
    flux = 2.511864d^mag
    wt = 1.0d0/err^2
    totalwt = total(wt,2,/nan)
    totalflux = total(flux*wt,2,/nan)
    totalerr = total((err^2)*wt,2,/nan) 
    newflux = totalflux/totalwt
    newmag = 2.50*alog10(newflux)
    newerr = sqrt(1.0/totalwt)
    newchi = total(chi*wt,2,/nan)/totalwt
    newsharp = total(sharp*wt,2,/nan)/totalwt
    newprob = total(prob*wt,2,/nan)/totalwt
    ; newflag, do logical OR across all detections
    ;  0s are essentially ignored
    newflag = flag[*,0]*0
    for k=0,nfiltind-1 do newflag=newflag OR flag[*,k]
    bd = where(finite(newmag) eq 0,nbd)
    if nbd gt 0 then begin
      newmag[bd] = 99.99
      newerr[bd] = 9.99
      newchi[bd] = 99.99
      newsharp[bd] = 99.99
      newflag[bd] = -1
      newprob[bd] = 99.99
    endif

    magind = where(lfinaltags eq ufilter[i])
    errind = where(lfinaltags eq ufilter[i]+'err')
    final.(magind) = newmag
    final.(errind) = newerr
    final.chi = newchi
    final.sharp = newsharp
    final.flag = newflag
    final.prob = newprob
  endelse  ; combine multiple exposures for this filter
endfor ; unique filter loop

print,'Getting SFD E(B-V)'
glactc,final.ra,final.dec,2000.0,lon,lat,1,/deg
ebv = dust_getval(lon,lat,/interp,/noloop)
; add EBV, G0 and I0 to the catalogs
final.ebv = ebv

;stop

; write final file out
outfile = outdir+field+'_combined'
print,'Writing combined file to ',outfile
mwrfits,fstr,outfile+'_exposures.fits',/create
mwrfits,fchstr,outfile+'_chips.fits',/create
mwrfits,sepall,outfile+'_sepall.fits',/create
mwrfits,final,outfile+'_final.fits',/create

;stop

end
