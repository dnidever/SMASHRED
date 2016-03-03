;+
;
; SMASHRED_CROSSMATCH
;
; This program crossmatches sources for a field and creates
; the ALLSTR and ALLOBJ structures.
;
; INPUTS:
;  fstr    The structure with information for each exposure.
;  chstr   The structure with information for each chip.
;  =dcr    The matching radius in arcsec.  The default is 0.5 arcsec.
;  /silent Don't print anything to the screen.
;
; OUTPUTS:
;  allstr  The structure with information for each source detection.
;  allobj  The structure with information for each unique object.
;  =error  The error message if one occurred.
;
; USAGE:
;  IDL>smashred_crossmatch,fstr,chstr,allstr,allobj
; 
; By D.Nidever  March 2016
;-

pro smashred_crossmatch,fstr,chstr,allstr,allobj,dcr=dcr,error=error,silent=silent

undefine,allstr,allobj

; Not enough inputs
if n_elements(fstr) eq 0 or n_elements(chstr) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_crossmatch,fstr,chstr,allstr,allobj'
  return
endif

nfstr = n_elements(fstr)
nchstr = n_elements(chstr)

; Defaults
if n_elements(dcr) eq 0 then dcr=0.5  ; 0.7

; Unique exposures from CHSTR
uiexp = uniq(chstr.expnum,sort(chstr.expnum))
uexp = chstr[uiexp].expnum
nuexp = n_elements(uexp)
; Unique exposures from FSTR
uifexp = uniq(fstr.expnum,sort(fstr.expnum))
ufexp = fstr[uifexp].expnum
nufexp = n_elements(ufexp)
; unqiue FSTR and CHSTR exposures don't match
if nuexp ne nfstr then begin
  error = 'Exposures do NOT match in FSTR and CHSTR'
  if not keyword_set(silent) then print,error
  return
endif
; duplicate exposures in FSTR
if nufexp ne nfstr then begin
  error = 'Duplicate exposures in FSTR'
  if not keyword_set(silent) then print,error
  return
endif


nan = !values.f_nan
dnan = !values.d_nan

; Make the ALLOBJ structure schema
fdum = {id:'',ra:0.0d0,dec:0.0d0,ndet:0L,depthflag:0B,srcindx:lonarr(nuexp)-1,srcfindx:lonarr(nuexp)-1,$
        u:99.99,uerr:9.99,g:99.99,gerr:9.99,r:99.99,rerr:9.99,i:99.99,ierr:9.99,z:99.99,zerr:9.99,chi:nan,sharp:nan,flag:-1,prob:nan,ebv:99.99}
; ALLSRC schema
allsrcdum = {cmbindx:-1L,chipindx:-1L,fid:'',id:-1L,x:0.0,y:0.0,mag:0.0,err:0.0,$
             chi:0.0,sharp:0.0,flag:-1,prob:-1.0,ra:0.0d0,dec:0.0d0}
allsrc = replicate(allsrcdum,5000000L)  ; initial allsrc with 5 million elements
cur_allsrc_indx = 0LL         ; next one starts from HERE
nallsrc = n_elements(allsrc)  ; current number of total Allsrc elements, NOT all filled
; SRCINDX has NDET indices at the front of the array
; SRCFINDX has them in the element that matches the frame they were
; detected in
allobjtags = tag_names(fdum)
lallobjtags = strlowcase(allobjtags)

; This is the CHIPSTR schema
;chipstrdum = {field:'NAN',file:'NAN',expnum:'NAN',chip:-1L,base:'NAN',filter:'NAN',exptime:nan,utdate:'NAN',uttime:'NAN',$
;              airmass:nan,gain:nan,rdnoise:nan,nx:-1L,ny:-1L,wcstype:'NAN',pixscale:nan,ra:dnan,dec:dnan,wcsrms:nan,fwhm:nan,$
;              skymode:nan,skysig:nan,dao_nsources:-1L,dao_depth:nan,dao_npsfstars:-1L,dao_psftype:'NAN',dao_psfboxsize:-1L,$
;              dao_psfvarorder:-1L,dao_psfchi:nan,alf_nsources:-1L,alf_depth:nan,calib_depth:nan,calib_color:'NAN',calib_zpterm:nan,$
;              calib_amterm:nan,calib_colorterm:nan,calib_magname:'NAN',apcor:nan,ebv:nan}
add_tag,chstr,'allsrcindx',-1L,chstr

; Loop through the exposures
;fstr0 = fstr
;undefine,fstr,fchstr,allstr,allobj
;for i=0,nuexp-1 do begin
for i=0,nfstr-1 do begin
  ;expind = where(chstr.expnum eq uexp[i],nexpind)
  ;fstr0ind = where(fstr0.expnum eq uexp[i],nfstr0ind)
  ;push,fstr,fstr0[fstr0ind]
  expind = where(chstr.expnum eq fstr[i].exposure,nexpind)
  if fstr[i].exptime lt 100 then depthbit=1 else depthbit=2  ; short or long

  print,strtrim(i+1,2),'/',strtrim(nfstr,2),' adding exposure ',fstr[i].exposure

  ; Get all chips for this exposure
  undefine,expnew,expnewallsrcindx
;dtall = 0
  for j=0,nexpind-1 do begin
    ;; Add chip-level information
    ;;  remove DATA, add allsrcindx, index to first element in ALLSRC
    ;chtemp = chipstrdum
    ;add_tag,chtemp,'magoffset',nan,chtemp
    ;add_tag,chtemp,'ndata',-1L,chtemp
    ;add_tag,chtemp,'allsrcindx',-1L,chtemp
    ;struct_assign,chstr[expind[j]],chtemp,/nozero
    ;chtemp.allsrcindx = cur_allsrc_indx
    ; DON'T REMOVE DATA HERE!!! Can do that later
    chstr[expind[j]].allsrcindx = cur_allsrc_indx

    ; Add catalog data for this chip
    temp = replicate(allsrcdum,chstr[expind[j]].ndata)
    struct_assign,*chstr[expind[j]].data,temp,/nozero
    ;temp.chipindx = n_elements(fchstr)
    temp.chipindx = expind[j]
    ntemp = n_elements(temp)
    temp_allsrcindx = lindgen(ntemp)+cur_allsrc_indx

    push,expnew,temp         ; add to exposure "new" structure
    ;push,fchstr,chtemp    ; add to ALLOBJ chstr
    push,expnewallsrcindx,temp_allsrcindx  ; allsrc index for this exposure

    ; New elements needed for ALLSRC
    if max(temp_allsrcindx) gt nallsrc-1 then begin
t0 = systime(1)
      print,'Adding new elements to ALLSRC'
      new = replicate(allsrcdum,nallsrc+5000000L)  ; add another 5 million elements
      new[0:nallsrc-1] = allsrc
      allsrc = new
      undefine,new
      nallsrc = n_elements(allsrc)
print,'dt=',systime(1)-t0,' sec.  to add new ALLSRC elements'
    endif

    ; Add TEMP structure to ALLSRC
    allsrc[temp_allsrcindx] = temp
    cur_allsrc_indx += ntemp
  endfor ; chips in exposure loo

  ;------------------------------
  ; PUT IN ALLOBJ MERGED CATALOG
  ;------------------------------
  ; Copy to new structure type
  If i eq 0 then begin
    allobj = replicate(fdum,n_elements(expnew))
    ;allobj.id = field+'_'+strtrim(all.chip,2)+'.'+strtrim(all.number,2)
    allobj.id = field+'.'+strtrim(lindgen(n_elements(expnew))+1,2)
    ;allobj.field = field
    allobj.ra = expnew.ra
    allobj.dec = expnew.dec

    ; Put ALLSRC index in ALLOBJ
    allobj.srcindx[0] = lindgen(n_elements(expnew))
    allobj.srcfindx[i] = lindgen(n_elements(expnew))
    allobj.ndet = 1
    allobj.depthflag OR= depthbit             ; OR combine to depthflag, 1-short, 2-long, 3-short+long
    ; Put ID, CMBINDX in ALLSRC
    allsrc[0:cur_allsrc_indx-1].fid = allobj.id
    allsrc[0:cur_allsrc_indx-1].cmbindx = lindgen(n_elements(allobj))

  ; 2nd and later exposures, check for repeats/overlap
  Endif else begin

    ; Match sources
t0 = systime(1)
    srcmatch,allobj.ra,allobj.dec,expnew.ra,expnew.dec,dcr,ind1,ind2,count=nmatch,/sph,/usehist  ; use faster histogram_nd method
print,'dt=',systime(1)-t0,' sec.  matching time'
    print,' ',strtrim(nmatch,2),' matched sources'
    ; Some matches, add data to existing record for these sources
    if nmatch gt 0 then begin
      for k=0LL,nmatch-1 do allobj[ind1[k]].srcindx[allobj[ind1[k]].ndet] = expnewallsrcindx[ind2[k]]     ; put SRCINDX in ALLOBJ
      allobj[ind1].srcfindx[i] = expnewallsrcindx[ind2]    ; put SRCFINDX in ALLOBJ
      allobj[ind1].ndet++
      allobj[ind1].depthflag OR= depthbit                  ; OR combine to depthflag  
      allsrc[expnewallsrcindx[ind2]].fid = allobj[ind1].id  ; put ID, CMBINDX in ALLSRC
      allsrc[expnewallsrcindx[ind2]].cmbindx = ind1
      ; Remove stars
      if nmatch lt n_elements(expnew) then remove,ind2,expnew,expnewallsrcindx else undefine,expnew,expnewallsrcindx
    endif

    ; Some left, add records for these sources
    if n_elements(expnew) gt 0 then begin
      print,' ',strtrim(n_elements(expnew),2),' sources left to add'
      newallobj = replicate(fdum,n_elements(expnew))
      newallobj.id = field+'.'+strtrim(lindgen(n_elements(expnew))+1+n_elements(allobj),2)
      ;newallobj.id = field+'_'+strtrim(all.chip,2)+'.'+strtrim(all.number,2)  ;
      ;newallobj.field = field
      newallobj.ra = expnew.ra
      newallobj.dec = expnew.dec
      newallobj.srcindx[0] = expnewallsrcindx      ; put SRCINDX in ALLOBJ
      newallobj.srcfindx[i] = expnewallsrcindx     ; put SRCFINDX in ALLOBJ
      newallobj.ndet = 1
      newallobj.depthflag OR= depthbit             ; OR combine to depthflag
      allsrc[expnewallsrcindx].fid = newallobj.id  ; put ID, CMBINDX in ALLSRC
      allsrc[expnewallsrcindx].cmbindx = lindgen(n_elements(expnew))+n_elements(allobj)
      ; concatenating two large structures causes lots to be zerod out
;t0 = systime(1)
      nold = n_elements(allobj)
      nnew = n_elements(newallobj)
      new = replicate(fdum,nold+nnew)
      new[0:nold-1] = allobj
      new[nold:*] = newallobj
      allobj = new
      undefine,new
;print,'dt=',systime(1)-t0
    endif

  Endelse  ; 2nd or later

bd = where(allobj.ra eq 0.0,nbd)
if nbd gt 0 then stop,'ALLOBJ Zerod out elements problem!!!'

    ;stop

Endfor  ; exposure loop

; Pruning extra ALLSRC elements
if nallsrc gt cur_allsrc_indx then begin
  print,'Pruning extra Allsrc elements'
  allsrc = allsrc[0:cur_allsrc_indx-1]
endif

;stop

end
