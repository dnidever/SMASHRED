pro compare_decam09m_calibration_new_corr,save=save

; Compare DECam and 0.9m data for the new
; 0.9m runs

; Fields
;fieldnum = [104,113,48,51,71]
fieldnum = [109,133,134,135,142,149,150,171,180,84,94,98,99]
fields = 'Field'+strtrim(fieldnum,2)
nfields = n_elements(fields)

version = 'v2'

setdisp
filters = ['u','g','r','i','z']

;restore,'/data/smash/cp/red/photred/0.9m/new/compare_decam09m_calibration_new_corr.dat'
;goto,makeplots

offstr = replicate({field:'',type:'',nmatch:0L,offset:fltarr(5)+!values.f_nan,$
                    erroff:fltarr(5)+!values.f_nan,nused:lonarr(5)},nfields)

for i=0,nfields-1 do begin

  ifield = fields[i]
  ifieldnum = strmid(ifield,5)
  allobj = mrdfits('/data/smash/cp/red/photred/catalogs/final/'+version+'/'+ifield+'_combined_allobj.fits.gz',1,/silent)
  file2 = file_search('/data/smash/cp/red/photred/0.9m/new/SMASH_'+ifieldnum+'_*phot.fits',count=nfile2)
  ;file2 = file_search('/data/smash/cp/red/photred/0.9m/'+ifield+'_*phot.fits',count=nfile2)
  if nfile2 gt 1 then stop,'two '+ifield+' 0.9m files'
  base2 = file_basename(file2,'.fits')
  lo2 = strpos(base2,'_')
  hi2 = strpos(base2,'phot')
  type = strmid(base2,lo2+1,hi2-lo2-1)

  str = mrdfits(file2[0],1,/silent)
  srcmatch,allobj.ra,allobj.dec,str.ra,str.dec,1.0,ind1,ind2,/sph,count=nmatch

  ; Remove astrometric offset and refit
  raoff = median(allobj[ind1].ra-str[ind2].ra)
  decoff = median(allobj[ind1].dec-str[ind2].dec)
  sig = sqrt(mean( ((allobj[ind1].ra-str[ind2].ra-raoff)*3600.*cos(median(str.dec)/!radeg))^2 + (allobj[ind1].dec-str[ind2].dec-decoff)^2 ))
  print,'RA offset = ',raoff,'  DEC offset = ',decoff,'  SIG = ',sig
  ; use smaller matchin radius
  dcr = 0.4 > 2.5*sig < 1.0
  srcmatch,allobj.ra,allobj.dec,str.ra+raoff,str.dec+decoff,dcr,ind1,ind2,/sph,count=nmatch
  print,strtrim(nmatch,2),' matches'
  allobj2 = allobj[ind1]
  str2 = str[ind2]

  ;plotc,(allobj2.ra-str2.ra)*3600.*cos(median(str2.dec)/!radeg),(allobj2.dec-str2.dec)*3600.,allobj2.g,ps=1,$
  ;      xr=[-1.5,1.5],yr=[-1.5,1.5],xs=1,ys=1,tit=ifield
  ;oplot,[-10,10],[0,0],linestyle=2
  ;oplot,[0,0],[-10,10],linestyle=2
  ;oplot,[raoff*3600*cos(median(str2.dec)/!radeg)],[decoff*3600],ps=1,co=200,symsize=4.0

  ; convert to common allobj schema, remove indices
  allobj_schema = {id:'',field:0,ra:0.0d0,dec:0.0d0,u:0.0,uerr:0.0,g:0.0,gerr:0.0,r:0.0,rerr:0.0,i:0.0,ierr:0.0,$
                   z:0.0,zerr:0.0,chi:0.0,sharp:0.0,flag:0,prob:0.0,ebv:0.0}
  newobj = replicate(allobj_schema,nmatch)
  struct_assign,allobj2,newobj
  ifieldnum = long(strmid(ifield,5))
  newobj.field = ifieldnum

  ; Add info to str2
  add_tag,str2,'rad',0.0,str2
  add_tag,str2,'udiff',99.0,str2
  add_tag,str2,'udifferr',99.0,str2
  add_tag,str2,'gdiff',99.0,str2
  add_tag,str2,'gdifferr',99.0,str2
  add_tag,str2,'rdiff',99.0,str2
  add_tag,str2,'rdifferr',99.0,str2
  add_tag,str2,'idiff',99.0,str2
  add_tag,str2,'idifferr',99.0,str2
  add_tag,str2,'zdiff',99.0,str2
  add_tag,str2,'zdifferr',99.0,str2
  gd = where(newobj.u lt 50 and finite(str2.u) eq 1,ngd)
  if ngd gt 0 then begin
    str2[gd].udiff = str2[gd].u-newobj[gd].u
    str2[gd].udifferr = sqrt(str2[gd].uerr^2+newobj[gd].uerr^2)
  endif
  gd = where(newobj.g lt 50 and finite(str2.g) eq 1,ngd)
  if ngd gt 0 then begin
    str2[gd].gdiff = str2[gd].g-newobj[gd].g
    str2[gd].gdifferr = sqrt(str2[gd].gerr^2+newobj[gd].gerr^2)
  endif
  gd = where(newobj.r lt 50 and finite(str2.r) eq 1,ngd)
  if ngd gt 0 then begin
    str2[gd].rdiff = str2[gd].r-newobj[gd].r
    str2[gd].rdifferr = sqrt(str2[gd].rerr^2+newobj[gd].rerr^2)
  endif
  gd = where(newobj.i lt 50 and finite(str2.i) eq 1,ngd)
  if ngd gt 0 then begin
    str2[gd].idiff = str2[gd].i-newobj[gd].i
    str2[gd].idifferr = sqrt(str2[gd].ierr^2+newobj[gd].ierr^2)
  endif
  gd = where(newobj.z lt 50 and finite(str2.z) eq 1,ngd)
  if ngd gt 0 then begin
    str2[gd].zdiff = str2[gd].z-newobj[gd].z
    str2[gd].zdifferr = sqrt(str2[gd].zerr^2+newobj[gd].zerr^2)
  endif
  cenra = mean(minmax(str.ra))
  cendec = mean(minmax(str.dec))
  rotsphcen,str2.ra,str2.dec,cenra,cendec,lon,lat,/gnomic
  rad = sqrt(lon^2+lat^2)
  str2.rad = rad


  push,mallobj,newobj
  push,mstr,str2


  ; Check the offsets for this field
  mallobj1 = newobj
  mstr1 = str2
  otags = tag_names(mallobj1)
  stags = tag_names(mstr1)
  offstr[i].field = ifield
  offstr[i].type = type
  offstr[i].nmatch = nmatch
  for k=0,4 do begin

    if filters[k] eq 'u' then maglim=19.0 else maglim=18.0
    omagind = where(otags eq strupcase(filters[k]),nomagind)
    oerrind = where(otags eq strupcase(filters[k])+'ERR',noerrind)
    smagind = where(stags eq strupcase(filters[k]),nsmagind)
    serrind = where(stags eq strupcase(filters[k])+'ERR',nserrind)
    err = sqrt( (mallobj1.(oerrind))^2 + (mstr1.(serrind))^2 )
    ;gd = where(mallobj1.(omagind) lt maglim and finite(mstr1.(smagind)) eq 1 and finite(err) eq 1 and $
    ;           mallobj1.g lt 50 and mallobj1.i lt 50,ngd)
    ;if ngd lt 3 then goto,BOMB1
    ;mallobj2 = mallobj1[gd]
    ;mstr2 = mstr1[gd]
    ;err2 = err[gd]     
    ; Initial cut
    gd1 = where(mallobj1.(omagind) lt maglim and finite(mstr1.(smagind)) eq 1 and mallobj1.g lt 50 and mallobj1.i lt 50 and $
                mallobj1.g-mallobj1.i gt 0.5 and mallobj1.g-mallobj1.i lt 3.0,ngd1)
    if ngd1 lt 3 then goto,BOMB1
    med = median(mallobj1[gd1].(omagind)-mstr1[gd1].(smagind))
    sig = mad(mallobj1[gd1].(omagind)-mstr1[gd1].(smagind))
    ; Second cut to remove outliers
    gd2 = where(mallobj1.(omagind) lt maglim and finite(mstr1.(smagind)) eq 1 and mallobj1.g lt 50 and mallobj1.i lt 50 and $
                mallobj1.g-mallobj1.i gt 0.5 and mallobj1.g-mallobj1.i lt 3.0 and $
                err lt 0.1 and abs(mallobj1.(omagind)-mstr1.(smagind)-med) lt 3*sig,ngd2)
    coef0 = dln_poly_fit(mallobj1[gd2].g-mallobj1[gd2].i,mallobj1[gd2].(omagind)-mstr1[gd2].(smagind),0,$
                         measure_errors=err[gd2],sigma=sigma0,/bootstrap)
    coef1 = dln_poly_fit(mallobj1[gd2].g-mallobj1[gd2].i,mallobj1[gd2].(omagind)-mstr1[gd2].(smagind),1,$
                         measure_errors=err[gd2],sigma=sigma1,/bootstrap)
    ;print,'Diff=',strtrim(coef0[0],2),'+/-',strtrim(sigma0[0],2),' mag'
    offstr[i].offset[k] = coef0[0]
    offstr[i].erroff[k] = sigma0[0]
    offstr[i].nused[k] = ngd2
    BOMB1:
    ;stop
  endfor
  print,offstr[i].field+'  ',offstr[i].type,offstr[i].offset

  ;stop

endfor

; save the matched data
;save,mallobj,mstr,offstr,file='/data/smash/cp/red/photred/0.9m/new/compare_decam09m_calibration_new_corr.dat'

stop

MAKEPLOTS:

if keyword_set(save) then !p.font=0
setdisp
;maglim = 18.0


; Removing Field48 because it has problems, possibly because of crowding
;gd = where(mallobj.field ne 48,ngd)
;mallobj = mallobj[gd]
;mstr = mstr[gd]
;print,'REMOVING FIELD 48!!!!'


; Loop through the filters and compare the data
otags = tag_names(mallobj)
stags = tag_names(mstr)
undefine,pdffiles
for i=0,4 do begin
  print,strtrim(i+1,2),' ',filters[i]
  if filters[i] eq 'u' then maglim=19.0 else maglim=18.0
  omagind = where(otags eq strupcase(filters[i]),nomagind)
  oerrind = where(otags eq strupcase(filters[i])+'ERR',noerrind)
  smagind = where(stags eq strupcase(filters[i]),nsmagind)
  serrind = where(stags eq strupcase(filters[i])+'ERR',nserrind)
  err = sqrt( (mallobj.(oerrind))^2 + (mstr.(serrind))^2 )
  gd = where(mallobj.(omagind) lt maglim and finite(mstr.(smagind)) eq 1 and finite(err) eq 1 and $
             mallobj.g lt 50 and mallobj.i lt 50,ngd)
  mallobj2 = mallobj[gd]
  mstr2 = mstr[gd]
  err2 = err[gd]

  if keyword_set(save) then begin
    file = '/data/smash/cp/red/photred/0.9m/new/plots/compare_decam09m_calibration_new_'+filters[i]
    ;file = '/data/smash/cp/red/photred/0.9m/plots/compare_decam09m_calibration_new_'+filters[i]
    ps_open,file,/color,thick=4,/encap
  endif

  zmax = max(err2) < 0.2
  plotc,mallobj2.g-mallobj2.i,mallobj2.(omagind)-mstr2.(smagind),err2,ps=1,xr=[-1,3.5],yr=[-0.3,0.3],xtit='g-i',$
        ytit=filters[i]+' residuals (DECam-0.9m)',max=zmax,tit='Color coded by errors'
  oplot,[-5,5],[0,0],linestyle=2,co=250

  ; Initial cut
  gd1 = where(mallobj.(omagind) lt maglim and finite(mstr.(smagind)) eq 1 and mallobj.g lt 50 and mallobj.i lt 50 and $
              mallobj.g-mallobj.i gt 0.5 and mallobj.g-mallobj.i lt 3.0,ngd1)
  med = median(mallobj[gd1].(omagind)-mstr[gd1].(smagind))
  sig = mad(mallobj[gd1].(omagind)-mstr[gd1].(smagind))
  ; Second cut to remove outliers
  gd2 = where(mallobj.(omagind) lt maglim and finite(mstr.(smagind)) eq 1 and mallobj.g lt 50 and mallobj.i lt 50 and $
              mallobj.g-mallobj.i gt 0.5 and mallobj.g-mallobj.i lt 3.0 and $
              err lt 0.1 and abs(mallobj.(omagind)-mstr.(smagind)-med) lt 3*sig,ngd2)
  coef0 = dln_poly_fit(mallobj[gd2].g-mallobj[gd2].i,mallobj[gd2].(omagind)-mstr[gd2].(smagind),0,measure_errors=err[gd2],sigma=sigma0,/bootstrap)
  coef1 = dln_poly_fit(mallobj[gd2].g-mallobj[gd2].i,mallobj[gd2].(omagind)-mstr[gd2].(smagind),1,measure_errors=err[gd2],sigma=sigma1,/bootstrap)
  print,'Diff=',strtrim(coef0[0],2),'+/-',strtrim(sigma0[0],2),' mag'
  oplot,[0.5,2.0],[0,0]+coef0[0],co=0
  x = scale_vector(findgen(100),0.5,2.0)
  oplot,x,poly(x,coef1),co=250

  al_legend,['Median offset='+stringize(coef0[0],ndec=3)+' mag'],/top,/right,charsize=1.3

  if keyword_set(save) then begin
    ps_close
    ps2png,file+'.eps',/eps
    spawn,['epstopdf',file+'.eps'],/noshell
    push,pdffiles,file+'.pdf'
  endif

  ; Plot residuals vs. mag, color coded by field.
  if keyword_set(save) then begin
    file = '/data/smash/cp/red/photred/0.9m/plots/compare_decam09m_calibration_new_'+filters[i]+'_mag'
    ps_open,file,/color,thick=4,/encap
  endif
  gd2 = where(mallobj.(omagind) lt 50 and finite(mstr.(smagind)) eq 1 and mallobj.g-mallobj.i gt 0.5 and mallobj.g-mallobj.i lt 3.0,ngd2)
  xr = minmax(mallobj[gd2].(omagind))
  xr = [xr[0], xr[1] < 21.0]
  plotc,mallobj[gd2].(omagind),mallobj[gd2].(omagind)-mstr[gd2].(smagind),mallobj[gd2].field,ps=1,xs=1,xr=xr,yr=[-0.3,0.3],xtit=filters[i],$
        ytit=filters[i]+' residuals (DECam-0.9m)',tit='Color coded by Field number'
  oplot,[0,30],[0,0],linestyle=2
  al_legend,['Median offset='+stringize(coef0[0],ndec=3)+' mag'],/top,/right,charsize=1.3
  if keyword_set(save) then begin
    ps_close
    ps2png,file+'.eps',/eps
    spawn,['epstopdf',file+'.eps'],/noshell
    push,pdffiles,file+'.pdf'
  endif

  stop
endfor

; g-band is off by ~0.1 mag but the rest seem to be within ~0.01-0.02
; but some have smallish color terms.

; Field48 has large scatter in g-band.  What's going on there.
; astrometric problem.

if keyword_set(save) then begin
  spawn,'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=/data/smash/cp/red/photred/0.9m/new/plots/compare_decam09m_calibration_new.pdf '+strjoin(pdffiles,' ')
endif

stop

end
