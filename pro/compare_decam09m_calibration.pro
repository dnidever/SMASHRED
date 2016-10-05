pro compare_decam09m_calibration,save=save

; Fields
fieldnum = [104,113,48,51,71]
fields = 'Field'+strtrim(fieldnum,2)
nfields = n_elements(fields)

restore,'/data/smash/cp/red/photred/0.9m/compare_decam09m_calibration.dat'
goto,makeplots

for i=0,nfields-1 do begin

  ifield = fields[i]
  allobj = mrdfits('/data/smash/cp/red/photred/catalogs/final/'+ifield+'_combined_allobj.fits.gz',1)
  str = mrdfits('/data/smash/cp/red/photred/0.9m/'+ifield+'_phot.fits',1)
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

  plotc,(allobj2.ra-str2.ra)*3600.*cos(median(str2.dec)/!radeg),(allobj2.dec-str2.dec)*3600.,allobj2.g,ps=1,$
        xr=[-1.5,1.5],yr=[-1.5,1.5],xs=1,ys=1,tit=ifield
  oplot,[-10,10],[0,0],linestyle=2
  oplot,[0,0],[-10,10],linestyle=2
  oplot,[raoff*3600*cos(median(str2.dec)/!radeg)],[decoff*3600],ps=1,co=200,symsize=4.0


  ; convert to common allobj schema, remove indices
  allobj_schema = {id:'',field:0,ra:0.0d0,dec:0.0d0,u:0.0,uerr:0.0,g:0.0,gerr:0.0,r:0.0,rerr:0.0,i:0.0,ierr:0.0,$
                   z:0.0,zerr:0.0,chi:0.0,sharp:0.0,flag:0,prob:0.0,ebv:0.0}
  newobj = replicate(allobj_schema,nmatch)
  struct_assign,allobj2,newobj
  ifieldnum = long(strmid(ifield,5))
  newobj.field = ifieldnum

  push,mallobj,newobj
  push,mstr,str2

  ;stop

endfor

; save the matched data
;save,mallobj,mstr,file='/data/smash/cp/red/photred/0.9m/compare_decam09m_calibration.dat'

stop

MAKEPLOTS:

if keyword_set(save) then !p.font=0
setdisp
;maglim = 18.0


; Removing Field48 because it has problems, possibly because of crowding
gd = where(mallobj.field ne 48,ngd)
mallobj = mallobj[gd]
mstr = mstr[gd]
print,'REMOVING FIELD 48!!!!'


; Loop through the filters and compare the data
filters = ['u','g','r','i','z']
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
    file = '/data/smash/cp/red/photred/0.9m/plots/compare_decam09m_calibration_'+filters[i]
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
    file = '/data/smash/cp/red/photred/0.9m/plots/compare_decam09m_calibration_'+filters[i]+'_mag'
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

  ;stop
endfor

; g-band is off by ~0.1 mag but the rest seem to be within ~0.01-0.02
; but some have smallish color terms.

; Field48 has large scatter in g-band.  What's going on there.
; astrometric problem.

if keyword_set(save) then begin
  spawn,'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=plots/compare_decam09m_calibration.pdf '+strjoin(pdffiles,' ')
endif

stop

end
