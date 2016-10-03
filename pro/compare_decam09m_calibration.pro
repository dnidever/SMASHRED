pro compare_decam09m_calibration

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
  print,strtrim(nmatch,2),' matches'
  allobj2 = allobj[ind1]
  str2 = str[ind2]

  ; convert to common allobj schema, remove indices
  allobj_schema = {id:'',ra:0.0d0,dec:0.0d0,u:0.0,uerr:0.0,g:0.0,gerr:0.0,r:0.0,rerr:0.0,i:0.0,ierr:0.0,$
                   z:0.0,zerr:0.0,chi:0.0,sharp:0.0,flag:0,prob:0.0,ebv:0.0}
  newobj = replicate(allobj_schema,nmatch)
  struct_assign,allobj2,newobj
  push,mallobj,newobj
  push,mstr,str2

endfor

; save the matched data
;save,mallobj,mstr,file='/data/smash/cp/red/photred/0.9m/compare_decam09m_calibration.dat'

stop

MAKEPLOTS:

; Loop through the filters and compare the data
filters = ['u','g','r','i','z']
otags = tag_names(mallobj)
stags = tag_names(mstr)
for i=0,4 do begin
  print,strtrim(i+1,2),' ',filters[i]
  omagind = where(otags eq strupcase(filters[i]),nomagind)
  oerrind = where(otags eq strupcase(filters[i])+'ERR',noerrind)
  smagind = where(stags eq strupcase(filters[i]),nsmagind)
  serrind = where(stags eq strupcase(filters[i])+'ERR',nserrind)
  err = sqrt( (mallobj.(oerrind))^2 + (mstr.(serrind))^2 )
  gd = where(mallobj.(omagind) lt 50 and finite(mstr.(smagind)) eq 1 and finite(err) eq 1,ngd)
  mallobj2 = mallobj[gd]
  mstr2 = mstr[gd]
  err2 = err[gd]
  plotc,mallobj2.g-mallobj2.i,mallobj2.(omagind)-mstr2.(smagind),err2,ps=1,xr=[-1,3.5],yr=[-0.3,0.3],xtit='g-i',ytit=filters[i],max=0.2
  oplot,[-5,5],[0,0],linestyle=2,co=250

  gd1 = where(mallobj.(omagind) lt 50 and finite(mstr.(smagind)) eq 1 and mallobj.g-mallobj.i gt 0.5 and mallobj.g-mallobj.i lt 3.0 and $
              err lt 0.1,ngd1)
  err = sqrt( (mallobj[gd1].(oerrind))^2 + (mstr[gd1].(serrind))^2 )
  med = median(mallobj[gd1].(omagind)-mstr[gd1].(smagind))
  sig = mad(mallobj[gd1].(omagind)-mstr[gd1].(smagind))
  gd2 = where(mallobj.(omagind) lt 50 and finite(mstr.(smagind)) eq 1 and mallobj.g-mallobj.i gt 0.5 and mallobj.g-mallobj.i lt 3.0 and $
              err lt 0.1 and abs(mallobj.(omagind)-mstr.(smagind)-med) lt 3*sig,ngd2)
  coef0 = dln_poly_fit(mallobj[gd2].g-mallobj[gd2].i,mallobj[gd2].(omagind)-mstr[gd2].(smagind),0,measure_errors=err[gd2],sigma=sigma0,/bootstrap)
  coef1 = dln_poly_fit(mallobj[gd2].g-mallobj[gd2].i,mallobj[gd2].(omagind)-mstr[gd2].(smagind),1,measure_errors=err[gd2],sigma=sigma1,/bootstrap)
  print,'Diff=',strtrim(coef0[0],2),'+/-',strtrim(sigma0[0],2),' mag'
  oplot,[0.5,2.0],[0,0]+coef0[0],co=80
  x = scale_vector(findgen(100),0.5,2.0)
  oplot,x,poly(x,coef1),co=150

  stop
endfor

; g-band is off by ~0.1 mag but the rest seem to be within ~0.01-0.02
; but some have smallish color terms.

stop

end
