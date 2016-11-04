pro correct_new09mphot

; Compare DECam and 0.9m data for the new
; 0.9m runs and CORRECT any offsets.

; Fields
;fieldnum = [109,133,134,135,142,149,150,171,180,84,94,98,99]
fieldnum = [109,133,134,135,149,150,171,180,84,94,98,99]
fields = 'Field'+strtrim(fieldnum,2)
nfields = n_elements(fields)

version = 'v2'

setdisp
filters = ['u','g','r','i','z']

restore,'/data/smash/cp/red/photred/0.9m/correct_new09mphot.dat'
goto,fitting

offstr = replicate({field:'',type:'',nmatch:0L,offset:fltarr(5)+!values.f_nan,$
                    erroff:fltarr(5)+!values.f_nan,nused:lonarr(5)},nfields)

for i=0,nfields-1 do begin
  ifield = fields[i]
  allobj = mrdfits('/data/smash/cp/red/photred/catalogs/final/'+version+'/'+ifield+'_combined_allobj.fits.gz',1,/silent)
  file2 = file_search('/data/smash/cp/red/photred/0.9m/'+ifield+'_*phot.fits',count=nfile2)
  if nfile2 gt 1 then stop,'two '+ifield+' 0.9m files'
  base2 = file_basename(file2,'.fits')
  lo2 = strpos(base2,'_')
  hi2 = strpos(base2,'phot')
  type = strmid(base2,lo2+1,hi2-lo2-1)

  ; TEST the corrected files to make sure they are now good
  ;file2 = file_search('/data/smash/cp/red/photred/0.9m/SMASH_'+strtrim(fieldnum[i],2)+'_*phot_corr.fits',count=nfile2)
  ;base2 = file_basename(file2,'.fits')
  ;dum = strsplit(base2,'_',/extract)
  ;hi2 = strpos(dum[2],'phot')
  ;type = strmid(dum[2],0,hi2)

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

;save,mallobj,mstr,offstr,file='/data/smash/cp/red/photred/0.9m/correct_new09mphot.dat'
stop

fitting:

; DES u-band bad on 171 and 94
bd = where((mstr.field eq 'SMASH_171' or mstr.field eq 'SMASH_94') and $
            mstr.system eq 'DES',nbd)
mstr[bd].u = !values.f_nan
mstr[bd].udiff = 99.99

; All bands show the same spatial/radial patterns
; fit all the points together
stags = tag_names(mstr)
undefine,rad,diff,mag
for i=0,4 do begin
  magind = where(stags eq strupcase(filters[i]),nmagind)
  diffind = where(stags eq strupcase(filters[i]+'DIFF'),ndiffind)
  push,rad,mstr.rad
  push,diff,mstr.(diffind)-median(mstr.(diffind))
  push,mag,mstr.(magind)
endfor
gd = where(finite(mag) eq 1 and mag lt 18,ngd,comp=bd,ncomp=nbd)
remove,bd,rad,diff,mag
bindata,rad,diff,xbin,ybin,binsize=0.01,min=0,max=0.16,/med,gdind=gdind
coef = robust_poly_fitq(xbin[gdind],ybin[gdind],3)
plotc,rad,diff,mag,ps=3,yr=[-0.4,0.4]
oplot,rad,poly(rad,coef),ps=8,co=150

; Now derive corrections for each band and system
magstr = replicate({filter:'',system:'',coef:fltarr(4)},10)
systems = ['sdss','DES']
for i=0,1 do begin
  isystem = systems[i]
  ind = where(mstr.system eq isystem,nind)
  mstr1 = mstr[ind]
  for j=0,4 do begin
    magind = where(stags eq strupcase(filters[j]),nmagind)
    diffind = where(stags eq strupcase(filters[j]+'DIFF'),ndiffind)

    gd = where(finite(mstr1.(magind)) eq 1 and mstr1.(magind) lt 18 and $
               abs(mstr1.(diffind)) lt 1 and mstr1.g-mstr1.i gt 0.5 and mstr1.g-mstr1.i lt 3.0,ngd)
    med = median(mstr1[gd].(diffind)-poly(mstr1[gd].rad,coef))
    coef2 = coef
    coef2[0] += med
    print,isystem,'  ',filters[j],'  ',med

    plotc,mstr1[gd].rad,mstr1[gd].(diffind),mstr1[gd].(magind),ps=3,yr=[-0.4,0.4],tit=isystem+' '+filters[j]
    oplot,mstr1[gd].rad,poly(mstr1[gd].rad,coef2),ps=8,co=250

    ; median and scatter in final residuals
    resid = mstr1[gd].(diffind)-poly(mstr1[gd].rad,coef2)
    medresid = median(resid)
    sigresid = mad(resid-medresid,/zero)
    print,'Median resid = ',stringize(medresid,ndec=4),'  scatter=',stringize(sigresid,ndec=4)

    magstr[i*5+j].system = isystem
    magstr[i*5+j].filter = filters[j]
    magstr[i*5+j].coef = coef2
stop
  endfor
endfor

stop

; Correct the new 0.9m photometry catalogs
;-----------------------------------------
print,'Correct the new 0.9m photometry catalogs'
dir = '/data/smash/cp/red/photred/0.9m/'
files = file_search(dir+'SMASH_*phot.fits',count=nfiles)
for i=0,nfiles-1 do begin
  phot = MRDFITS(files[i],1,/silent)
  phot0 = phot
  tags = tag_names(phot)

  base = file_basename(files[i],'.fits')
  dum = strsplit(base,'_',/extract)
  hi = strpos(dum[2],'phot')
  type = strmid(dum[2],0,hi)
  fieldnum = strtrim(dum[1],2)
  print,strtrim(i+1,2),' ',fieldnum,' ',type

  cenra = mean(minmax(phot.ra))
  cendec = mean(minmax(phot.dec))
  rotsphcen,phot.ra,phot.dec,cenra,cendec,lon,lat,/gnomic
  rad = sqrt(lon^2+lat^2)

  for j=0,4 do begin
    magind = where(tags eq strupcase(filters[j]),nmagind)
    ind = where(strupcase(magstr.system) eq strupcase(type) and magstr.filter eq filters[j],nind)
    corr = poly(rad,magstr[ind[0]].coef)
    gd = where(finite(phot.(magind)) eq 1,ngd)
    if ngd gt 0 then phot[gd].(magind) -= corr[gd]
  endfor

  newfile = dir+base+'_corr.fits'
  print,'Writing corrected file to ',newfile
  MWRFITS,phot,newfile,/create
  ;stop
endfor

stop

end
