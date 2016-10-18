pro qacheck_spatialtransvar,field,redo=redo,stp=stp

; This checks for spatial variations in the
; transparency of exposures due to patchy clouds

if n_elements(field) eq 0 then begin
  print,'Syntax - qacheck_spatialtransvar,field,redo=redo'
  return
endif

setdisp,/silent

reduxdir = '/data/smash/cp/red/photred/'
gaiadir = reduxdir+'gaia/'
catdir = reduxdir+'catalogs/final/'
version = 'v2'
if n_elements(version) gt 0 then catdir+=version+'/'
plotsdir = catdir+'plots/'
if file_test(plotsdir,/directory) eq 0 then file_mkdir,plotsdir

print,'Running QACHECK_SPATIALTRANSVAR on field ',field

outcatfile = plotsdir+field+'_spatialtransvar.fits'
if file_test(outcatfile) eq 1 and not keyword_set(redo) then begin
  print,outcatfile,' already EXISTS and /redo not set'
  return
endif

; Restoring the files
allobj = mrdfits(catdir+field+'_combined_allobj_bright.fits',1)
allobj.id = strtrim(allobj.id,2)
allsrc = mrdfits(catdir+field+'_combined_allsrc.fits.gz',1)
chstr = mrdfits(catdir+field+'_combined_chips.fits.gz',1)
gaia = mrdfits(gaiadir+field+'_gaia.fits',1)
; need to restore the photred final file

otags = tag_names(allobj)

; Match up ALLOBJ and GAIA sources
srcmatch,allobj.ra,allobj.dec,gaia.ra_icrs,gaia.de_icrs,0.5,ind1,ind2,/sph,count=nmatch
print,strtrim(nmatch,2),' between GAIA and bright ALLOBJ'
mallobj = allobj[ind1]
mgaia = gaia[ind2]

; Exposure loop
ui = uniq(chstr.expnum,sort(chstr.expnum))
uexpnum = chstr[ui].expnum
nexpnum = n_elements(uexpnum)
expstr = replicate({expnum:'',filter:'',exptime:0.0,totexpsrc:0L,ngaiamatch:0L,nused:0L,med:0.0,sig:0.0},nexpnum)
undefine,psfiles
for i=0,nexpnum-1 do begin

  ; Get all of the chips and sources for this exposure
  chind = where(chstr.expnum eq uexpnum[i],nchind)
  exptotsrc = total(chstr[chind].nsrc)
  ; Get blank allsrc schema
  schema = allsrc[0]
  struct_assign,{dum:''},schema
  expsrc = replicate(schema,exptotsrc)
  cnt = 0LL
  for j=0,nchind-1 do begin
    chstr1 = chstr[chind[j]]
    srcind = lindgen(chstr1.nsrc)+chstr1.allsrcindx
    src = allsrc[srcind]
    src.fid = strtrim(src.fid,2)
    ; Correct for exptime and aperture correction
    src.cmag = src.mag+2.5*alog10(chstr1.exptime)-chstr1.apcor
    expsrc[cnt:cnt+chstr1.nsrc-1] = src
    cnt += chstr1.nsrc
  endfor

  ; Match to MALLOBJ
  MATCH,mallobj.id,expsrc.fid,ind1b,ind2b,/sort,count=nmatchb
  print,strtrim(nmatchb,2),' matches of this exposure to MALLOBJ'
  expsrc2 = expsrc[ind2b]
  mallobj2 = mallobj[ind1b]
  mgaia2 = mgaia[ind1b]

  ; Getting ALLOBJ mag column for this filter
  omagind = where(otags eq strupcase(chstr1.filter),nomagind)

  ; Make the image
  maglim = 21.0  ; 22.0
  gd = where(mallobj2.(omagind) lt maglim and mallobj2.g-mallobj2.i gt 0.4 and mallobj2.g-mallobj2.i lt 1.5,ngd)
  cenra = mean(minmax(mallobj2[gd].ra))
  cendec = mean(minmax(mallobj2[gd].dec))
  rotsphcen,mallobj2[gd].ra,mallobj2[gd].dec,cenra,cendec,lon,lat,/gnomic
  hess,lon,lat,expsrc2[gd].cmag-mgaia2[gd]._gmag_,im,dx=0.05,dy=0.05,/mean,xarr=xarr,yarr=yarr,/noplot

  ; Create a smoothed version of the image
  smim = gsmooth(im,2)
  wt = gsmooth(float(im ne 0),2)
  bwt = where(im eq 0,nbwt)
  if nbwt gt 0 then wt[bwt]=1.0
  smim /= wt
  smim[bwt] = 0

  gg = where(smim ne 0.0,ngg)
  med = median(smim[gg])
  sig = mad(smim[gg])
  rng = range(smim[gg])


  print,strtrim(i+1,2),'/',strtrim(nexpnum,2),' ',uexpnum[i],med,sig

  ; Save the figure
  !p.font = 0
  file = plotsdir+field+'_'+uexpnum[i]+'_spatialtransvar'
  ps_open,file,/color,thick=4,/encap
  device,/inches,xsize=9.5,ysize=9.5
  z1 = med-2.5*sig
  z2 = med+2.5*sig
  displayc,smim,xarr,yarr,maskv=0.0,maskc=0,min=z1,max=z2,tit=field+' '+uexpnum[i],xtit=textoidl('\zeta (deg)'),ytit=textoidl('\eta (deg)')
  al_legend,[chstr1.filter+' '+stringize(chstr1.exptime,ndec=1)],/top,/right,charsize=1.5,textcolor=255
  al_legend,['med='+stringize(med,ndec=2),'sig='+stringize(sig,ndec=3)],/top,/left,charsize=1.5,textcolor=255
  ps_close
  ps2jpg,file+'.eps',/eps
  spawn,['epstopdf',file+'.eps'],/noshell
  push,psfiles,file

  ; Save the image and stats
  savefile = plotsdir+field+'_'+uexpnum[i]+'_spatialtransvar.dat'
  print,'Writing information to ',savefile
  save,im,smim,xarr,yarr,med,sig,file=savefile

  ; Put in exposure structure
  expstr[i].expnum = uexpnum[i]
  expstr[i].filter = chstr1.filter
  expstr[i].exptime = chstr1.exptime
  expstr[i].totexpsrc = exptotsrc
  expstr[i].ngaiamatch = nmatchb
  expstr[i].nused = ngd
  expstr[i].med = med
  expstr[i].sig = sig
  expstr[i].range = rng

  ;stop

endfor

; Make plot of med/sigma
file = plotsdir+field+'_spatialtransvar'
ps_open,file,/color,thick=4,/encap
filtnum = intarr(nexpnum)
gu = where(expstr.filter eq 'u',ngu)
gus = where(expstr.filter eq 'u' and expstr.exptime lt 100,ngus)
gul = where(expstr.filter eq 'u' and expstr.exptime gt 100,ngul)
gg = where(expstr.filter eq 'g',ngg)
ggs = where(expstr.filter eq 'g' and expstr.exptime lt 100,nggs)
ggl = where(expstr.filter eq 'g' and expstr.exptime gt 100,nggl)
gr = where(expstr.filter eq 'r',ngr)
grs = where(expstr.filter eq 'r' and expstr.exptime lt 100,ngrs)
grl = where(expstr.filter eq 'r' and expstr.exptime gt 100,ngrl)
gi = where(expstr.filter eq 'i',ngi)
gis = where(expstr.filter eq 'i' and expstr.exptime lt 100,ngis)
gil = where(expstr.filter eq 'i' and expstr.exptime gt 100,ngil)
gz = where(expstr.filter eq 'z',ngz)
gzs = where(expstr.filter eq 'z' and expstr.exptime lt 100,ngzs)
gzl = where(expstr.filter eq 'z' and expstr.exptime gt 100,ngzl)
xr = [min(expstr.med)-0.05*range(expstr.med),max(expstr.med)+0.05*range(expstr.med)]
yr = [min(expstr.range)-0.05*range(expstr.range),max(expstr.range)+0.05*range(expstr.range)]
plot,expstr.med,expstr.range,/nodata,ps=3,xr=xr,yr=yr,xs=1,ys=1,xtit='Median offset (mag)',ytit='Spatial range in offset (mag)',tit=field
symsize = alog10(expstr.exptime)
if ngus gt 0 then oplot,[expstr[gus].med],[expstr[gus].range],ps=1,co=80
if nggs gt 0 then oplot,[expstr[ggs].med],[expstr[ggs].range],ps=1,co=120
if ngrs gt 0 then oplot,[expstr[grs].med],[expstr[grs].range],ps=1,co=150
if ngis gt 0 then oplot,[expstr[gis].med],[expstr[gis].range],ps=1,co=200
if ngzs gt 0 then oplot,[expstr[gzs].med],[expstr[gzs].range],ps=1,co=250
if ngul gt 0 then oplot,[expstr[gul].med],[expstr[gul].range],ps=4,co=80
if nggl gt 0 then oplot,[expstr[ggl].med],[expstr[ggl].range],ps=4,co=120
if ngrl gt 0 then oplot,[expstr[grl].med],[expstr[grl].range],ps=4,co=150
if ngil gt 0 then oplot,[expstr[gil].med],[expstr[gil].range],ps=4,co=200
if ngzl gt 0 then oplot,[expstr[gzl].med],[expstr[gzl].range],ps=4,co=250
robust_mean,expstr[gu].med,umn
;umed = median(expstr[gu].med)
oplot,[0,0]+umn,[-100,100],linestyle=2,col=80,thick=3
robust_mean,expstr[gg].med,gmn
;gmed = median(expstr[gg].med)
oplot,[0,0]+gmn,[-100,100],linestyle=2,col=120,thick=3
robust_mean,expstr[gr].med,rmn
;rmed = median(expstr[gr].med)
oplot,[0,0]+rmn,[-100,100],linestyle=2,col=150,thick=3
robust_mean,expstr[gi].med,imn
;imed = median(expstr[gi].med)
oplot,[0,0]+imn,[-100,100],linestyle=2,col=200,thick=3
robust_mean,expstr[gz].med,zmn
;zmed = median(expstr[gz].med)
oplot,[0,0]+zmn,[-100,100],linestyle=2,col=250,thick=3
al_legend,['u ('+strtrim(ngu,2)+')','g ('+strtrim(ngg,2)+')','r ('+strtrim(ngr,2)+')',$
           'i ('+strtrim(ngi,2)+')','z ('+strtrim(ngz,2)+')'],textcolor=[80,120,150,200,250],/top,/left,charsize=1.5
al_legend,['Short','Long'],psym=[1,4],/top,/right,charsize=1.5
ps_close
ps2jpg,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
psfiles = [file,psfiles]

; Create combined PDF file
combpdffile = field+'_spatialtransvar_comb.pdf'
cd,current=curdir
cd,plotsdir
cmd = 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+combpdffile+' '
cmd += strjoin(file_basename(psfiles)+'.pdf',' ')
spawn,cmd,out,errout
cd,curdir

; Save the exposure structure for this field
MWRFITS,expstr,outcatfile,/create

stop
if keyword_set(stp) then stop

end
