pro solve_transphot_aug2017,file,fitcolam=fitcolam,errlim=errlim,arr=arr,mstr=mstr_fixcolzpam,resid=resid_fixcolzpam,$
                         rejected=rejected_fixcolzpam,noreject=noreject,fitchipzpmjdterm=fitchipzpmjdterm,$
                         fixchipzpterm=fixchipzpterm,nosmith=nosmith,noexpreject=noexpreject,stp=stp

RESOLVE_ROUTINE,'solve_transphot',/compile_full_file

; Derive transformation equations for the SMASH August 4, 2017
; DD half-night.

if n_elements(file) eq 0 then begin
  print,'Syntax - solve_transphot_aug2017,file'
  return
endif

; Outline of the program
; 1.) Removing outlier exposures
; 2.) Run on all the nights fitting everything individually
;  -set photometric tag
;  -flag bad solutions with not enough stars
; 3.) Determine median zpterm and colterm per chip
; 4.) Fixing colterm. Redetermining zpterm and amterm
;  -set photometric tag
;  -flag bad solutions with not enough stars
; 5.) Redetermine median zpterm per chip
; 6.) Fixing chip-level relative zpterm and colterm.  Retermining
;       nightly zpterm and amterm
; 7.) Averaging airmass terms
; 8.) Redetermine nightly zpterm holding everything else fixed
; 9.) Making final chip and night-level structures and write output files

; Need to use the same color terms as previously derived

print,'RUNNING SOLVE_TRANSPHOT_AUG2017 ON  >> ',file,' <<'

rootdir = '/dl1/users/dnidever/'
reduxdir = rootdir+'smash/cp/red/photred/'
plotsdir = rootdir+'smash/cp/red/photred/stdred_20170804/plots/'

; Defaults
if n_elements(noreject) eq 0 then noreject=0         ; measurement-level rejection
if n_elements(noexpreject) eq 0 then noexpreject=1   ; exposure-level rejection
if n_elements(nosmith) eq 0 then nosmith=1
if n_elements(fixchipzpterm) eq 0 then fixchipzpterm=0
if n_elements(fitchipzpmjdterm) eq 0 then fitchipzpmjdterm=0
; Setting
print,'NOREJECT = ',keyword_set(noreject)
print,'NOEXPREJECT = ',keyword_set(noexpreject)
print,'NOSMITH = ',keyword_set(nosmith)
print,'FIXCHIPZPTERM = ',keyword_set(fixchipzpterm)
print,'FITCHIPZPMJDTERM = ',keyword_set(fitchipzpmjdterm)


; Load the data
print,'Loading the data'
arr = MRDFITS(file,1)
narr = n_elements(arr)

; Add EXPNUM if necessary
if tag_exist(arr,'expnum') eq 0 then begin
  ; Adding EXPNUM
  print,'Adding EXPNUM'
  add_tag,arr,'expnum','',arr
  dum = strsplitter(strtrim(arr.frame,2),'_',/extract)
  expnum = reform(dum[0,*])
  arr.expnum = expnum
endif

tags = tag_names(arr)
filter = strlowcase(tags[6])  ; this should be the filter name
colorname = repstr(strlowcase(tags[8]),'_','-')

; Suffix tag for plots
pstag = ''
;if keyword_set(noreject) then pstag='_norej'

; Trim stars
if n_elements(errlim) eq 0 then errlim = 0.05
;if filter eq 'u' then errlim=0.08
print,'Imposing error cut, ERR<=',strtrim(errlim,2),' mag'
gd = where(arr.err lt errlim,ngd)
arr = arr[gd]
narr = n_elements(arr)

; Remove the Smith data
if keyword_set(nosmith) then begin
  print,'#############################'
  print,'REMOVING SMITH DATA!!!!!'
  print,'#############################'
  bd = where(strmid(strtrim(arr.id,2),0,5) eq 'ugriz',nbd)
  if nbd gt 0 then remove,bd,arr
  narr = n_elements(arr)
endif

; u-band, removing very blue and very red stars
;  and removing the quadratic color shape
if filter eq 'u' then begin
  print,'Removing very blue/red stars'
  gd = where(arr.(8) ge 1.00 and arr.u_g lt 2.5,ngd)  ; u-g
  ;gd = where(arr.(8) ge 0.73 and arr.u_g lt 2.5,ngd)  ; u-g
  ;gd = where(arr.(8) ge 1.3 and arr.u_g lt 2.5,ngd)  ; u-g
  arr = arr[gd]
  narr = n_elements(arr)

  ; Removing quadratic color shape
  print,'Removing quadratic color shape'
  ucoef = [ -0.19919922d0, 0.27770109d0, -0.088840717d0]
  arr.mag -= poly(arr.(8),ucoef)
endif

; r-band, removing very red stars
if filter eq 'r' then begin
  print,'Removing very red stars'
  gd = where(arr.(8) lt 1.20,ngd) ; g-r
  arr = arr[gd]
  narr = n_elements(arr)
endif

; Adding expnum to ARR
if tag_exist(arr,'expnum') eq 0 then begin
  print,'Adding EXPNUM column to ARR'
  add_tag,arr,'EXPNUM','',arr
  dum = strsplitter(arr.frame,'_',/extract)
  expnum = reform(dum[0,*])
  arr.expnum = expnum
endif
; Ading blank seeing tag
add_tag,arr,'seeing',0.0,arr

; Loading the weather conditions table
conditions = importascii('~/projects/SMASHRED/obslog/smash_observing_conditions.txt',/header)
print,'Weather conditions'
print,'  NUM    MJD       DATE  PHOTOMETRIC     COMMENTS'
WRITECOL,-1,conditions.num,conditions.mjd,conditions.date,conditions.photometric,'   '+conditions.comments,fmt='(I5,I8,I12,I5,A-50)'


; 1.) Removing outlier exposures
;---------------------------------
; there is also an exposure outlier rejection in solve_transphot_night.pro
print,'' & print,'' & print,'1.) Removing outlier exposures'
ui = uniq(arr.expnum,sort(arr.expnum))
uexp = arr[ui].expnum
nexp = n_elements(uexp)
colcoef = robust_poly_fitq(arr.(8),arr.mag-arr.(6),1,/silent)
amcoef = robust_poly_fitq(arr.airmass,arr.mag-arr.(6),1,/silent)
expstr = replicate({expnum:'',night:0L,mjd:0L,exptime:0.0,airmass:0.0,nsources:0L,zeropt:0.0,residzeropt:0.0,residzeroptnsig:0.0},nexp)
; Get indices for each group
expsi = sort(arr.expnum)
expbrk = where(shift(arr[expsi].expnum,1) ne arr[expsi].expnum,nexpbrk)
explo = expbrk
exphi = [expbrk[1:nexpbrk-1]-1,n_elements(arr)-1]
for i=0,nexp-1 do begin
  ;MATCH,arr.expnum,uexp[i],ind1,ind2,/sort,count=nmatch
  ind = expsi[explo[i]:exphi[i]]
  nind = exphi[i]-explo[i]+1
  medoff = median([arr[ind].mag-arr[ind].(6)-arr[ind].(8)*colcoef[1]-arr[ind].airmass*amcoef[1]])
  expstr[i].expnum = uexp[i]
  ;expstr[i].night = arr[ind[0]].night
  expstr[i].mjd = arr[ind[0]].mjd
  expstr[i].exptime = arr[ind[0]].exptime
  expstr[i].airmass = median([arr[ind].airmass])
  expstr[i].nsources = nind
  expstr[i].zeropt = medoff
  ;print,i+1,uexp[i],medoff,format='(I5,A12,F10.4)'
endfor
; Compute medoff per MJD
uimjd = uniq(expstr.mjd,sort(expstr.mjd))
umjd = expstr[uimjd].mjd
nmjd = n_elements(umjd)
mjdzeropt = fltarr(nmjd)
mjdzpsig = fltarr(nmjd)
for i=0,nmjd-1 do begin
  ind = where(expstr.mjd eq umjd[i],nind)
  expstr[ind].night = i+1
  mjdzeropt[i] = median([expstr[ind].zeropt])
  mjdzpsig[i] = mad([expstr[ind].zeropt])
  expstr[ind].residzeropt = expstr[ind].zeropt - mjdzeropt[i]
endfor
; Measure sigma for each MJD
bdsig = where(mjdzpsig lt 1e-8,nbdsig,comp=gdsig)
zpsig = median(mjdzpsig[gdsig])
if nbdsig gt 0 then mjdzpsig[bdsig] = zpsig
; Measure deviation in Nsigma
;  constrain the sigma to [0.7,1.5]*median(sigma)
for i=0,nmjd-1 do begin
  ind = where(expstr.mjd eq umjd[i],nind)
  expstr[ind].residzeroptnsig = abs(expstr[ind].residzeropt) / (0.7*zpsig > mjdzpsig[i] < zpsig*1.5)
endfor
; Find bad exposures
bdexp = where(expstr.residzeroptnsig gt 5,nbdexp)
if nbdexp gt 0 then begin
  for i=0,nmjd-1 do begin
    match,umjd[i],expstr[bdexp].mjd,ind1,ind2,/sort,count=nmatch1
    if nmatch1 gt 0 then print,strtrim(umjd[i],2),'  ',strtrim(nmatch1,2),' outlier exposures: ',expstr[bdexp[ind2]].expnum
  endfor
  ;print,'Removing ',strtrim(nbdexp,2),' bad exposures: ',expstr[bdexp].expnum
endif
undefine,bdind
for i=0,nbdexp-1 do begin
  ind = expsi[explo[bdexp[i]]:exphi[bdexp[i]]]
  push,bdind,ind
endfor
if not keyword_set(noexpreject) then REMOVE,bdind,arr else print,'/NOEXPREJECT.  Not rejecting any exposures.'
narr = n_elements(arr)

; 57097  huge airmass dependence, maybe only use low-airmass exps to
;          get zpterm
; 57099  ~4 expnums are way off the others


; Removing BAD exposure data
print,'Removing data from BAD exposures'
badexp = ['00187882','00187882',$   ; g, 56372, issues with a number of chips (881 not as bad)
          '00279765','00279770','00279775',$  ; g, 56688, beg of night issues
          '00380580',$   ; g, 56984, high-airmass exposure, way off, looks bad
          '00279764','00279769','00279774',$  ; r, 56688, beg of night issues
          '00187875','00187876',$   ; r, 56372, huge scatter and causing issues
          '00187878','00187877',$   ; i, 56372, huge scatter and causing issues (877 not as bad)
          '00279763','00279768','00279773',$  ; i, 56688, beg of night issues
          '00187880',$   ; z, 56372, large scatter
          '00279762','00279767','00279772',$  ; z, 56688, beg of night issues
          '00380576',$   ; u, 56984, one exp way off, g/z-band have issues to
          '00279535','00279884','00279950','00279997']  ; u, 56687/56688, too few stars
         ; '00272500','00272501','00272504','00272505','00272508','00272509']   ; 56666 issues??
;badexp = ['00187877','00187878',$  ; there are ~10 chips that have problems, night 5+6, i-band
;          '00422929','00423077','00423026','00422153','00422939','00422158']   ; 57097, 57099 bad data, non-photometric
nbadexp = n_elements(badexp)
for i=0,nbadexp-1 do begin
  MATCH,arr.expnum,badexp[i],ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then begin
    remove,ind1,arr
    print,'Removing bad exposure: ',badexp[i]
  endif
endfor


; 2.) Run on all the nights fitting everything individually
;------------------------------------------------------------
print,'' & print,'2.) Run on all the nights fitting everything individually'
SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr0,fitstr0,resid=resid0,expstr=expstr0,rejected=rejected0,/bootstrap,noreject=noreject
nmjds = n_elements(mstr0)
uimjd = uniq(mstr0.mjd,sort(mstr0.mjd))
umjd = mstr0[uimjd].mjd
nfitstr = n_elements(fitstr0)


; Set PHOTOMETRIC tag
;--------------------
;   2   56370  20130318   24080     0.1266,  zpterm high by 0.5 mag
;        log says there were intermittent clouds
;  40   57099  20150317   64949     0.6022,  zpterm WAY HIGH by 1.0 mag
;        log says non-photometric
;  50   57435  20160216    5579     0.0409,  zpterm high by 0.3 mag
;        log says thin cirrus, non-photometric
MATCH,mstr0.mjd,conditions.mjd,ind1,ind2
mstr0[ind1].photometric = conditions[ind2].photometric
for i=0,n_elements(conditions)-1 do begin
  MATCH,fitstr0.mjd,conditions[i].mjd,ind1,ind2,count=nmatch,/sort
  if nmatch gt 0 then fitstr0[ind1].photometric=conditions[i].photometric
endfor



; Flag "bad" solution with not Not enough stars
;--------------------------------
; 8   56664  20140106, only 16 stars,  SO FEW STARS, DON'T TRUST ANY OF THE RESULTS
;      maybe determine single global zpterm for the night
; 9   56665  20140107, only 51 stars,  not much better
bdmjd = where(mstr0.nstars lt 100,nbdmjd,comp=gdmjd,ncomp=ngdmjd)
if nbdmjd gt 0 then mstr0[bdmjd].badsoln = 1    ; not enough stars
if ngdmjd gt 0 then mstr0[gdmjd].badsoln = 0    ; enough stars
fitstr0.badsoln = 0
for i=0,nbdmjd-1 do begin
  MATCH,fitstr0.mjd,mstr0[bdmjd[i]].mjd,ind1,ind2,count=nmatch,/sort
  fitstr0[ind1].badsoln = 1
endfor

; Get the previously-derived color terms per chip
chipstr = mrdfits(rootdir+'smash/cp/red/photred/stdred/smashred_transphot_eqns.fits',2)
gdchip = where(chipstr.filter eq filter,ngdchip)
chipstr = chipstr[gdchip]


; 2.) Fixing colterm, Redetermine zpterm and airmass term
;---------------------------------------------------------
print,'' & print,'2.) Fixing colterm.  Redetermining zpterm and amterm'
SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr_fixcolr,fitstr_fixcolr,fixcolr=chipstr,resid=resid_fixcolr,rejected=rejected_fixcolr,$
                          /bootstrap,noreject=noreject

; Set PHOTOMETRIC tag
MATCH,mstr_fixcolr.mjd,conditions.mjd,ind1,ind2
mstr_fixcolr[ind1].photometric = conditions[ind2].photometric
for i=0,n_elements(conditions)-1 do begin
  MATCH,fitstr_fixcolr.mjd,conditions[i].mjd,ind1,ind2,count=nmatch,/sort
  if nmatch gt 0 then fitstr_fixcolr[ind1].photometric=conditions[i].photometric
endfor

; Flag "bad" solution with not Not enough stars
bdmjd = where(mstr_fixcolr.nstars lt 100,nbdmjd,comp=gdmjd,ncomp=ngdmjd)
if nbdmjd gt 0 then mstr_fixcolr[bdmjd].badsoln = 1    ; not enough stars
if ngdmjd gt 0 then mstr_fixcolr[gdmjd].badsoln = 0    ; enough stars
fitstr_fixcolr.badsoln = 0
for i=0,nbdmjd-1 do begin
  MATCH,fitstr_fixcolr.mjd,mstr_fixcolr[bdmjd[i]].mjd,ind1,ind2,count=nmatch,/sort
  fitstr_fixcolr[ind1].badsoln = 1
endfor


;stop

; Make final residual figures
print,'Making final residual plots'
undefine,residfigs
for i=0,nmjds-1 do begin

  MATCH,arr.mjd,mstr_fixcolr[i].mjd,ind1,ind2,count=nmatch,/sort
  gd = where(rejected_fixcolr[ind1] eq 0,ngd)
  if ngd gt 0 then ind = ind1[gd] else ind=-1
  nind = ngd

  ; Make resid plots, one per night
  ;  resid vs. color, resid vs. airmass, resid vs. chip number

  setdisp,/silent
  psfile = plotsdir+'transphot_'+filter+'_'+strtrim(mstr_fixcolr[i].mjd,2)+'_resid'+pstag
  ps_open,psfile,/color,thick=4,/encap
  device,/inches,xsize=9.5,ysize=9.5
  x0 = 0.08 & x1=0.97
  y0 = 0.05 & dy = 0.30
  sig = mad(resid_fixcolr[ind])
  ;yr = minmax(resid_fixcolr[ind])
  yr = [-1,1]*(3.5*sig > 0.2)
  yr = [yr[0] > (-1.0), yr[1] < 1.0]
  zr = [0.0, max(arr.err)] 
  psym = 3
  symsize = 1.0
  if nind lt 10000 then begin
    psym = 8
    symsize = 0.1
  endif
  if nind lt 1000 then begin
    psym = 8
    symsize = 0.4
  endif
  ; Resid vs. color
  xr = minmax(arr[ind].(8))
  xr = [-0.5,2.0 < (xr[1]+0.05) > 3.5]
  plotc,[arr[ind].(8)],[resid_fixcolr[ind]],[arr[ind].err],ps=psym,symsize=symsize,xr=xr,yr=yr,xs=1,ys=1,xtit=colorname,ytit='Residuals',$
        pos=[x0,y0,x1,dy],colpos=[x0,0.97,x1,0.99],xticklen=0.04,min=zr[0],max=zr[1]
  xyouts,xr[0]+0.03*range(xr),yr[1]-0.12*range(yr),'<COLTERM> = '+stringize(mstr_fixcolr[i].colterm,ndec=4)+'+/-'+$
         stringize(mstr_fixcolr[i].coltermsig,ndec=4),align=0,charsize=1.3
  oplot,xr,[0,0],linestyle=2
  ; Resid vs. airmass
  xr = minmax(arr[ind].airmass)
  xr = [0.9,(xr[1]+0.05) > 2.0]
  plotc,[arr[ind].airmass],[resid_fixcolr[ind]],[arr[ind].err],ps=psym,symsize=symsize,xr=xr,yr=yr,xs=1,ys=1,xtit='Airmass',ytit='Residuals',$
        pos=[x0,y0+dy,x1,2*dy],/nocolorbar,xticklen=0.04,/noerase,min=zr[0],max=zr[1]
  xyouts,xr[0]+0.03*range(xr),yr[1]-0.12*range(yr),'AMTERM = '+stringize(mstr_fixcolr[i].amterm,ndec=4)+'+/'+$
         stringize(mstr_fixcolr[i].amtermsig,ndec=4),align=0,charsize=1.3
  if mstr_fixcolr[i].photometric eq 1 then photcom='PHOTOMETRIC' else photcom='NON-PHOTOMETRIC'
  if mstr_fixcolr[i].photometric eq 1 then col=0 else col=250
  xyouts,xr[1]-0.03*range(xr),yr[1]-0.12*range(yr),photcom,align=1,charsize=1.3,color=col
  oplot,xr,[0,0],linestyle=2
  ; Resid vs. chip number
  if nind gt 0 then rnd = randomu(seed,nind)*0.6-0.3 else rnd=0.0
  xr = [0,63]
  plotc,[arr[ind].chip+rnd],[resid_fixcolr[ind]],[arr[ind].err],ps=psym,symsize=symsize,xr=xr,yr=yr,xs=1,ys=1,xtit='Chip',ytit='Residuals',/noerase,$
        xticklen=0.04,tit=strtrim(mstr_fixcolr[i].mjd,2)+' - Residuals vs. color/airmass/chip (color-coded by ERR)',pos=[x0,y0+2*dy,x1,3*dy],$
        min=zr[0],max=zr[1],/nocolorbar
  xyouts,xr[0]+0.03*range(xr),yr[1]-0.12*range(yr),'ZPTERM = '+stringize(mstr_fixcolr[i].zpterm,ndec=4)+'+/-'+$
         stringize(mstr_fixcolr[i].zptermsig,ndec=4),align=0,charsize=1.3
  xyouts,xr[1]-0.03*range(xr),yr[1]-0.12*range(yr),'RMS = '+stringize(mstr_fixcolr[i].rms,ndec=4),align=1,charsize=1.3
  oplot,xr,[0,0],linestyle=2

  ps_close
  ps2png,psfile+'.eps',/eps
  push,residfigs,psfile

  ;stop

endfor

; Convert png to pdf and combine
print,'Converting residual plots to PDF'
cd,current=curdir
cd,plotsdir
for i=0,n_elements(residfigs)-1 do begin
  residfigs1 = file_basename(residfigs[i])
  spawn,['convert',residfigs1+'.png',residfigs1+'.pdf'],/noshell
endfor
; Combine
print,'Writing combined residual plots to transphot_'+filter+'_resid_all'+pstag+'.pdf'
cmd = 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=transphot_'+filter+'_resid_all'+pstag+'.pdf '
cmd += strjoin(file_basename(residfigs)+'.pdf',' ')
spawn,cmd,out,errout
cd,curdir

; --- Creating summary statistics for each standard star field ---
ui = uniq(arr.stdfield,sort(arr.stdfield))
ustdfield = arr[ui].stdfield
nstdfield = n_elements(ustdfield)
fieldstr = replicate({stdfield:'',stdnum:0L,ra:0.0d0,dec:0.0d0,sdss:0,nexp:0L,nmjd:0L,$
                      nsources:0L,nrejected:0L,medresid:0.0,sigresid:0.0},nstdfield)
; Get indices for each group
stdsi = sort(arr.stdfield)
stdbrk = where(shift(arr[stdsi].stdfield,1) ne arr[stdsi].stdfield,nstdbrk)
stdlo = stdbrk
stdhi = [stdbrk[1:nstdbrk-1]-1,n_elements(arr)-1]
for i=0,nstdfield-1 do begin
  ind = stdsi[stdlo[i]:stdhi[i]]
  nind = stdhi[i]-stdlo[i]+1
  if strmid(strupcase(ustdfield[i]),0,4) eq 'SDSS' or strmid(strupcase(ustdfield[i]),0,2) eq 'SA' then fieldstr[i].sdss=1
  uiexp = uniq(arr[ind].expnum,sort(arr[ind].expnum))
  uimjd = uniq(arr[ind].mjd,sort(arr[ind].mjd))
  fieldstr[i].stdfield = ustdfield[i]
  fieldstr[i].stdnum = i+1
  fieldstr[i].ra = median([arr[ind].ra])
  fieldstr[i].dec = median([arr[ind].dec])
  fieldstr[i].nexp = n_elements(uiexp)
  fieldstr[i].nmjd = n_elements(uimjd)
  fieldstr[i].nsources = nind
  fieldstr[i].nrejected = total(rejected_fixcolr[ind])
  fieldstr[i].medresid = median(resid_fixcolr[ind])
  fieldstr[i].sigresid = mad(resid_fixcolr[ind])
endfor


; --- Creating summary statistics for each exposure ---
ui = uniq(arr.expnum,sort(arr.expnum))
uexp = arr[ui].expnum
nexp = n_elements(uexp)
expstr = replicate({expnum:'',mjd:0L,nightnum:0L,exptime:0.0,airmass:0.0,seeing:0.0,ra:0.0d0,dec:0.0d0,$
                    stdfield:'',stdnum:0L,sdss:0,nsources:0L,nrejected:0L,medresid:0.0,sigresid:0.0},nexp)
; Get indices for each group
expsi = sort(arr.expnum)
expbrk = where(shift(arr[expsi].expnum,1) ne arr[expsi].expnum,nexpbrk)
explo = expbrk
exphi = [expbrk[1:nexpbrk-1]-1,n_elements(arr)-1]
for i=0,nexp-1 do begin
  ind = expsi[explo[i]:exphi[i]]
  nind = exphi[i]-explo[i]+1
  expstr[i].expnum = uexp[i]
  expstr[i].mjd = arr[ind[0]].mjd
  MATCH,mstr0.mjd,expstr[i].mjd,ind1,ind2,/sort
  expstr[i].nightnum = mstr0[ind1[0]].nightnum
  expstr[i].exptime = arr[ind[0]].exptime
  expstr[i].airmass = median([arr[ind].airmass])
  expstr[i].seeing = median([arr[ind].seeing])
  expstr[i].ra = median([arr[ind].ra])
  expstr[i].dec = median([arr[ind].dec])
  expstr[i].stdfield = arr[ind[0]].stdfield
  if strmid(strupcase(expstr[i].stdfield),0,4) eq 'SDSS' or strmid(strupcase(expstr[i].stdfield),0,2) eq 'SA'  then expstr[i].sdss=1
  MATCH,fieldstr.stdfield,arr[ind[0]].stdfield,ind1,ind2,/sort
  expstr[i].stdnum = fieldstr[ind1[0]].stdnum
  expstr[i].nsources = nind
  expstr[i].nrejected = total(rejected_fixcolr[ind])
  expstr[i].medresid = median(resid_fixcolr[ind])
  expstr[i].sigresid = mad(resid_fixcolr[ind])
endfor

;stop


; 3.) Making final output transformation equation structures
;------------------------------------------------------------
print,'' & print,'3.) Making final transformation equation structures'

; Information for each chip+night combination
; most of the info is in FITSTR_FIXCOLR
fitstr = replicate({mjd:-1L,chip:-1L,nightnum:0,filter:'',color:'',colband:'',colsign:0,nstars:0L,nrejected:0L,seeing:99.9,zpterm:999.0,zptermerr:999.0,$
                     colterm:999.0,coltermerr:999.0,amterm:0.0,amtermerr:999.0,colamterm:0.0,colamtermerr:999.0,$
                     rms:999.0,sig:999.0,chisq:999.0,medresid:999.0,nbrt:0L,brtrms:999.0,brtsig:999.0,brtchisq:999.0,$
                     badsoln:-1,photometric:-1,amavgflag:0,namavg:0},n_elements(fitstr_fixcolr))
struct_assign,fitstr_fixcolr,fitstr
fitstr.filter = filter
fitstr.colamtermerr = 0.0
; Copy amtermerr over
for i=0,nmjds-1 do begin
  MATCH,fitstr.mjd,mstr_fixcolr[i].mjd,ind1,ind2,/sort,count=nmatch
  fitstr[ind1].amtermerr = mstr_fixcolr[ind2[0]].amtermsig
  ;fitstr[ind1].amavgflag = mstr_fixcolr[ind2[0]].amavgflag
  ;fitstr[ind1].namavg = mstr_fixcolr[ind2[0]].namavg
endfor
fitstr.color = colorname
dum = strsplit(colorname,'-',/extract)
if dum[0] eq filter then begin  ; filter - colband
  colband = dum[1]
  colsign = 1
endif else begin                ; colband - filter
  colband = dum[0]
  colsign = -1
endelse
fitstr.colband = colband
fitstr.colsign = colsign

; Chip-level color and zero-point terms
;   CHIPSTR has everything
fchipstr = chipstr

; Nightly zeropoint and airmass terms
;  NTSTR_FIXCOLZP has the fixed airmass terms
;  MSTR_FIXCOLR has the nightly zero-point terms
fntstr = replicate({mjd:0L,nightnum:0,date:'',nchips:0L,filter:'',nstars:0L,seeing:99.9,rms:99.0,brtrms:99.0,airmass0:0.0,airmass1:0.0,amrange:0.0,$
                    zpterm:99.0,zptermsig:99.0,amterm:99.0,amtermsig:99.0,colamterm:99.0,colamtermsig:99.0,$
                    badsoln:-1,photometric:-1,amavgflag:0,namavg:0},nmjds)
struct_assign,mstr_fixcolr,fntstr
fntstr.filter = filter
; Copy the AMTERM info from NTSTR_FIXCOLZP
fntstr.amterm = mstr_fixcolr.amterm
fntstr.amtermsig = mstr_fixcolr.amtermsig
;fntstr.amavgflag = mstr_fixcolr.amavgflag
;fntstr.namavg = mstr_fixcolr.namavg
add_tag,fntstr,'color',colorname,fntstr
add_tag,fntstr,'colband',colband,fntstr
add_tag,fntstr,'colsign',colsign,fntstr

; -- FINAL OUTPUT --
; 1st extension: chip-level structure with color term (+err) and
;     zpterm (+err)
; 2nd extension: night-level structure with zpter (+err) and
;     amterm (+err) and photometric flag
outfile = 'transphot_'+filter+'_eqns_20170804.fits'
print,'Writing final photometric transformation equations to ',outfile
print,'HDU1: chip+night information'
print,'HDU2: unique chip-level information'
print,'HDU3: unique night-level information'
MKHDR,head,0
sxaddhist,'SMASH photometric transformation equuations for filter='+filter,head
sxaddhist,'HDU1: chip+night information (colterm, zeropoint)',head
sxaddhist,'HDU2: unique chip-level information (colterm, relative zeropoint)',head
sxaddhist,'HDU3: unique night-level information (amterm, nightly zeropoint)',head
sxaddpar,head,'filter',filter
sxaddhist,'NOREJECT = '+strtrim(keyword_set(noreject),2),head
sxaddhist,'NOEXPREJECT = '+strtrim(keyword_set(noexpreject),2),head
sxaddhist,'NOSMITH = '+strtrim(keyword_set(nosmith),2),head
sxaddhist,'FITCHIPZPMJDTERM = '+strtrim(keyword_set(fitchipzpmjdterm),2),head
sxaddhist,'FIXCHIPZPTERM = '+strtrim(keyword_set(fixchipzpterm),2),head
FITS_WRITE,outfile,0,head
MWRFITS,fitstr,outfile,/silent
MWRFITS,fchipstr,outfile,/silent
MWRFITS,fntstr,outfile,/silent

save,arr,fitstr0,fitstr_fixcolr,mstr_fixcolr,resid_fixcolr,rejected_fixcolr,$
     expstr,fieldstr,file='solve_transphot_'+filter+'_resid_20170804.dat'

if keyword_set(stp) then stop

end
