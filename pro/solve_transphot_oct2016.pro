pro solve_transphot_oct2016,file,fitcolam=fitcolam,errlim=errlim,arr=arr,mstr=mstr_fixcolzpam,resid=resid_fixcolzpam,$
                         rejected=rejected_fixcolzpam,noreject=noreject,fitchipzpmjdterm=fitchipzpmjdterm,$
                         fixchipzpterm=fixchipzpterm,nosmith=nosmith,noexpreject=noexpreject,stp=stp

RESOLVE_ROUTINE,'solve_trans',/compile_full_file

; Derive transformation equations for the SMASH Oct 29-31, 2016
; observing run.

if n_elements(file) eq 0 then begin
  print,'Syntax - solve_transphot_oct2016,file'
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

print,'RUNNING SOLVE_TRANSPHOT_OCT2016 ON  >> ',file,' <<'

reduxdir = '/data/smash/cp/red/photred/'
plotsdir = reduxdir+'stdred/plots/'

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

; Loading the weather conditions table
conditions = importascii('smash_observing_conditions.txt',/header)
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
SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr0,fitstr0,resid=resid0,expstr=expstr0,rejected=rejected0,/bootstrap,$
                          noreject=noreject,noexpreject=noexpreject
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
colstr = mrdfits('/data/smash/cp/red/photred/stdred/smashred_transphot_eqns.fits',2)
gdcol = where(colstr.filter eq filter,ngdcol)
colstr = colstr[gdcol]

stop


; 3.) Determine median zpterm and color terms per chip
;------------------------------------------------------
print,'' & print,'3.) Determining median zpterm and colterm per chip'
ui = uniq(fitstr0.chip,sort(fitstr0.chip))  ; unique chips
uchips = fitstr0[ui].chip
nuchips = n_elements(uchips)
; Measure median and RMS per chip
if keyword_set(fitchipzpmjdterm) then begin
  chipstr = replicate({chip:0L,zpterm:99.0,zptermerr:99.0,zpmjdterm:99.0,zpmjdtermerr:99.0,zptermsig:99.0,colterm:99.0,coltermerr:99.0,$
                       coltermsig:99.0,rms:99.0,medrms:99.0},nuchips)
endif else begin
  chipstr = replicate({chip:0L,zpterm:99.0,zptermerr:99.0,zptermsig:99.0,colterm:99.0,coltermerr:99.0,coltermsig:99.0,rms:99.0,medrms:99.0},nuchips)
endelse
fitstr0.relzpterm = fitstr0.zpterm        ; initialize with zpterm
print,' CHIP   ZPTERM   ZPTERMERR  COLTERM  COLTERMER    RMS'
for i=0,nuchips-1 do begin
  chipstr[i].chip = uchips[i]
  MATCH,fitstr0.chip,uchips[i],ind,ind0,count=nind,/sort
  if nind gt 0 then begin
    ; Remove median zpterm for each night
    MATCH,mstr0.mjd,fitstr0[ind].mjd,ind1,ind2
    fitstr0[ind[ind2]].relzpterm -= mstr0[ind1].zpterm  ; remove the nightly median zeropoint
    ; Only use "good" value to estimate median and scatter
    gdphot = where(fitstr0[ind].photometric eq 1 and fitstr0[ind].badsoln eq 0,ngdphot)
    ; Fit linear temporal variations in chip relative zeropoint term
    if keyword_set(fitchipzpmjdterm) then begin
      zpmjdcoef1 = robust_poly_fitq(fitstr0[ind[gdphot]].mjd,fitstr0[ind[gdphot]].relzpterm,1)
      ydiff1 = fitstr0[ind[gdphot]].relzpterm-poly(fitstr0[ind[gdphot]].mjd,zpmjdcoef1)
      relzptermsig1 = mad(ydiff1)
      gd2 = where(abs(ydiff1) lt 5*relzptermsig1,ngd2)
      zpmjdcoef = dln_poly_fit(fitstr0[ind[gdphot[gd2]]].mjd,fitstr0[ind[gdphot[gd2]]].relzpterm,1,$
                               measure_errors=fitstr0[ind[gdphot[gd2]]].zptermerr,sigma=zpmjdcoeferr,yerror=yerror,status=status,yfit=yfit1,/bootstrap)
      ydiff = fitstr0[ind[gdphot]].relzpterm-poly(fitstr0[ind[gdphot]].mjd,zpmjdcoef)
      relzptermsig = mad(ydiff)
      chipstr[i].zpterm = zpmjdcoef[0]
      chipstr[i].zptermerr = zpmjdcoeferr[0]
      chipstr[i].zpmjdterm = zpmjdcoef[1]
      chipstr[i].zpmjdtermerr = zpmjdcoeferr[1]
      chipstr[i].zptermsig = relzptermsig
      fitstr0[ind].chiprelzpterm = fitstr0[ind].relzpterm - (chipstr[i].zpterm+fitstr0[ind].mjd*chipstr[i].zpmjdterm) ; relative to night and chip median
    endif else begin
      chipstr[i].zpterm = median(fitstr0[ind[gdphot]].relzpterm)
      chipstr[i].zptermsig = mad(fitstr0[ind[gdphot]].relzpterm)
      chipstr[i].zptermerr = chipstr[i].zptermsig / sqrt(ngdphot)
      fitstr0[ind].chiprelzpterm = fitstr0[ind].relzpterm - chipstr[i].zpterm   ; relative to night and chip medians
    endelse
    chipstr[i].colterm = median([fitstr0[ind[gdphot]].colterm])
    chipstr[i].coltermsig = mad([fitstr0[ind[gdphot]].colterm])
    chipstr[i].coltermerr = chipstr[i].coltermsig / sqrt(ngdphot)
    totresid = fitstr0[ind[gdphot]].nbrt*fitstr0[ind[gdphot]].brtrms^2
    rms = sqrt(mean(total(totresid)/total(fitstr0[ind[gdphot]].nbrt)))  ; brt rms
    chipstr[i].rms = rms
    chipstr[i].medrms = median([fitstr0[ind[gdphot]].rms])
    fitstr0[ind].chiprelcolterm = fitstr0[ind].colterm - chipstr[i].colterm
    print,chipstr[i].chip,chipstr[i].zpterm,chipstr[i].zptermerr,chipstr[i].colterm,chipstr[i].coltermerr,chipstr[i].rms,format='(I5,5F10.5)'
  endif
endfor

;stop

; Making some figures
print,'Making some figures'
setdisp,/silent
!p.font = 0
; Relative ZPTERM vs. chip
psfile = plotsdir+'transphot_'+filter+'_zpterm_chip'+pstag
ps_open,psfile,/color,thick=4,/encap
device,/inches,xsize=9.5,ysize=9.5
rnd = randomu(seed,nfitstr)*0.6-0.3
yr = [-0.1,0.1]
gd = where(fitstr0.photometric eq 1 and fitstr0.badsoln eq 0,ngd)
plotc,fitstr0[gd].chip+rnd[gd],fitstr0[gd].relzpterm,fitstr0[gd].mjd,ps=8,sym=0.5,xr=[0,63],yr=yr,xs=1,ys=1,xtit='Chip',ytit='Relative ZPTERM',$
      tit=filter+'-band zero-points   (color-coded by MJD)',pos=[0.08,0.08,0.98,0.40],colpos=[0.08,0.47,0.98,0.48]
oplot,chipstr.chip,chipstr.zpterm,ps=-1
xyouts,2,yr[1]-0.08*range(yr),'RMS = '+stringize(mad(fitstr0[gd].chiprelzpterm),ndec=4),align=0,charsize=1.1,charthick=4
plotc,fitstr0[gd].chip+rnd[gd],fitstr0[gd].relzpterm,fitstr0[gd].zptermerr,ps=8,sym=0.5,xr=[0,63],yr=yr,xs=1,ys=1,max=0.01,xtit='Chip',ytit='Relative ZPTERM',$
      tit=filter+'-band zero-points   (color-coded by ZPTERMERR)',pos=[0.08,0.58,0.98,0.90],colpos=[0.08,0.97,0.98,0.98],/noerase
oplot,chipstr.chip,chipstr.zpterm,ps=-1
xyouts,2,yr[1]-0.08*range(yr),'RMS = '+stringize(mad(fitstr0[gd].chiprelzpterm),ndec=4),align=0,charsize=1.1,charthick=4
ps_close
ps2png,psfile+'.eps',/eps
; COLTERM vs. chip
psfile = plotsdir+'transphot_'+filter+'_colterm_chip'+pstag
ps_open,psfile,/color,thick=4,/encap
device,/inches,xsize=9.5,ysize=9.5
rnd = randomu(seed,nfitstr)*0.6-0.3
;yr = [-0.3,0.3]
yr = [min(-3*chipstr.coltermsig+chipstr.colterm),max(3*chipstr.coltermsig+chipstr.colterm)]
plotc,fitstr0[gd].chip+rnd[gd],fitstr0[gd].colterm,fitstr0[gd].mjd,ps=8,sym=0.5,xr=[0,63],yr=yr,xs=1,ys=1,xtit='Chip',ytit='COLTERM',$
      tit=filter+'-band color terms   (color-coded by MJD)',pos=[0.08,0.08,0.98,0.40],colpos=[0.08,0.47,0.98,0.48]
oplot,chipstr.chip,chipstr.colterm,ps=-1
xyouts,2,yr[1]-0.08*range(yr),'RMS = '+stringize(mad(fitstr0[gd].chiprelcolterm),ndec=4),align=0,charsize=1.1,charthick=4
plotc,fitstr0[gd].chip+rnd[gd],fitstr0[gd].colterm,fitstr0[gd].coltermerr,ps=8,sym=0.5,xr=[0,63],yr=yr,xs=1,ys=1,max=0.01,xtit='Chip',ytit='COLTERM',$
      tit=filter+'-band color terms   (color-coded by COLTERMERR)',pos=[0.08,0.58,0.98,0.90],colpos=[0.08,0.97,0.98,0.98],/noerase
oplot,chipstr.chip,chipstr.colterm,ps=-1
xyouts,2,yr[1]-0.08*range(yr),'RMS = '+stringize(mad(fitstr0[gd].chiprelcolterm),ndec=4),align=0,charsize=1.1,charthick=4
ps_close
ps2png,psfile+'.eps',/eps
; CHIPRELZPTERM vs. NIGHTNUM
psfile = plotsdir+'transphot_'+filter+'_chiprelzpterm_nightnum'+pstag
ps_open,psfile,/color,thick=4,/encap
device,/inches,xsize=9.5,ysize=9.5
rnd = randomu(seed,nfitstr)*0.6-0.3
gg = where(abs(fitstr0.chiprelzpterm) lt 2)
yr = [ min(fitstr0[gg].chiprelzpterm) > (-3*stddev(fitstr0[gg].chiprelzpterm)),$
       max(fitstr0[gg].chiprelzpterm) < 3*stddev(fitstr0[gg].chiprelzpterm)]
plotc,fitstr0.nightnum+rnd,fitstr0.chiprelzpterm,fitstr0.chip,ps=8,sym=0.5,xr=[0,53],yr=yr,xs=1,ys=1,xtit='Night Number',ytit='Relative ZPTERM',$
      tit=filter+'-band relative zero-point term   (color-coded by CHIP)',pos=[0.08,0.08,0.98,0.40],colpos=[0.08,0.47,0.98,0.48]
plotc,fitstr0.nightnum+rnd,fitstr0.chiprelzpterm,fitstr0.zptermerr,ps=8,sym=0.5,xr=[0,53],yr=yr,xs=1,ys=1,max=0.01,xtit='Night Number',ytit='Relative ZPTERM',$
      tit=filter+'-band relative zero-point term   (color-coded by ZPTERMERR)',pos=[0.08,0.58,0.98,0.90],colpos=[0.08,0.97,0.98,0.98],/noerase
ps_close
ps2png,psfile+'.eps',/eps
; CHIPRELCOLTERM vs. NIGHTNUM
psfile = plotsdir+'transphot_'+filter+'_chiprelcolterm_nightnum'+pstag
ps_open,psfile,/color,thick=4,/encap
device,/inches,xsize=9.5,ysize=9.5
rnd = randomu(seed,nfitstr)*0.6-0.3
gg = where(abs(fitstr0.chiprelcolterm) lt 2)
yr = [ min(fitstr0[gg].chiprelcolterm) > (-5*stddev(fitstr0[gg].chiprelcolterm)),$
       max(fitstr0[gg].chiprelcolterm) < 5*stddev(fitstr0[gg].chiprelcolterm)]
plotc,fitstr0.nightnum+rnd,fitstr0.chiprelcolterm,fitstr0.chip,ps=8,sym=0.5,xr=[0,53],yr=yr,xs=1,ys=1,xtit='Night Number',ytit='Relative COLTERM',$
      tit=filter+'-band relative color term   (color-coded by CHIP)',pos=[0.08,0.08,0.98,0.40],colpos=[0.08,0.47,0.98,0.48]
plotc,fitstr0.nightnum+rnd,fitstr0.chiprelcolterm,fitstr0.coltermerr,ps=8,sym=0.5,xr=[0,53],yr=yr,xs=1,ys=1,max=0.01,xtit='Night Number',ytit='Relative COLTERM',$
      tit=filter+'-band relative color term   (color-coded by COLTERMERR)',pos=[0.08,0.58,0.98,0.90],colpos=[0.08,0.97,0.98,0.98],/noerase
ps_close
ps2png,psfile+'.eps',/eps


; 4.) Fixing colterm, Redetermine zpterm and airmass term
;---------------------------------------------------------
print,'' & print,'4.) Fixing colterm.  Redetermining zpterm and amterm'
SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr_fixcolr,fitstr_fixcolr,fixcolr=chipstr,resid=resid_fixcolr,rejected=rejected_fixcolr,$
                          /bootstrap,noreject=noreject,noexpreject=noexpreject

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


; 5.) Re-Determine median zpterm per chip
;------------------------------------------
print,'' & print,'5.) Re-determine median zpterm per chip'
; Measure median chip zpterm again
if keyword_set(fitchipzpmjdterm) then begin
  chipstr_fixcolr = replicate({chip:0L,zpterm:99.0,zptermerr:99.0,zpmjdterm:99.0,zpmjdtermerr:99.0,zptermsig:99.0,colterm:99.0,coltermerr:99.0,$
                               coltermsig:99.0,rms:99.0,medrms:99.0},nuchips)
endif else begin
  chipstr_fixcolr = replicate({chip:0L,zpterm:99.0,zptermerr:99.0,zptermsig:99.0,colterm:99.0,coltermerr:99.0,coltermsig:99.0,rms:99.0,medrms:99.0},nuchips)
endelse
fitstr_fixcolr.relzpterm = fitstr_fixcolr.zpterm        ; initialize with zpterm
print,' CHIP   ZPTERM   ZPTERMERR  COLTERM  COLTERMER    RMS'
for i=0,nuchips-1 do begin
  chipstr_fixcolr[i].chip = uchips[i]
  MATCH,fitstr_fixcolr.chip,uchips[i],ind,ind0,count=nind,/sort
  if nind gt 0 then begin
    ; Remove median zpterm for each night
    MATCH,mstr_fixcolr.mjd,fitstr_fixcolr[ind].mjd,ind1,ind2
    fitstr_fixcolr[ind[ind2]].relzpterm -= mstr_fixcolr[ind1].zpterm  ; remove the nightly median zeropoint
    ; Only use "good" value to estimate median and scatter
    gdphot = where(fitstr_fixcolr[ind].photometric eq 1 and fitstr_fixcolr[ind].badsoln eq 0,ngdphot)
    ; Fit linear temporal variations in chip relative zeropoint term
    if keyword_set(fitchipzpmjdterm) then begin
      zpmjdcoef1 = robust_poly_fitq(fitstr_fixcolr[ind[gdphot]].mjd,fitstr_fixcolr[ind[gdphot]].relzpterm,1)
      ydiff1 = fitstr_fixcolr[ind[gdphot]].relzpterm-poly(fitstr_fixcolr[ind[gdphot]].mjd,zpmjdcoef1)
      relzptermsig1 = mad(ydiff1)
      gd2 = where(abs(ydiff1) lt 5*relzptermsig1,ngd2)
      zpmjdcoef = dln_poly_fit(fitstr_fixcolr[ind[gdphot[gd2]]].mjd,fitstr_fixcolr[ind[gdphot[gd2]]].relzpterm,1,$
                               measure_errors=fitstr_fixcolr[ind[gdphot[gd2]]].zptermerr,sigma=zpmjdcoeferr,yerror=yerror,status=status,yfit=yfit1,/bootstrap)
      ydiff = fitstr_fixcolr[ind[gdphot]].relzpterm-poly(fitstr_fixcolr[ind[gdphot]].mjd,zpmjdcoef)
      relzptermsig = mad(ydiff)
      chipstr_fixcolr[i].zpterm = zpmjdcoef[0]
      chipstr_fixcolr[i].zptermerr = zpmjdcoeferr[0]
      chipstr_fixcolr[i].zpmjdterm = zpmjdcoef[1]
      chipstr_fixcolr[i].zpmjdtermerr = zpmjdcoeferr[1]
      chipstr_fixcolr[i].zptermsig = relzptermsig
      fitstr_fixcolr[ind].chiprelzpterm = fitstr_fixcolr[ind].relzpterm - (chipstr_fixcolr[i].zpterm+fitstr_fixcolr[ind].mjd*chipstr_fixcolr[i].zpmjdterm) ; relative to night and chip median
    endif else begin
      chipstr_fixcolr[i].zpterm = median(fitstr_fixcolr[ind[gdphot]].relzpterm)
      chipstr_fixcolr[i].zptermsig = mad(fitstr_fixcolr[ind[gdphot]].relzpterm)
      chipstr_fixcolr[i].zptermerr = chipstr_fixcolr[i].zptermsig / sqrt(ngdphot)
      fitstr_fixcolr[ind].chiprelzpterm = fitstr_fixcolr[ind].relzpterm - chipstr_fixcolr[i].zpterm   ; relative to night and chip medians
    endelse
    chipstr_fixcolr[i].colterm = median([fitstr_fixcolr[ind[gdphot]].colterm])
    chipstr_fixcolr[i].coltermsig = mad([fitstr_fixcolr[ind[gdphot]].colterm])
    chipstr_fixcolr[i].coltermerr = chipstr_fixcolr[i].coltermsig / sqrt(ngdphot)
    totresid = fitstr_fixcolr[ind[gdphot]].nbrt*fitstr_fixcolr[ind[gdphot]].brtrms^2
    rms = sqrt(mean(total(totresid)/total(fitstr_fixcolr[ind[gdphot]].nbrt)))  ; brt rms
    chipstr_fixcolr[i].rms = rms
    chipstr_fixcolr[i].medrms = median([fitstr_fixcolr[ind[gdphot]].rms])
    fitstr_fixcolr[ind].chiprelcolterm = fitstr_fixcolr[ind].colterm - chipstr_fixcolr[i].colterm
    print,chipstr_fixcolr[i].chip,chipstr_fixcolr[i].zpterm,chipstr_fixcolr[i].zptermerr,chipstr_fixcolr[i].colterm,$
          chipstr_fixcolr[i].coltermerr,chipstr_fixcolr[i].rms,format='(I5,5F10.5)'
    ;; Checking for temporal trends in relative zpterm
    ;zpcoef = robust_poly_fitq(fitstr_fixcolr[ind[gdphot]].mjd,fitstr_fixcolr[ind[gdphot]].chiprelzpterm,1)
    ;plotc,fitstr_fixcolr[ind[gdphot]].mjd,fitstr_fixcolr[ind[gdphot]].chiprelzpterm,ps=8,yerr=fitstr_fixcolr[ind[gdphot]].zptermerr,xs=1,tit='chip = '+strtrim(uchips[i],2)
    ;xx = scale_vector(findgen(100),min(umjd),max(umjd))
    ;oplot,xx,poly(xx,zpcoef),co=250
    ;print,i+1,zpcoef[1],zpcoef[1]*range(umjd)
    ;wait,1
  endif
endfor
; COLTERMERR = 0.0 because it was fixed

;stop


; 6.) Redetermine NIGHTLY zpterm and amterm with chip-level relative zpterm 
;      and colterm fixed
;---------------------------------------------------------------------------
print,'' & print,'6.) Fixing chip-level relative zpterm and colterm.  Redetermining NIGHTLY zpterm and amterm'
if keyword_set(fixchipzpterm) then begin
  SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr_fixcolzp,fitstr_fixcolzp,fixcolr=chipstr,fixchipzp=chipstr_fixcolr,resid=resid_fixcolzp,$
                            rejected=rejected_fixcolzp,/bootstrap,noreject=noreject,noexpreject=noexpreject
endif else begin
  print,'NOT FIXING CHIP-LEVEL ZERO-POINT TERMS'
  mstr_fixcolzp = mstr_fixcolr
  fitstr_fixcolzp = fitstr_fixcolr
  resid_fixcolzp = resid_fixcolr
  rejected_fixcolzp = rejected_fixcolr
endelse

; Reset PHOTOMETRIC tag
MATCH,mstr_fixcolzp.mjd,conditions.mjd,ind1,ind2
mstr_fixcolzp[ind1].photometric = conditions[ind2].photometric
for i=0,n_elements(conditions)-1 do begin
  MATCH,fitstr_fixcolzp.mjd,conditions[i].mjd,ind1,ind2,count=nmatch,/sort
  if nmatch gt 0 then fitstr_fixcolzp[ind1].photometric=conditions[i].photometric
endfor

; Flag "bad" solution with not Not enough stars
bdmjd = where(mstr_fixcolzp.nstars lt 100,nbdmjd,comp=gdmjd,ncomp=ngdmjd)
if nbdmjd gt 0 then mstr_fixcolzp[bdmjd].badsoln = 1    ; not enough stars
if ngdmjd gt 0 then mstr_fixcolzp[gdmjd].badsoln = 0    ; enough stars
fitstr_fixcolzp.badsoln = 0
for i=0,nbdmjd-1 do begin
  MATCH,fitstr_fixcolzp.mjd,mstr_fixcolzp[bdmjd[i]].mjd,ind1,ind2,count=nmatch,/sort
  fitstr_fixcolzp[ind1].badsoln = 1
endfor

;stop

; 7.) Averaging airmass terms
;----------------------------
; Bad for non-photometric nights 2, 40, 50
; Bad for bad solution nights 8, 9
; Bad for 38
; 6, 15, 22, 24 are a bit high

; i-band: 6, 14, 23, 37, 39, 49  these are off
;      56510 (clear/phot), 56681 (cloudy), 56702 (clear), 57097 (a few
;      cloud, crazy standard data), 57099 (nonphot), 57435 (some clouds)
; most are non-phot or bad soln, but only one that's "good" is off 56510 (clear/phot)

averageamterm:

; the airmass terms DO get more consistent/better after the various
; parameter "fixings"

; Fix airmass terms for nights with low airmass range
print,'' & print,'7.) Averaging airmass terms'
ntstr_fixcolzp = mstr_fixcolzp
add_tag,ntstr_fixcolzp,'amavgflag',0,ntstr_fixcolzp  ; amterm averaged
add_tag,ntstr_fixcolzp,'namavg',0,ntstr_fixcolzp  ; amterm averaged
print,'  NUM   MJD    AMTERM   AMTERMSIG FLAG  NAVG   COMMENT'
for i=0,nmjds-1 do begin
  ; Only look at photometric nights
  if mstr_fixcolzp[i].photometric eq 1 then begin
    ; Low airmass range or bad solution
    if mstr_fixcolzp[i].amrange lt 0.35 or mstr_fixcolzp[i].badsoln eq 1 then begin
      ; Get nights close in time with good airmass range
      si = sort(abs(mstr_fixcolzp.mjd-mstr_fixcolzp[i].mjd))
      gdnei = where(mstr_fixcolzp[si].amrange gt 0.35 and mstr_fixcolzp[si].badsoln eq 0,ngdnei)
      nei = si[gdnei[0:3]]
      ; Weights, time difference and amterm uncertainty
      wt = 1.0 / abs(mstr_fixcolzp[nei].mjd-mstr_fixcolzp[i].mjd)
      wt *= 1.0/mstr_fixcolzp[nei].amtermsig^2
      wt /= total(wt)   ; normalize
      ; Computed weighted mean
      xmn = total(mstr_fixcolzp[nei].amterm * wt) / total(wt)
      ; error of weighted mean
      ;  https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
      ;  "Statistical Properties" section
      xsig = sqrt( total( wt^2 * mstr_fixcolzp[nei].amtermsig^2 ) )
      ; Stuff it in the structure
      ntstr_fixcolzp[i].amterm = xmn
      ntstr_fixcolzp[i].amtermsig = xsig
      ntstr_fixcolzp[i].amavgflag = 1      ; we averaged it the AMTERM values
      ntstr_fixcolzp[i].namavg = n_elements(nei)    ; number averaged
      comment = '  low amrange'

    ; Refit AMTERM with neighboring nights's data
    endif else begin ; low amrange

      ; Get nights close in time with good airmass range
      diff_mjd = abs(mstr_fixcolzp.mjd-mstr_fixcolzp[i].mjd)
      gdnei1 = where(mstr_fixcolzp.amrange gt 0.35 and diff_mjd le 30 and mstr_fixcolzp.badsoln eq 0 and mstr_fixcolzp.photometric eq 1,ngdnei1)

      ; Enough to average
      if ngdnei1 gt 1 then begin

        ; Take the nearest 4 neighbor nights (and the current one)
        si = sort(diff_mjd[gdnei1])
        gdnei = gdnei1[si[0:(4<(ngdnei1-1))]]
        ngdnei = n_elements(gdnei)

        undefine,useind,resid
        for j=0,ngdnei-1 do begin
          MATCH,arr.mjd,mstr_fixcolzp[gdnei[j]].mjd,ind1,ind2,count=nmatch,/sort
          if nmatch gt 0 then begin
            push,useind,ind1
            ; Add AMTERM back into residuals, but subtract 1
            ;   so that airmass=1 points will be around zero
            push,resid,resid_fixcolzp[ind1]+(arr[ind1].airmass-1)*mstr_fixcolzp[gdnei[j]].amterm
          endif
        endfor
        coef1 = robust_poly_fitq(arr[useind].airmass,resid,1)
        sig1 = mad(resid-poly(arr[useind].airmass,coef1))
        ; cut outliers and stars previously rejected
        gdind = where(abs(resid-poly(arr[useind].airmass,coef1)) lt 4.0*sig1 and rejected_fixcolzp[useind] eq 0,ngdind)
        coef = dln_poly_fit(arr[useind[gdind]].airmass,resid[gdind],1,measure_errors=arr[useind[gdind]].err,sigma=coeferr,yerror=yerror,status=status,/bootstrap)
        ntstr_fixcolzp[i].amterm = coef[1]
        ntstr_fixcolzp[i].amtermsig = coeferr[1]
        ntstr_fixcolzp[i].amavgflag = 2      ; we determine AMTERM from actual points
        ntstr_fixcolzp[i].namavg = ngdnei    ; number averaged
        comment = ''
      endif else comment='  not enough nearby nights to average'
    endelse  ; average

  ; Non-photometric
  ;   Do NOT determine/set airmass or zpterm for NON-PHOTOMETRIC nights
  endif else begin
    ntstr_fixcolzp[i].amterm = 0.0
    ntstr_fixcolzp[i].amtermsig = 0.0
    comment = '  non-photometric'
  endelse
  print,i+1,ntstr_fixcolzp[i].mjd,ntstr_fixcolzp[i].amterm,ntstr_fixcolzp[i].amtermsig,ntstr_fixcolzp[i].amavgflag,$
        ntstr_fixcolzp[i].namavg,comment,format='(I5,I7,2F10.5,I5,I5,A-30)'
endfor

;plotc,mstr_fixcolzp.nightnum,mstr_fixcolzp.amterm,mstr_fixcolzp.amrange,ps=1
;oplot,mstr_fixcolr.nightnum,mstr_fixcolr.amterm,ps=4,co=250
;oplot,mstr0.nightnum,mstr0.amterm,ps=6,co=150
;
;plotc,mstr_fixcolzp.mjd,mstr_fixcolzp.amterm,mstr_fixcolzp.amrange,ps=1,xs=1
;oplot,mstr_fixcolr.mjd,mstr_fixcolr.amterm,ps=4,co=250
;oplot,mstr0.mjd,mstr0.amterm,ps=6,co=150


; Make a figure of AMTERM
psfile = plotsdir+'transphot_'+filter+'_amterm'+pstag
ps_open,psfile,/color,thick=4,/encap
device,/inches,xsize=9.5,ysize=9.5
gphot = where(ntstr_fixcolzp.photometric eq 1,ngphot)
xr = minmax(ntstr_fixcolzp[gphot].mjd)
xr = [xr[0]-10,xr[1]+10]
plotc,ntstr_fixcolzp[gphot].mjd,ntstr_fixcolzp[gphot].amterm,ntstr_fixcolzp[gphot].amrange,ps=7,sym=1.5,pos=[0.08,0.52,0.97,0.91],colpos=[0.08,0.97,0.97,0.99],$
      xr=xr,xs=1,yerr=ntstr_fixcolzp[gphot].amtermsig,xtit='MJD',ytit='AMTERM',tit='AMTERM (photometric nights, color-coded by airmass range)'
oplot,mstr_fixcolzp[gphot].mjd,mstr_fixcolzp[gphot].amterm,ps=1
al_legend,['Original value','Averaged value'],psym=[1,7],color=[0,250],/top,/left,charsize=1.0
xr2 = minmax(ntstr_fixcolzp[gphot].nightnum)
xr2 = [xr2[0]-1,xr2[1]+1]
plotc,ntstr_fixcolzp[gphot].nightnum,ntstr_fixcolzp[gphot].amterm,ntstr_fixcolzp[gphot].amrange,ps=7,sym=1.5,pos=[0.08,0.07,0.97,0.46],/nocolorbar,/noerase,$
      xr=xr2,xs=1,yerr=ntstr_fixcolzp[gphot].amtermsig,xtit='Night Number',ytit='AMTERM'
oplot,mstr_fixcolzp[gphot].nightnum,mstr_fixcolzp[gphot].amterm,ps=1
ps_close
ps2png,psfile+'.eps',/eps

; Combine filter-level summary plots
; transphot_g_zpterm_chip.png
; transphot_g_colterm_chip.png
; transphot_g_chiprelzpterm_nightnum.png
; transphot_g_chiprelcolterm_nightnum.png
; transphot_g_amterm.png
cd,current=curdir
cd,plotsdir
figfiles = 'transphot_'+filter+'_'+['zpterm_chip','colterm_chip','chiprelzpterm_nightnum','chiprelcolterm_nightnum','amterm']
for i=0,n_elements(figfiles)-1 do spawn,['epstopdf',figfiles[i]+'.eps'],/noshell
; Combine
print,'Writing combined filter-level summary plots to transphot_'+filter+'_comb'+pstag+'.pdf'
cmd = 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=transphot_'+filter+'_comb'+pstag+'.pdf '
cmd += strjoin(figfiles+'.pdf',' ')
spawn,cmd,out,errout
cd,curdir


; 8.) Redetermine NIGHTLY zpterm holding everything else fixed
;--------------------------------------------------------------
print,'' & print,'8.) Redetermine NIGHTLY zpterm holding everything else fixed'
if keyword_set(fixchipzpterm) then begin
  SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr_fixcolzpam,fitstr_fixcolzpam,fixcolr=chipstr,fixchipzp=chipstr_fixcolr,$
                          fixam=ntstr_fixcolzp,resid=resid_fixcolzpam,rejected=rejected_fixcolzpam,/bootstrap,$
                          noreject=noreject,noexpreject=noexpreject
endif else begin
  print,'NOT FIXING CHIP-LEVEL ZERO-POINT TERMS'
  SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr_fixcolzpam,fitstr_fixcolzpam,fixcolr=chipstr,$
                          fixam=ntstr_fixcolzp,resid=resid_fixcolzpam,rejected=rejected_fixcolzpam,/bootstrap,$
                          noreject=noreject,noexpreject=noexpreject
endelse

; Reset PHOTOMETRIC tag
MATCH,mstr_fixcolzpam.mjd,conditions.mjd,ind1,ind2
mstr_fixcolzpam[ind1].photometric = conditions[ind2].photometric
for i=0,n_elements(conditions)-1 do begin
  MATCH,fitstr_fixcolzpam.mjd,conditions[i].mjd,ind1,ind2,count=nmatch,/sort
  if nmatch gt 0 then fitstr_fixcolzpam[ind1].photometric=conditions[i].photometric
endfor

; Flag "bad" solution with Not enough stars
bdmjd = where(mstr_fixcolzpam.nstars lt 100,nbdmjd,comp=gdmjd,ncomp=ngdmjd)
if nbdmjd gt 0 then mstr_fixcolzpam[bdmjd].badsoln = 1    ; not enough stars
if ngdmjd gt 0 then mstr_fixcolzpam[gdmjd].badsoln = 0    ; enough stars
fitstr_fixcolzpam.badsoln = 0
for i=0,nbdmjd-1 do begin
  MATCH,fitstr_fixcolzpam.mjd,mstr_fixcolzpam[bdmjd[i]].mjd,ind1,ind2,count=nmatch,/sort
  fitstr_fixcolzpam[ind1].badsoln = 1
endfor
bdfitstr = where(fitstr_fixcolzpam.nstars lt 3 and fitstr_fixcolzpam.zpterm gt 100,nbdfitstr)
if nbdfitstr gt 0 then fitstr_fixcolzpam[bdfitstr].badsoln = 1

; Make final residual figures
print,'Making final residual plots'
undefine,residfigs
for i=0,nmjds-1 do begin

  MATCH,arr.mjd,mstr_fixcolzpam[i].mjd,ind1,ind2,count=nmatch,/sort
  gd = where(rejected_fixcolzp[ind1] eq 0,ngd)
  if ngd gt 0 then ind = ind1[gd] else ind=-1
  nind = ngd

  ; Make resid plots, one per night
  ;  resid vs. color, resid vs. airmass, resid vs. chip number

  setdisp,/silent
  psfile = plotsdir+'transphot_'+filter+'_'+strtrim(mstr_fixcolzpam[i].mjd,2)+'_resid'+pstag
  ps_open,psfile,/color,thick=4,/encap
  device,/inches,xsize=9.5,ysize=9.5
  x0 = 0.08 & x1=0.97
  y0 = 0.05 & dy = 0.30
  sig = mad(resid_fixcolzpam[ind])
  ;yr = minmax(resid_fixcolzpam[ind])
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
  plotc,[arr[ind].(8)],[resid_fixcolzpam[ind]],[arr[ind].err],ps=psym,symsize=symsize,xr=xr,yr=yr,xs=1,ys=1,xtit=colorname,ytit='Residuals',$
        pos=[x0,y0,x1,dy],colpos=[x0,0.97,x1,0.99],xticklen=0.04,min=zr[0],max=zr[1]
  xyouts,xr[0]+0.03*range(xr),yr[1]-0.12*range(yr),'<COLTERM> = '+stringize(mstr_fixcolzpam[i].colterm,ndec=4)+'+/-'+$
         stringize(mstr_fixcolzpam[i].coltermsig,ndec=4),align=0,charsize=1.3
  oplot,xr,[0,0],linestyle=2
  ; Resid vs. airmass
  xr = minmax(arr[ind].airmass)
  xr = [0.9,(xr[1]+0.05) > 2.0]
  plotc,[arr[ind].airmass],[resid_fixcolzpam[ind]],[arr[ind].err],ps=psym,symsize=symsize,xr=xr,yr=yr,xs=1,ys=1,xtit='Airmass',ytit='Residuals',$
        pos=[x0,y0+dy,x1,2*dy],/nocolorbar,xticklen=0.04,/noerase,min=zr[0],max=zr[1]
  xyouts,xr[0]+0.03*range(xr),yr[1]-0.12*range(yr),'AMTERM = '+stringize(mstr_fixcolzpam[i].amterm,ndec=4)+'+/'+$
         stringize(ntstr_fixcolzp[i].amtermsig,ndec=4),align=0,charsize=1.3
  if mstr_fixcolzpam[i].photometric eq 1 then photcom='PHOTOMETRIC' else photcom='NON-PHOTOMETRIC'
  if mstr_fixcolzpam[i].photometric eq 1 then col=0 else col=250
  xyouts,xr[1]-0.03*range(xr),yr[1]-0.12*range(yr),photcom,align=1,charsize=1.3,color=col
  oplot,xr,[0,0],linestyle=2
  ; Resid vs. chip number
  if nind gt 0 then rnd = randomu(seed,nind)*0.6-0.3 else rnd=0.0
  xr = [0,63]
  plotc,[arr[ind].chip+rnd],[resid_fixcolzpam[ind]],[arr[ind].err],ps=psym,symsize=symsize,xr=xr,yr=yr,xs=1,ys=1,xtit='Chip',ytit='Residuals',/noerase,$
        xticklen=0.04,tit=strtrim(mstr_fixcolzpam[i].mjd,2)+' - Residuals vs. color/airmass/chip (color-coded by ERR)',pos=[x0,y0+2*dy,x1,3*dy],$
        min=zr[0],max=zr[1],/nocolorbar
  xyouts,xr[0]+0.03*range(xr),yr[1]-0.12*range(yr),'ZPTERM = '+stringize(mstr_fixcolzpam[i].zpterm,ndec=4)+'+/-'+$
         stringize(mstr_fixcolzpam[i].zptermsig,ndec=4),align=0,charsize=1.3
  xyouts,xr[1]-0.03*range(xr),yr[1]-0.12*range(yr),'RMS = '+stringize(mstr_fixcolzpam[i].rms,ndec=4),align=1,charsize=1.3
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
  fieldstr[i].nrejected = total(rejected_fixcolzpam[ind])
  fieldstr[i].medresid = median(resid_fixcolzpam[ind])
  fieldstr[i].sigresid = mad(resid_fixcolzpam[ind])
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
  expstr[i].nrejected = total(rejected_fixcolzpam[ind])
  expstr[i].medresid = median(resid_fixcolzpam[ind])
  expstr[i].sigresid = mad(resid_fixcolzpam[ind])
endfor

;stop


; 9.) Making final output transformation equation structures
;------------------------------------------------------------
print,'' & print,'9.) Making final transformation equation structures'

; Information for each chip+night combination
; most of the info is in FITSTR_FIXCOLZPAM
fitstr = replicate({mjd:-1L,chip:-1L,nightnum:0,filter:'',color:'',colband:'',colsign:0,nstars:0L,nrejected:0L,seeing:99.9,zpterm:999.0,zptermerr:999.0,$
                     colterm:999.0,coltermerr:999.0,amterm:0.0,amtermerr:999.0,colamterm:0.0,colamtermerr:999.0,$
                     rms:999.0,sig:999.0,chisq:999.0,medresid:999.0,nbrt:0L,brtrms:999.0,brtsig:999.0,brtchisq:999.0,$
                     badsoln:-1,photometric:-1,amavgflag:0,namavg:0},n_elements(fitstr_fixcolzpam))
struct_assign,fitstr_fixcolzpam,fitstr
fitstr.filter = filter
fitstr.colamtermerr = 0.0
; Copy amtermerr over
for i=0,nmjds-1 do begin
  MATCH,fitstr.mjd,ntstr_fixcolzp[i].mjd,ind1,ind2,/sort,count=nmatch
  fitstr[ind1].amtermerr = ntstr_fixcolzp[ind2[0]].amtermsig
  fitstr[ind1].amavgflag = ntstr_fixcolzp[ind2[0]].amavgflag
  fitstr[ind1].namavg = ntstr_fixcolzp[ind2[0]].namavg
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
;   CHIPSTR has the fixed color terms
;   CHIPSTR_FIXCOLR has the relative zeropoint terms
fchipstr = chipstr_fixcolr
fchipstr.colterm = chipstr.colterm
fchipstr.coltermsig = chipstr.coltermsig
fchipstr.coltermerr = chipstr.coltermerr
add_tag,fchipstr,'filter',filter,fchipstr
add_tag,fchipstr,'color',colorname,fchipstr
add_tag,fchipstr,'colband',colband,fchipstr
add_tag,fchipstr,'colsign',colsign,fchipstr

; Nightly zeropoint and airmass terms
;  NTSTR_FIXCOLZP has the fixed airmass terms
;  MSTR_FIXCOLZPAM has the nightly zero-point terms
fntstr = replicate({mjd:0L,nightnum:0,date:'',nchips:0L,filter:'',nstars:0L,seeing:99.9,rms:99.0,brtrms:99.0,airmass0:0.0,airmass1:0.0,amrange:0.0,$
                    zpterm:99.0,zptermsig:99.0,amterm:99.0,amtermsig:99.0,colamterm:99.0,colamtermsig:99.0,$
                    badsoln:-1,photometric:-1,amavgflag:0,namavg:0},nmjds)
struct_assign,mstr_fixcolzpam,fntstr
fntstr.filter = filter
; Copy the AMTERM info from NTSTR_FIXCOLZP
fntstr.amterm = ntstr_fixcolzp.amterm
fntstr.amtermsig = ntstr_fixcolzp.amtermsig
fntstr.amavgflag = ntstr_fixcolzp.amavgflag
fntstr.namavg = ntstr_fixcolzp.namavg
add_tag,fntstr,'color',colorname,fntstr
add_tag,fntstr,'colband',colband,fntstr
add_tag,fntstr,'colsign',colsign,fntstr

; NEED TO ADD IN INFO FOR NIGHTS WITH NO STANDARD STAR DATA!!!
;  or do we?  how should we handle these??
; 20140105, could take nightly zero-point and amterm from neighboring
; nights, average them.


; DO I NEED A LINE FOR EVERY NIGHT+CHIP???

; Deal with nights that don't have enough stars for their own solutions


; -- FINAL OUTPUT --
; 1st extension: chip-level structure with color term (+err) and
;     zpterm (+err)
; 2nd extension: night-level structure with zpter (+err) and
;     amterm (+err) and photometric flag
outfile = 'transphot_'+filter+'_eqns.fits'
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

save,arr,fitstr0,fitstr_fixcolr,fitstr_fixcolzp,fitstr_fixcolzpam,$
     mstr_fixcolzpam,resid_fixcolzpam,rejected_fixcolzpam,$
     expstr,fieldstr,file='solve_transphot_'+filter+'_resid.dat'

if keyword_set(stp) then stop

end
