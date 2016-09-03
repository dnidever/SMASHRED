pro qacheck_optfwexp,expnum,chstr,par,verbose=verbose,pl=pl,stp=stp,$
                     nsigthresh=nsigthresh,nabsthresh=nabsthresh,offstr=offstr

; Check that all of the FWHM values in the daophot opt
; files are okay and consistent

;setdisp,/silent
;loadct,39,/silent   ; this takes too long
undefine,chstr

; Not enough inputs
if n_elements(expnum) eq 0 then begin
  print,'Syntax - qacheck_optfwexp,expnum,chstr,verbose=verbose,pl=pl'
  return
endif

; Defaults
if n_elements(nsigthresh) eq 0 then nsigthresh=5.0
if n_elements(nabsthresh) eq 0 then nabsthresh=2.0

dir = file_dirname(expnum)+'/'
base = file_basename(expnum)

optfiles = file_search(dir+base+'_??.opt',count=noptfiles)
if noptfiles eq 0 then begin
  print,'NO .opt files for expnum=',expnum
  return
endif

chstr = replicate({dir:'',file:'',expnum:'',chip:-1L,fwhm:-1.0,va:-1,fitradius_fwhm:-1.0,$
                   fitfwhm:-1.0,xoff:0.0,yoff:0.0,bad:0,redo:0},noptfiles)
for i=0,noptfiles-1 do begin
  chip = long( (strsplit(file_basename(optfiles[i],'.opt'),'_',/extract))[1] )
  chstr[i].dir = dir
  chstr[i].file = optfiles[i]
  chstr[i].expnum = base
  chstr[i].chip = chip
  READCOL,optfiles[i],key,equal,value,format='A,A,F',/silent
  fwind = where(key eq 'FW',nfwind)
  if nfwind gt 0 then begin
    fwhm = value[fwind[0]]
    chstr[i].fwhm = fwhm 
  endif else fwhm=-1.0
  vaind = where(key eq 'VA',nvaind)
  if nvaind gt 0 then begin
    va = value[vaind[0]]
    chstr[i].va = va
  endif else va=-1
  fiind = where(key eq 'FI',nfiind)
  if nfiind gt 0 and nfwind gt 0 then begin
    fi = value[fiind[0]]
    fitradius_fwhm = fi / fwhm
    chstr[i].fitradius_fwhm = fitradius_fwhm
  endif else fitradius_fwhm=-1

  if keyword_set(verbose) then print,strtrim(i+1,2),' ',file_basename(optfiles[i]),' ',fwhm
endfor  ; opt file loop

med_fwhm = median(chstr.fwhm)
sig_fwhm = mad(chstr.fwhm)

; Fit linear plane to the FWHM value wrt RA/DEC
if n_elements(offstr) eq 0 then begin
  restore,'/data/smash/cp/red/photred/decam_chip_xyoff.dat'
  offstr = replicate({chip:0L,xoff:0.0,yoff:0.0},n_elements(chip))
  offstr.chip = chip
  offstr.xoff = xoff
  offstr.yoff = yoff
endif
MATCH,chstr.chip,offstr.chip,ind1,ind2,/sort
chstr[ind1].xoff = offstr[ind2].xoff
chstr[ind1].yoff = offstr[ind2].yoff

; Get estimates for X/Y slopes
xcoef = robust_poly_fitq(chstr.xoff,chstr.fwhm,1)
ycoef = robust_poly_fitq(chstr.yoff,chstr.fwhm,1)

initpars = dblarr(3)
initpars[*] = [ med_fwhm, xcoef[1], ycoef[0] ]
parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},3)
;  we can change the order by fixing certain parameters
; 0-constant, 1-x, 2-x^2, 3-x^3
; 4-x*y, 5-x^2*y, 6-x*y^2, 7-x^2*y^2
; 8-y, 9-y^2, 10-y^3
gdpts = where(abs(chstr.fwhm-med_fwhm)/sig_fwhm lt 5,ngdpts) ; outlier rejection for fitting
pars = MPFIT2DFUN('func_poly2d',chstr[gdpts].xoff,chstr[gdpts].yoff,chstr[gdpts].fwhm,$
                  chstr.fwhm*0+0.1,initpars,status=status,dof=dof,$
                  bestnorm=chisq,parinfo=parinfo,perror=perror,yfit=fitfwhm,/quiet)      
if status lt 1 then begin
  print,'Error in the fitting.  Using median value'
  pars = [med_fwhm, 0.0, 0.0]
  perror = -1
  chisq = -1
  dof = 0
  rchisq = -1
  fitfwhm = chstr.fwhm*0+med_fwhm
endif else begin
  fitfwhm = FUNC_POLY2D(chstr.xoff,chstr.yoff,pars)
  rchisq = chisq/dof
endelse
chstr.fitfwhm = fitfwhm

; Find outliers
diff_fwhm = chstr.fwhm-chstr.fitfwhm
sig_fitfwhm = mad(diff_fwhm)
bd_fwhm = where( ( abs(diff_fwhm)/sig_fitfwhm gt nsigthresh and abs(diff_fwhm) gt 0.5 ) or $
                 abs(diff_fwhm) gt nabsthresh,nbd_fwhm)
if nbd_fwhm gt 0 then begin
  if keyword_set(verbose) then print,' bad chips. ',strjoin(chstr[bd_fwhm].chip,' ')
  bdchstr = chstr[bd_fwhm]
  writecol,-1,bdchstr.file,bdchstr.fwhm,bdchstr.fwhm-med_fwhm,bdchstr.fwhm-bdchstr.fitfwhm,$
              abs(bdchstr.fwhm-bdchstr.fitfwhm)/sig_fitfwhm,fmt='(A20,4F10.3)'
  chstr[bd_fwhm].bad = 1
endif

; Check if there are any problems
if keyword_set(pl) then begin
  plot,chstr.chip,chstr.fwhm,ps=1,xr=[0,63],xs=1,tit=expnum[0]
  oplot,chstr.chip,chstr.fitfwhm,ps=4,co=80
  oplot,[0,63],[0,0]+med_fwhm,linestyle=2,co=80
  if nbd_fwhm gt 0 then oplot,[chstr[bd_fwhm].chip],[chstr[bd_fwhm].fwhm],ps=4,co=250,sym=2,thick=3
endif

if not keyword_set(silent) then print,'   ',expnum,' ',stringize(med_fwhm,ndec=2),' ',stringize(sig_fitfwhm,ndec=3),$
       '  ',strtrim(nbd_fwhm,2),' bad exposures'

if keyword_set(stp) then stop

end
