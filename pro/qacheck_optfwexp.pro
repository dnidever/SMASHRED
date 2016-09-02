pro qacheck_optfwexp,expnum,chstr,par,verbose=verbose,pl=pl,stp=stp,$
                     nsigthresh=nsigthresh,nabsthresh=nabsthresh

; Check that all of the FWHM values in the daophot opt
; files are okay and consistent

setdisp
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

chstr = replicate({dir:'',expnum:'',chip:-1L,fwhm:-1.0,fitfwhm:-1.0,bad:0},noptfiles)
for i=0,noptfiles-1 do begin
  chip = long( (strsplit(file_basename(optfiles[i],'.opt'),'_',/extract))[1] )
  chstr[i].dir = dir
  chstr[i].expnum = base
  chstr[i].chip = chip
  READCOL,optfiles[i],key,equal,value,format='A,A,F',/silent
  fwind = where(key eq 'FW',nfwind)
  if nfwind gt 0 then begin
    fwhm = value[fwind[0]]
    chstr[i].fwhm = fwhm 
  endif else fwhm=-1.0
  if keyword_set(verbose) then print,strtrim(i+1,2),' ',file_basename(optfiles[i]),' ',fwhm
endfor  ; opt file loop

med_fwhm = median(chstr.fwhm)
sig_fwhm = mad(chstr.fwhm)

; Fit linear plane to the FWHM value wrt RA/DEC
restore,'/data/smash/cp/red/photred/decam_chip_xyoff.dat'
offstr = replicate({chip:0L,xoff:0.0,yoff:0.0},n_elements(chip))
offstr.chip = chip
offstr.xoff = xoff
offstr.yoff = yoff
MATCH,chstr.chip,offstr.chip,ind1,ind2,/sort
initpars = dblarr(3)
initpars[0] = med_fwhm
parinfo = REPLICATE({limited:[0,0],limits:[0.0,0.0],fixed:0},3)
;  we can change the order by fixing certain parameters
; 0-constant, 1-x, 2-x^2, 3-x^3
; 4-x*y, 5-x^2*y, 6-x*y^2, 7-x^2*y^2
; 8-y, 9-y^2, 10-y^3
pars = MPFIT2DFUN('func_poly2d',offstr[ind2].xoff,offstr[ind2].yoff,chstr[ind1].fwhm,
                  chstr[ind1].fwhm*0+0.1,initpars,status=status,dof=dof,$
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
  rchisq = chisq/dof
endelse
chstr[ind1].fitfwhm = fitfwhm

; Find outliers
diff_fwhm = chstr.fwhm-chstr.fitfwhm
sig_fitfwhm = mad(diff_fwhm)
bd_fwhm = where(abs(diff_fwhm)/sig_fitfwhm gt nsigthresh or abs(diff_fwhm) gt nabsthresh,nbd_fwhm)
if keyword_set(verbose) and
if nbd_fwhm gt 0 then begin
  if keyword_set(verbose) then print,' bad chips. ',strjoin(chstr[bd_fwhm].chip,' ')
  bdchstr = chstr[bd_fwhm]
  writecol,-1,bdchstr.chip,bdchstr.fwhm,bdchstr.fwhm-med_fwhm,bdchstr.fwhm-bdchstr.fitfwhm,$
              abs(bdchstr.fwhm-bdchstr.med_fwhm)/sig_fitfwhm,fmt='(I4,4F10.3)'
  chstr[bd_fwhm].bad = 1
endif

; Check if there are any problems
if keyword_set(pl) then begin
  plot,chstr.chip,chstr.fwhm,ps=1,xr=[0,63],xs=1,tit=expnum[0]
  oplot,chstr.chip,chstr.fitfwhm,ps=4,co=80
  oplot,[0,63],[0,0]+med_fwhm,linestyle=2,co=80
  if nbd_fwhm gt 0 then oplot,[chstr[bd_fwhm].chip],[chstr[bd_fwhm].fwhm],ps=4,co=250,sym=2,thick=3
endif

if keyword_set(stp) then stop

end
