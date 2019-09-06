pro smashred_calibrate_healpix_cmdplot,pix
  
; Defaults
if n_elements(version) eq 0 then version='v6'
if n_elements(reduxdir) eq 0 then reduxdir=SMASHRED_ROOTDIR()+'cp/red/photred/'
outputdir = reduxdir+'catalogs/final/'
if n_elements(version) gt 0 then outputdir+=version+'/'

allobj = mrdfits(outputdir+strtrim(pix,2)+'_combined_allobj.fits.gz',1)
  
;; Make a CMD figure
setdisp
!p.font = 0
plotdir = outputdir+'/plots/'
if file_test(plotdir,/directory) eq 0 then file_mkdir,plotdir
file = plotdir+strtrim(pix,2)+'_cmd'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=9.5
undefine,dum,im
hess,allobj.g-allobj.i,allobj.g,dum,im,dx=0.02,dy=0.05,xr=[-1,3],yr=[26,13],/noplot,xarr=xarr,yarr=yarr
displayc,im,xarr,yarr,/yflip,xtit='g-i',ytit='g',tit=strtrim(pix,2),charsize=1.2,/log
ps_close
ps2png,file+'.eps',/eps

;stop

end
