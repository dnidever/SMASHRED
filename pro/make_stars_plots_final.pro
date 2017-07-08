pro make_stars_plots_final,version,sversion,showcalib=showcalib,redo=redo

dir = smashred_rootdir()+'cp/red/photred/catalogs/'
;dir = '/data/smash/cp/red/photred/catalogs/'
undefine
if n_elements(version) eq 0 then begin
  print,'Need to input the version, i.e. "v1"'
  return
endif
if n_elements(sversion) eq 0 then sversion='1'
files = file_search(dir+'final/'+version+'/stars'+sversion+'/*_allobj_stars.fits.gz',count=nfiles)
;files = file_search(dir+'final/'+version+'/*_combined_allobj.fits.gz',count=nfiles)
;files1 = file_search(dir+'inst/comb/*_combined_final_roughcal.fits',count=nfiles1)
;if nfiles1 gt 0 then push,files,files1
;files2 = file_search(dir+'inst/comb/*_combined_final_roughcal.fits.gz',count=nfiles2)
;if nfiles2 gt 0 then push,files,files2
;nfiles = n_elements(files)

smash = importascii('~/projects/SMASHRED/data/smash_fields_final.txt',/header)

print,strtrim(nfiles,2),' files to process'

; Make the directory if it doesn't exist
catdir = dir+'final/'+version+'/'
outdir = dir+'plots/final/'+version+'/'
if file_test(outdir,/directory) eq 0 then file_mkdir,outdir

; Defaults
if n_elements(showcalib) eq 0 then showcalib=1

;stop

for i=0,nfiles-1 do begin
;for i=38,nfiles-1 do begin

  print,strtrim(i+1,2),' ',files[i]
  dum = strsplit(file_basename(files[i]),'_',/extract)
  field = reform(dum[0])
  ;field = file_basename(files[i],'_combined_final_roughcal.fits')

  ; Check if plots already exist
  file1 = outdir+field+'_hess_final'+version+'_stars'+sversion
  file2 = outdir+field+'_hess_final'+version+'_stars'+sversion+'_dered'
  if file_test(file1+'.jpg') eq 1 and file_test(file2+'.jpg') eq 1 and not keyword_set(redo) then goto,BOMB

  str = mrdfits(files[i],1)

  ; Load the CHIPS structure as well
  ;chipfile = file_dirname(files[i])+'/'+field+'_combined_chips.fits.gz'
  chipfile = catdir+'/'+field+'_combined_chips.fits.gz'
  chipstr = mrdfits(chipfile,1)

  ; Figuring out how the zero-point was set
  orange = 200
  ;green = fsc_color('forest green',1)
  green = 1
  lgreen = 2
  red = 250
  blue = 80
  colarr = [0,green,lgreen,orange,blue]
  caltxt = ['none','cal','olap','0.9m','gaia']
  ucaltxt = 'none'
  uchip = where(chipstr.filter eq 'u',nuchips)
  if nuchips gt 0 then begin
    uzptermflag = long(median([chipstr[uchip].zpcalibflag]))
    ucolor = colarr[uzptermflag]
    ucaltxt = caltxt[uzptermflag]
  endif else begin
    uzptermflag = 0
    ucolor = 0
    ucaltxt= 'none'
  endelse
  gchip = where(chipstr.filter eq 'g',ngchips)
  if ngchips gt 0 then begin
    gzptermflag = long(median([chipstr[gchip].zpcalibflag]))
    gcolor = colarr[gzptermflag]
    gcaltxt = caltxt[gzptermflag]
  endif else begin
    gzptermflag = 0
    gcolor = 0
    gcaltxt = 'none'
  endelse
  rchip = where(chipstr.filter eq 'r',nrchips)
  if nrchips gt 0 then begin
    rzptermflag = long(median([chipstr[rchip].zpcalibflag]))
    rcolor = colarr[rzptermflag]
    rcaltxt = caltxt[rzptermflag]
  endif else begin
    rzptermflag = 0
    rcolor = 0
    rcaltxt = 'none'
  endelse
  ichip = where(chipstr.filter eq 'i',nichips)
  if nichips gt 0 then begin
    izptermflag = long(median([chipstr[ichip].zpcalibflag]))
    icolor = colarr[izptermflag]
    icaltxt = caltxt[izptermflag]
  endif else begin
    izptermflag = 0
    icolor = 0
    icaltxt = 'none'
  endelse
  zchip = where(chipstr.filter eq 'z',nzchips)
  if nzchips gt 0 then begin
    zzptermflag = long(median([chipstr[zchip].zpcalibflag]))
    zcolor = colarr[zzptermflag]
    zcaltxt = caltxt[zzptermflag]
  endif else begin
    zzptermflag = 0
    zcolor = 0
    zcaltxt = 'none'
  endelse

  ;; Figuring out what data are calibrated
  ;orange = 200
  ;;green = fsc_color('forest green',1)
  ;green = 1
  ;red = 250
  ;ucolor=red & gcolor=red & rcolor=red & icolor=red & zcolor=red
  ;uchip = where(chipstr.filter eq 'u',nuchips)
  ;ucalib = where(chipstr.filter eq 'u' and chipstr.calibrated eq 1,nucalib)
  ;ucalibfrac = nucalib/float(nuchips)
  ;if ucalibfrac gt 0 then ucolor=orange
  ;if ucalibfrac gt 0.70 then ucolor = green
  ;gchip = where(chipstr.filter eq 'g',ngchips)
  ;gcalib = where(chipstr.filter eq 'g' and chipstr.calibrated eq 1,ngcalib)
  ;gcalibfrac = ngcalib/float(ngchips)
  ;if gcalibfrac gt 0 then gcolor=orange
  ;if gcalibfrac gt 0.70 then gcolor = green
  ;rchip = where(chipstr.filter eq 'r',nrchips)
  ;rcalib = where(chipstr.filter eq 'r' and chipstr.calibrated eq 1,nrcalib)
  ;rcalibfrac = nrcalib/float(nrchips)
  ;if rcalibfrac gt 0 then rcolor=orange
  ;if rcalibfrac gt 0.70 then rcolor = green
  ;ichip = where(chipstr.filter eq 'i',nichips)
  ;icalib = where(chipstr.filter eq 'i' and chipstr.calibrated eq 1,nicalib)
  ;icalibfrac = nicalib/float(nichips)
  ;if icalibfrac gt 0 then icolor=orange
  ;if icalibfrac gt 0.70 then icolor = green
  ;zchip = where(chipstr.filter eq 'z',nzchips)
  ;zcalib = where(chipstr.filter eq 'z' and chipstr.calibrated eq 1,nzcalib)
  ;zcalibfrac = nzcalib/float(nzchips)
  ;if zcalibfrac gt 0 then zcolor=orange
  ;if zcalibfrac gt 0.70 then zcolor = green

  ; morphology cuts
  ;gd = where(str.chi lt 3 and abs(str.sharp) lt 1 and str.prob gt 0.8,ngd)
  ;gd = where(str.chi lt 3 and abs(str.sharp) lt 0.3,ngd)
  gd = where(str.chi lt 3 and abs(str.sharp) lt 1.0,ngd)
  ; the PROB values for combined_final_roughcal seem HORRIBLE
  ; I think something's wrong
  str = str[gd]

  tags = tag_names(str)
  uind = where(tags eq 'U',nuind)
  gind = where(tags eq 'G',ngind)
  rind = where(tags eq 'R',nrind)
  iind = where(tags eq 'I',niind)
  zind = where(tags eq 'Z',nzind)

  if nuind gt 0 then dum=where(str.(uind[0]) lt 50,nu) else nu=0
  if ngind gt 0 then dum=where(str.(gind[0]) lt 50,ng) else ng=0
  if nrind gt 0 then dum=where(str.(rind[0]) lt 50,nr) else nr=0
  if niind gt 0 then dum=where(str.(iind[0]) lt 50,ni) else ni=0
  if nzind gt 0 then dum=where(str.(zind[0]) lt 50,nz) else nz=0
  ndet = [nu,ng,nr,ni,nz]
  magind = [uind,gind,rind,iind,zind]
  magnames = ['u','g','r','i','z']
  ;magext = [4.239, 3.303, 2.285, 1.698, 1.263]
  magext = [5.155, 3.793, 2.751, 2.086, 1.479]  ; from SFD98

  ;nbands = (nu gt 0)+(ng gt 0)+(nr gt 0)+(ni gt 0)+(nz gt 0)
  ;if nbands eq 1 then begin
  if ng eq 0 or ni eq 0 then begin
    ;print,'Only 1 band.  Need at least 2 to make a CMD!'
    print,'Need g and i'
    goto, BOMB
  endif

  ;get_stellar_locus,str,locstr
  ;cut_stellar_locus,str,locstr,ind
  
  ; Extinction coefficients
  ;u  4.239
  ;g  3.303
  ;r  2.285
  ;i  1.698
  ;z  1.263

  ;undefine,bind
  ;if ni gt 0 and ng gt 0 then begin
    bind = gind
    rind = iind
    mind = gind
    xtit = 'g-i'
    ytit = 'g'
    xtit0 = '(g-i)!d0!n'
    ytit0 = 'g!d0!n'
    xr = [-1,4]
    col = str.g-str.i
    mag = str.g
    ;col_locus = str[ind].gmag-str[ind].imag
    ;mag_locus = str[ind].gmag
    col0 = (str.g-str.ebv*magext[1])-(str.i-str.ebv*magext[3])
    mag0 = str.g-str.ebv*magext[1]
    ;col0 = (str.g-str.ebv*3.303)-(str.i-str.ebv*1.698)
    ;mag0 = str.g-str.ebv*3.303
    ;col_locus0 = (str[ind].gmag-str[ind].ebv*3.303)-(str[ind].imag-str[ind].ebv*1.698)
    ;mag_locus0 = str[ind].gmag-str[ind].ebv*3.303
  ;endif
  ;;;;if ni eq 0 and ng gt 0 and nu gt 0 then begin
  ;;;;  bind = uind
  ;  rind = gind
  ;  mind = gind
  ;  xtit = 'u-g'
  ;  ytit = 'g'
  ;  xtit0 = '(u-g)!d0!n'
  ;  ytit0 = 'g!d0!n'
  ;  xr = [-1,4]
  ;  col = str.umag-str.gmag
  ;  mag = str.gmag
  ;  ;col_locus = str[ind].umag-str[ind].gmag
  ;  ;mag_locus = str[ind].gmag
  ;  col0 = (str.umag-str.ebv*4.239)-(str.gmag-str.ebv*3.303)
  ;  mag0 = str.gmag-str.ebv*3.303
  ;  ;col_locus0 = (str[ind].umag-str[ind].ebv*4.239)-(str[ind].gmag-str[ind].ebv*3.303)
  ;  ;mag_locus0 = str[ind].gmag-str[ind].ebv*3.303
  ;endif
  ;if ng gt 0 and nr gt 0 and ni eq 0 and nu eq 0 then begin
  ;  bind = gind
  ;  rind = rind
  ;  mind = gind
  ;  xtit = 'g-r'
  ;  ytit = 'g'
  ;  xtit0 = '(g-r)!d0!n'
  ;  ytit0 = 'g!d0!n'
  ;  xr = [-1,4]
  ;  col = str.gmag-str.rmag
  ;  mag = str.gmag
  ;  ;col_locus = str[ind].umag-str[ind].gmag
  ;  ;mag_locus = str[ind].gmag
  ;  col0 = (str.gmag-str.ebv*3.303)-(str.rmag-str.ebv*2.285)
  ;  mag0 = str.gmag-str.ebv*3.303
  ;  ;col_locus0 = (str[ind].umag-str[ind].ebv*4.239)-(str[ind].gmag-str[ind].ebv*3.303)
  ;  ;mag_locus0 = str[ind].gmag-str[ind].ebv*3.303
  ;endif
  ;if nu eq 0 and ng eq 0 and ni eq 0 then begin
  ;  bind = rind
  ;  rind = zind
  ;  mind = rind
  ;  xtit = 'r-z'
  ;  ytit = 'r'
  ;  xtit0 = '(r-z)!d0!n'
  ;  ytit0 = 'r!d0!n'
  ;  xr = [-1,2]
  ;  col = str.rmag-str.zmag
  ;  mag = str.rmag
  ;  ;col_locus = str[ind].rmag-str[ind].zmag
  ;  ;mag_locus = str[ind].rmag
  ;  col0 = (str.rmag-str.ebv*2.285)-(str.zmag-str.ebv*1.263)
  ;  mag0 = str.rmag-str.ebv*2.285
  ;  ;col_locus0 = (str[ind].rmag-str[ind].ebv*2.285)-(str[ind].zmag-str[ind].ebv*1.263)
  ;  ;mag_locus0 = str[ind].rmag-str[ind].ebv*2.285
  ;endif
  ;
  ;; No matches, take the reddest and bluest bands
  ;if n_elements(bind) eq 0 then begin
  ;  g = where(ndet gt 0,ng)
  ;  rmagind = max(g)
  ;  bmagind = min(g)
  ;  mmagind = bmagind
  ;  rind = magind[rmagind]
  ;  bind = magind[bmagind]
  ;  mind = bind
  ;  xtit = magnames[bmagind]+'-'+magnames[rmagind]
  ;  ytit = magnames[bmagind]
  ;  xtit0 = '('+xtit+')!d0!n'
  ;  ytit0 = ytit+'!d0!n'
  ;  xr = [-1,4]
  ;  col = str.(bind)-str.(rind)
  ;  mag = str.(mind)
  ;  ;col_locus = str[ind].gmag-str[ind].imag
  ;  ;mag_locus = str[ind].gmag
  ;  col0 = (str.(bind)-str.ebv*magext[bmagind])-(str.(rind)-str.ebv*magext[rmagind])
  ;  mag0 = str.(mind)-str.ebv*magext[mmagind]
  ;endif


;  ;setdisp
;  !p.font = 0
;  file = '../plots/'+field+'_cmd'
;  if file_test(file+'.eps') eq 0 or keyword_set(redo) then begin
;    ps_open,file,/color,thick=4,/encap
;    ;device,/inches,xsize=8.5,ysize=10
;    device,/inches,xsize=8.5,ysize=11
;    ;xr = [-1,4]
;    ;yr = [25,12]
;    yr = [25,11.5]
;    plot,str.(bind)-str.(rind),str.(mind),ps=3,xr=xr,yr=yr,xs=1,ys=1,xtit=xtit,ytit=ytit,$
;         charsize=1.5,tit=field+' Instrumental PSF Photometry'
;    ra = median(str.ra)
;    dec = median(str.dec)
;    xyouts,xr[0]+0.05*range(xr),yr[1]+0.05*range(yr),'('+textoidl('\alpha')+','+textoidl('\delta')+$
;           ')=('+stringize(ra,ndec=1)+','+stringize(dec,ndec=1)+')',align=0,charsize=1.4
;    glactc,ra,dec,2000.,glon,glat,1,/deg
;    xyouts,xr[0]+0.05*range(xr),yr[1]+0.09*range(yr),'(l,b)=('+stringize(glon,ndec=1)+','+stringize(glat,ndec=1)+')',align=0,charsize=1.4
;    gal2mag,glon,glat,mlon,mlat
;    xyouts,xr[0]+0.05*range(xr),yr[1]+0.13*range(yr),'(L!dMS!n,B!dMS!n)=('+stringize(mlon,ndec=1)+','+stringize(mlat,ndec=1)+')',align=0,charsize=1.4
;    cel2lmc,ra,dec,palmc,radlmc  
;    xyouts,xr[0]+0.55*range(xr),yr[1]+0.05*range(yr),'('+textoidl('\rho')+'!dLMC!n,PA!dLMC!n)=('+stringize(radlmc[0],ndec=1)+','+$
;           stringize(palmc[0],ndec=1)+')',align=0,charsize=1.4
;    cel2smc,ra,dec,pasmc,radsmc  
;    xyouts,xr[0]+0.55*range(xr),yr[1]+0.09*range(yr),'('+textoidl('\rho')+'!dSMC!n,PA!dSMC!n)=('+stringize(radsmc[0],ndec=1)+','+$
;            stringize(pasmc[0],ndec=1)+')',align=0,charsize=1.4
;    ps_close
;    ps2jpg,file+'.eps',/eps
;  endif

  ;setdisp
  !p.font = 0
  file = outdir+field+'_hess_final'+version+'_stars'+sversion
  if file_test(file+'.eps') eq 0 or file_test(file+'.jpg') eq 0 or keyword_set(redo) then begin
  ;if file_test(file+'.eps') eq 0 or keyword_set(redo) then begin
    ps_open,file,/color,thick=4,/encap
    ;device,/inches,xsize=8.5,ysize=10
    device,/inches,xsize=8.5,ysize=11
    ;xr = [-1,4]
    ;yr = [25,12]
    yr = [26,14.0]
    loadcol,3
    black = fsc_color('black',0)
    ;hess,str.(bind)-str.(rind),str.(mind),dx=0.02,dy=0.05,xr=xr,yr=yr,xtit=xtit,ytit=ytit,$
    ;     charsize=1.5,/log,tit=field+' Deep Calibrated PSF Photometry'
    undefine,dum,im
    hess,str.(bind)-str.(rind),str.(mind),dum,im,dx=0.02,dy=0.05,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
    loadcol,3
    black = fsc_color('black',0)
    display,im,xarr,yarr,xtit=xtit,ytit=ytit,/yflip,position=[0.12,0.10,0.98,0.94],$
         charsize=2.2,/log,tit=field+' stars'
    ra = median(str.ra)
    dec = median(str.dec)
    xyouts,xr[0]+0.05*range(xr),yr[1]+0.05*range(yr),'('+textoidl('\alpha')+','+textoidl('\delta')+$
           ' )=('+stringize(ra,ndec=1)+','+stringize(dec,ndec=1)+')',align=0,charsize=1.4
    glactc,ra,dec,2000.,glon,glat,1,/deg
    xyouts,xr[0]+0.05*range(xr),yr[1]+0.09*range(yr),'(l,b)=('+stringize(glon,ndec=1)+','+stringize(glat,ndec=1)+')',align=0,charsize=1.4
    gal2mag,glon,glat,mlon,mlat
    xyouts,xr[0]+0.05*range(xr),yr[1]+0.13*range(yr),'(L!dMS!n,B!dMS!n)=('+stringize(mlon,ndec=1)+','+stringize(mlat,ndec=1)+')',align=0,charsize=1.4
    cel2lmc,ra,dec,palmc,radlmc  
    xyouts,xr[0]+0.55*range(xr),yr[1]+0.05*range(yr),'('+textoidl('\rho')+'!dLMC!n,PA!dLMC!n)=('+stringize(radlmc[0],ndec=1)+','+$
           stringize(palmc[0],ndec=1)+')',align=0,charsize=1.4
    cel2smc,ra,dec,pasmc,radsmc  
    xyouts,xr[0]+0.55*range(xr),yr[1]+0.09*range(yr),'('+textoidl('\rho')+'!dSMC!n,PA!dSMC!n)=('+stringize(radsmc[0],ndec=1)+','+$
           stringize(pasmc[0],ndec=1)+')',align=0,charsize=1.4
    if keyword_set(showcalib) then begin
      loadct,39,/silent
      orange = 200
      green = fsc_color('forest green',1)
      lgreen = fsc_color('light sea green',2)  ; green yellow
      red = 250
      xyouts,xr[1]-0.22*range(xr),yr[1]+0.15*range(yr),'u: '+strtrim(nuchips,2)+' ('+ucaltxt+')',align=0,charsize=1.4,color=ucolor
      xyouts,xr[1]-0.22*range(xr),yr[1]+0.18*range(yr),'g: '+strtrim(ngchips,2)+' ('+gcaltxt+')',align=0,charsize=1.4,color=gcolor
      xyouts,xr[1]-0.22*range(xr),yr[1]+0.21*range(yr),'r: '+strtrim(nrchips,2)+' ('+rcaltxt+')',align=0,charsize=1.4,color=rcolor
      xyouts,xr[1]-0.22*range(xr),yr[1]+0.24*range(yr),'i: '+strtrim(nichips,2)+' ('+icaltxt+')',align=0,charsize=1.4,color=icolor
      xyouts,xr[1]-0.22*range(xr),yr[1]+0.27*range(yr),'z: '+strtrim(nzchips,2)+' ('+zcaltxt+')',align=0,charsize=1.4,color=zcolor
      ;xyouts,xr[1]-0.25*range(xr),yr[1]+0.15*range(yr),'u: '+strtrim(nucalib,2)+'/'+strtrim(nuchips,2)+' calib',align=0,charsize=1.4,color=ucolor
      ;xyouts,xr[1]-0.25*range(xr),yr[1]+0.18*range(yr),'g: '+strtrim(ngcalib,2)+'/'+strtrim(ngchips,2)+' calib',align=0,charsize=1.4,color=gcolor
      ;xyouts,xr[1]-0.25*range(xr),yr[1]+0.21*range(yr),'r: '+strtrim(nrcalib,2)+'/'+strtrim(nrchips,2)+' calib',align=0,charsize=1.4,color=rcolor
      ;xyouts,xr[1]-0.25*range(xr),yr[1]+0.24*range(yr),'i: '+strtrim(nicalib,2)+'/'+strtrim(nichips,2)+' calib',align=0,charsize=1.4,color=icolor
      ;xyouts,xr[1]-0.25*range(xr),yr[1]+0.27*range(yr),'z: '+strtrim(nzcalib,2)+'/'+strtrim(nzchips,2)+' calib',align=0,charsize=1.4,color=zcolor
    endif
    ps_close
    ps2jpg,file+'.eps',/eps
  endif

  ; DEREDDENED CMDs
  ;-----------------
  ;setdisp
 ; !p.font = 0
 ; file = '../plots/'+field+'_cmd_dered'
 ; if file_test(file+'.eps') eq 0 or keyword_set(redo) then begin
 ;   ps_open,file,/color,thick=4,/encap
 ;   ;device,/inches,xsize=8.5,ysize=10
 ;   device,/inches,xsize=8.5,ysize=11
 ;   ;xr = [-1,4]
 ;   ;yr = [25,12]
 ;   yr = [25,11.5]
 ;   plot,col0,mag0,ps=3,xr=xr,yr=yr,xs=1,ys=1,xtit=xtit0,ytit=ytit0,$
 ;        charsize=1.5,tit=field+' Instrumental PSF Photometry (dered)'
 ;   ra = median(str.ra)
 ;   dec = median(str.dec)
 ;   xyouts,xr[0]+0.05*range(xr),yr[1]+0.05*range(yr),'('+textoidl('\alpha')+','+textoidl('\delta')+$
 ;          ')=('+stringize(ra,ndec=1)+','+stringize(dec,ndec=1)+')',align=0,charsize=1.4
 ;   glactc,ra,dec,2000.,glon,glat,1,/deg
 ;   xyouts,xr[0]+0.05*range(xr),yr[1]+0.09*range(yr),'(l,b)=('+stringize(glon,ndec=1)+','+stringize(glat,ndec=1)+')',align=0,charsize=1.4
 ;   gal2mag,glon,glat,mlon,mlat
 ;   xyouts,xr[0]+0.05*range(xr),yr[1]+0.13*range(yr),'(L!dMS!n,B!dMS!n)=('+stringize(mlon,ndec=1)+','+stringize(mlat,ndec=1)+')',align=0,charsize=1.4
 ;   cel2lmc,ra,dec,palmc,radlmc  
 ;   xyouts,xr[0]+0.55*range(xr),yr[1]+0.05*range(yr),'('+textoidl('\rho')+'!dLMC!n,PA!dLMC!n)=('+stringize(radlmc[0],ndec=1)+','+$
 ;          stringize(palmc[0],ndec=1)+')',align=0,charsize=1.4
 ;   cel2smc,ra,dec,pasmc,radsmc  
 ;   xyouts,xr[0]+0.55*range(xr),yr[1]+0.09*range(yr),'('+textoidl('\rho')+'!dSMC!n,PA!dSMC!n)=('+stringize(radsmc[0],ndec=1)+','+$
 ;           stringize(pasmc[0],ndec=1)+')',align=0,charsize=1.4
 ;   ps_close
 ;   ps2jpg,file+'.eps',/eps
 ; endif

  ;setdisp
  !p.font = 0
  file = outdir+field+'_hess_final'+version+'_stars'+sversion+'_dered'
  if file_test(file+'.eps') eq 0 or file_test(file+'.jpg') eq 0 or keyword_set(redo) then begin
  ;if file_test(file+'.eps') eq 0 or keyword_set(redo) then begin
    ps_open,file,/color,thick=4,/encap
    ;device,/inches,xsize=8.5,ysize=10
    device,/inches,xsize=8.5,ysize=11
    ;xr = [-1,4]
    ;yr = [25,12]
    ;yr = [25,11.5]
    yr = [26,14.0]
    loadcol,3
    black = fsc_color('black',0)
    ;hess,col0,mag0,dx=0.02,dy=0.05,xr=xr,yr=yr,xtit=xtit0,ytit=ytit0,$
    ;     charsize=1.5,/log,tit=field+' Instrumental PSF Photometry (dered)'
    hess,col0,mag0,dum,im,dx=0.02,dy=0.05,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
    loadcol,3
    black = fsc_color('black',0)
    display,im,xarr,yarr,xtit=xtit0,ytit=ytit0,/yflip,position=[0.12,0.10,0.98,0.94],$
         charsize=2.2,/log,tit=field+' stars'
    ra = median(str.ra)
    dec = median(str.dec)
    xyouts,xr[0]+0.05*range(xr),yr[1]+0.05*range(yr),'('+textoidl('\alpha')+','+textoidl('\delta')+$
           ' )=('+stringize(ra,ndec=1)+','+stringize(dec,ndec=1)+')',align=0,charsize=1.4
    glactc,ra,dec,2000.,glon,glat,1,/deg
    xyouts,xr[0]+0.05*range(xr),yr[1]+0.09*range(yr),'(l,b)=('+stringize(glon,ndec=1)+','+stringize(glat,ndec=1)+')',align=0,charsize=1.4
    gal2mag,glon,glat,mlon,mlat
    xyouts,xr[0]+0.05*range(xr),yr[1]+0.13*range(yr),'(L!dMS!n,B!dMS!n)=('+stringize(mlon,ndec=1)+','+stringize(mlat,ndec=1)+')',align=0,charsize=1.4
    cel2lmc,ra,dec,palmc,radlmc  
    xyouts,xr[0]+0.55*range(xr),yr[1]+0.05*range(yr),'('+textoidl('\rho')+'!dLMC!n,PA!dLMC!n)=('+stringize(radlmc[0],ndec=1)+','+$
           stringize(palmc[0],ndec=1)+')',align=0,charsize=1.4
    cel2smc,ra,dec,pasmc,radsmc  
    xyouts,xr[0]+0.55*range(xr),yr[1]+0.09*range(yr),'('+textoidl('\rho')+'!dSMC!n,PA!dSMC!n)=('+stringize(radsmc[0],ndec=1)+','+$
           stringize(pasmc[0],ndec=1)+')',align=0,charsize=1.4
    if keyword_set(showcalib) then begin
      loadct,39,/silent
      orange = 200
      green = fsc_color('forest green',1)
      lgreen = fsc_color('light sea green',2)  ; green yellow
      red = 250
      xyouts,xr[1]-0.22*range(xr),yr[1]+0.15*range(yr),'u: '+strtrim(nuchips,2)+' ('+ucaltxt+')',align=0,charsize=1.4,color=ucolor
      xyouts,xr[1]-0.22*range(xr),yr[1]+0.18*range(yr),'g: '+strtrim(ngchips,2)+' ('+gcaltxt+')',align=0,charsize=1.4,color=gcolor
      xyouts,xr[1]-0.22*range(xr),yr[1]+0.21*range(yr),'r: '+strtrim(nrchips,2)+' ('+rcaltxt+')',align=0,charsize=1.4,color=rcolor
      xyouts,xr[1]-0.22*range(xr),yr[1]+0.24*range(yr),'i: '+strtrim(nichips,2)+' ('+icaltxt+')',align=0,charsize=1.4,color=icolor
      xyouts,xr[1]-0.22*range(xr),yr[1]+0.27*range(yr),'z: '+strtrim(nzchips,2)+' ('+zcaltxt+')',align=0,charsize=1.4,color=zcolor
      ;xyouts,xr[1]-0.25*range(xr),yr[1]+0.15*range(yr),'u: '+strtrim(nucalib,2)+'/'+strtrim(nuchips,2)+' calib',align=0,charsize=1.4,color=ucolor
      ;xyouts,xr[1]-0.25*range(xr),yr[1]+0.18*range(yr),'g: '+strtrim(ngcalib,2)+'/'+strtrim(ngchips,2)+' calib',align=0,charsize=1.4,color=gcolor
      ;xyouts,xr[1]-0.25*range(xr),yr[1]+0.21*range(yr),'r: '+strtrim(nrcalib,2)+'/'+strtrim(nrchips,2)+' calib',align=0,charsize=1.4,color=rcolor
      ;xyouts,xr[1]-0.25*range(xr),yr[1]+0.24*range(yr),'i: '+strtrim(nicalib,2)+'/'+strtrim(nichips,2)+' calib',align=0,charsize=1.4,color=icolor
      ;xyouts,xr[1]-0.25*range(xr),yr[1]+0.27*range(yr),'z: '+strtrim(nzcalib,2)+'/'+strtrim(nzchips,2)+' calib',align=0,charsize=1.4,color=zcolor
    endif
    ps_close
    ps2jpg,file+'.eps',/eps
  endif

;  ; Stellar Locus CUTS and DEREDDENED
;  ;----------------------------------
;  ;setdisp
;  !p.font = 0
;  file = '../plots/'+field+'_cmd_locus_dered'
;  if file_test(file+'.eps') eq 0 or keyword_set(redo) then begin
;    ps_open,file,/color,thick=4,/encap
;    ;device,/inches,xsize=8.5,ysize=10
;    device,/inches,xsize=8.5,ysize=11
;    ;xr = [-1,4]
;    ;yr = [25,12]
;    yr = [25,11.5]
;    plot,col_locus0,mag_locus0,ps=3,xr=xr,yr=yr,xs=1,ys=1,xtit=xtit0,ytit=ytit0,$
;         charsize=1.5,tit=field+' Instrumental PSF Photometry (locus, dered)'
;    ra = median(str.ra)
;    dec = median(str.dec)
;    xyouts,xr[0]+0.05*range(xr),yr[1]+0.05*range(yr),'('+textoidl('\alpha')+','+textoidl('\delta')+$
;           ')=('+stringize(ra,ndec=1)+','+stringize(dec,ndec=1)+')',align=0,charsize=1.4
;    glactc,ra,dec,2000.,glon,glat,1,/deg
;    xyouts,xr[0]+0.05*range(xr),yr[1]+0.09*range(yr),'(l,b)=('+stringize(glon,ndec=1)+','+stringize(glat,ndec=1)+')',align=0,charsize=1.4
;    gal2mag,glon,glat,mlon,mlat
;    xyouts,xr[0]+0.05*range(xr),yr[1]+0.13*range(yr),'(L!dMS!n,B!dMS!n)=('+stringize(mlon,ndec=1)+','+stringize(mlat,ndec=1)+')',align=0,charsize=1.4
;    cel2lmc,ra,dec,palmc,radlmc  
;    xyouts,xr[0]+0.55*range(xr),yr[1]+0.05*range(yr),'('+textoidl('\rho')+'!dLMC!n,PA!dLMC!n)=('+stringize(radlmc[0],ndec=1)+','+$
;           stringize(palmc[0],ndec=1)+')',align=0,charsize=1.4
;    cel2smc,ra,dec,pasmc,radsmc  
;    xyouts,xr[0]+0.55*range(xr),yr[1]+0.09*range(yr),'('+textoidl('\rho')+'!dSMC!n,PA!dSMC!n)=('+stringize(radsmc[0],ndec=1)+','+$
;            stringize(pasmc[0],ndec=1)+')',align=0,charsize=1.4
;    ps_close
;    ps2jpg,file+'.eps',/eps
;  endif
;
;  ;setdisp
;  !p.font = 0
;  file = '../plots/'+field+'_hess_locus_dered'
;  if file_test(file+'.eps') eq 0 or keyword_set(redo) then begin
;    ps_open,file,/color,thick=4,/encap
;    device,/inches,xsize=8.5,ysize=10
;    device,/inches,xsize=8.5,ysize=11
;    ;xr = [-1,4]
;    ;yr = [25,12]
;    yr = [25,11.5]
;    loadcol,3
;    black = fsc_color('black',0)
;    hess,col_locus0,mag_locus0,dx=0.02,dy=0.05,xr=xr,yr=yr,xtit=xtit0,ytit=ytit0,$
;         charsize=1.5,/log,tit=field+' Instrumental PSF Photometry (locus, dered)'
;    ra = median(str.ra)
;    dec = median(str.dec)
;    xyouts,xr[0]+0.05*range(xr),yr[1]+0.05*range(yr),'('+textoidl('\alpha')+','+textoidl('\delta')+$
;           ')=('+stringize(ra,ndec=1)+','+stringize(dec,ndec=1)+')',align=0,charsize=1.4
;    glactc,ra,dec,2000.,glon,glat,1,/deg
;    xyouts,xr[0]+0.05*range(xr),yr[1]+0.09*range(yr),'(l,b)=('+stringize(glon,ndec=1)+','+stringize(glat,ndec=1)+')',align=0,charsize=1.4
;    gal2mag,glon,glat,mlon,mlat
;    xyouts,xr[0]+0.05*range(xr),yr[1]+0.13*range(yr),'(L!dMS!n,B!dMS!n)=('+stringize(mlon,ndec=1)+','+stringize(mlat,ndec=1)+')',align=0,charsize=1.4
;    cel2lmc,ra,dec,palmc,radlmc  
;    xyouts,xr[0]+0.55*range(xr),yr[1]+0.05*range(yr),'('+textoidl('\rho')+'!dLMC!n,PA!dLMC!n)=('+stringize(radlmc[0],ndec=1)+','+$
;           stringize(palmc[0],ndec=1)+')',align=0,charsize=1.4
;    cel2smc,ra,dec,pasmc,radsmc  
;    xyouts,xr[0]+0.55*range(xr),yr[1]+0.09*range(yr),'('+textoidl('\rho')+'!dSMC!n,PA!dSMC!n)=('+stringize(radsmc[0],ndec=1)+','+$
;           stringize(pasmc[0],ndec=1)+')',align=0,charsize=1.4
;    ps_close
;    ps2jpg,file+'.eps',/eps
;  endif

  ;stop

  BOMB:

endfor


stop

end
