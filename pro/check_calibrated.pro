pro check_calibrated,version

dir = smashred_rootdir()+'cp/red/photred/catalogs/'
;dir = '/data/smash/cp/red/photred/catalogs/'
undefine
if n_elements(version) eq 0 then version='v5'
files = file_search(dir+'final/'+version+'/*_combined_allobj.fits.gz',count=nfiles)
;files1 = file_search(dir+'inst/comb/*_combined_final_roughcal.fits',count=nfiles1)
;if nfiles1 gt 0 then push,files,files1
;files2 = file_search(dir+'inst/comb/*_combined_final_roughcal.fits.gz',count=nfiles2)
;if nfiles2 gt 0 then push,files,files2
;nfiles = n_elements(files)

smash = importascii('~/projects/SMASHRED/data/smash_fields_final.txt',/header)

print,strtrim(nfiles,2),' files to process'

;stop
calarr = lonarr(nfiles)
ucalibfracarr = fltarr(nfiles)
gcalibfracarr = fltarr(nfiles)
rcalibfracarr = fltarr(nfiles)
icalibfracarr = fltarr(nfiles)
zcalibfracarr = fltarr(nfiles)
info = replicate({field:'',ra:0.0d,dec:0.0d0,mslon:0.0d0,mslat:0.0d0,ucalibfrac:0.0,nuchips:0L,nucalib:0L,ucalibflag:0,$
                  gcalibfrac:0.0,ngchips:0L,ngcalib:0L,gcalibflag:0,$
                  rcalibfrac:0.0,nrchips:0L,nrcalib:0L,rcalibflag:0,icalibfrac:0.0,nichips:0L,nicalib:0L,$
                  icalibflag:0,zcalibfrac:0.0,nzchips:0L,nzcalib:0L,zcalibflag:0,calflag:0L},nfiles)


for i=0,nfiles-1 do begin
;for i=38,nfiles-1 do begin

  ;print,strtrim(i+1,2),' ',files[i]
  dum = strsplit(file_basename(files[i]),'_',/extract)
  field = first_el(dum[0])
  ;field = file_basename(files[i],'_combined_final_roughcal.fits')

  match,'Field'+strtrim(smash.num,2),field,ind1,ind2,/sort
  info[i].ra = smash[ind1[0]].radeg
  info[i].dec = smash[ind1[0]].dedeg
  info[i].mslon = smash[ind1[0]].mslon
  info[i].mslat = smash[ind1[0]].mslat

  ; Load the CHIPS structure as well
  chipfile = file_dirname(files[i])+'/'+field+'_combined_chips.fits.gz'
  chipstr = mrdfits(chipfile,1,/silent)

  ; Figuring out what data are calibrated
  orange = 200
  ;green = fsc_color('forest green',1)
  green = 1
  red = 250
  ucolor=red & gcolor=red & rcolor=red & icolor=red & zcolor=red
  uchip = where(chipstr.filter eq 'u',nuchips)
  ucalib = where(chipstr.filter eq 'u' and chipstr.calibrated eq 1,nucalib)
  if nuchips gt 0 then ucalibfrac = nucalib/float(nuchips) else ucalibfrac=1.0
  if ucalibfrac gt 0 then ucolor=orange
  if ucalibfrac gt 0.70 then ucolor = green
  if nuchips gt 0 then ucalibflag=chipstr[ucalib[0]].zpcalibflag else ucalibflag=0
  gchip = where(chipstr.filter eq 'g',ngchips)
  gcalib = where(chipstr.filter eq 'g' and chipstr.calibrated eq 1,ngcalib)
  if ngchips gt 0 then gcalibfrac = ngcalib/float(ngchips) else gcalibfrac=1.0
  if gcalibfrac gt 0 then gcolor=orange
  if gcalibfrac gt 0.70 then gcolor = green
  if ngchips gt 0 then gcalibflag=chipstr[gcalib[0]].zpcalibflag else gcalibflag=0
  rchip = where(chipstr.filter eq 'r',nrchips)
  rcalib = where(chipstr.filter eq 'r' and chipstr.calibrated eq 1,nrcalib)
  if nrchips gt 0 then rcalibfrac = nrcalib/float(nrchips) else rcalibfrac=1.0
  if rcalibfrac gt 0 then rcolor=orange
  if rcalibfrac gt 0.70 then rcolor = green
  if nrchips gt 0 then rcalibflag=chipstr[rcalib[0]].zpcalibflag else rcalibflag=0
  ichip = where(chipstr.filter eq 'i',nichips)
  icalib = where(chipstr.filter eq 'i' and chipstr.calibrated eq 1,nicalib)
  if nichips gt 0 then icalibfrac = nicalib/float(nichips) else icalibfrac=1.0
  if icalibfrac gt 0 then icolor=orange
  if icalibfrac gt 0.70 then icolor = green
  if nichips gt 0 then icalibflag=chipstr[icalib[0]].zpcalibflag else icalibflag=0
  zchip = where(chipstr.filter eq 'z',nzchips)
  zcalib = where(chipstr.filter eq 'z' and chipstr.calibrated eq 1,nzcalib)
  if nzchips gt 0 then zcalibfrac = nzcalib/float(nzchips) else zcalibfrac=1.0
  if zcalibfrac gt 0 then zcolor=orange
  if zcalibfrac gt 0.70 then zcolor = green
  if nzchips gt 0 then zcalibflag=chipstr[zcalib[0]].zpcalibflag else zcalibflag=0

  if ucalibfrac gt 0 or gcalibfrac gt 0 or rcalibfrac gt 0 or icalibfrac gt 0 or zcalibfrac gt 0 then calarr[i]=1      ; partial
  if ucalibfrac eq 1 and gcalibfrac eq 1 and rcalibfrac eq 1 and icalibfrac eq 1 and zcalibfrac eq 1 then calarr[i]=2  ; full

  info[i].field = field
  info[i].ucalibfrac = ucalibfrac
  info[i].nuchips = nuchips
  info[i].nucalib = nucalib
  info[i].ucalibflag = ucalibflag
  info[i].gcalibfrac = gcalibfrac
  info[i].ngchips = ngchips
  info[i].ngcalib = ngcalib
  info[i].gcalibflag = gcalibflag
  info[i].rcalibfrac = rcalibfrac
  info[i].nrchips = nrchips
  info[i].nrcalib = nrcalib
  info[i].rcalibflag = rcalibflag
  info[i].icalibfrac = icalibfrac
  info[i].nichips = nichips
  info[i].nicalib = nicalib
  info[i].icalibflag = icalibflag
  info[i].zcalibfrac = zcalibfrac
  info[i].nzchips = nzchips
  info[i].nzcalib = nzcalib
  info[i].zcalibflag = zcalibflag
  info[i].calflag = calarr[i]  

  ;  xyouts,xr[1]-0.25*range(xr),yr[1]+0.15*range(yr),'u: '+strtrim(nucalib,2)+'/'+strtrim(nuchips,2)+' calib',align=0,charsize=1.4,color=ucolor
  ;  xyouts,xr[1]-0.25*range(xr),yr[1]+0.18*range(yr),'g: '+strtrim(ngcalib,2)+'/'+strtrim(ngchips,2)+' calib',align=0,charsize=1.4,color=gcolor
  ;  xyouts,xr[1]-0.25*range(xr),yr[1]+0.21*range(yr),'r: '+strtrim(nrcalib,2)+'/'+strtrim(nrchips,2)+' calib',align=0,charsize=1.4,color=rcolor
  ;  xyouts,xr[1]-0.25*range(xr),yr[1]+0.24*range(yr),'i: '+strtrim(nicalib,2)+'/'+strtrim(nichips,2)+' calib',align=0,charsize=1.4,color=icolor
  ;  xyouts,xr[1]-0.25*range(xr),yr[1]+0.27*range(yr),'z: '+strtrim(nzcalib,2)+'/'+strtrim(nzchips,2)+' calib',align=0,charsize=1.4,color=zcolor

  print,i+1,'  '+field+'  ',calarr[i],'  '+strtrim(nucalib,2)+'/'+strtrim(nuchips,2),'  '+strtrim(ngcalib,2)+'/'+strtrim(ngchips,2),$
        '  '+strtrim(nrcalib,2)+'/'+strtrim(nrchips,2),'  '+strtrim(nicalib,2)+'/'+strtrim(nichips,2),'  '+strtrim(nzcalib,2)+'/'+strtrim(nzchips,2)

  BOMB:

endfor

setdisp
plotc,info.mslon,info.mslat,info.calflag,ps=1,/xflip
g=where(info.gcalibfrac eq 1 and info.icalibfrac eq 1,ng)
oplot,info[g].mslon,info[g].mslat,ps=6,sym=1.5,co=200

;mwrfits,info,dir+'final/'+version+'/check_calibrated_'+version+'.fits',/create

stop

end
