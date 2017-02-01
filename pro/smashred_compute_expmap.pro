;+
;
; SMASHRED_COMPUTE_EXPMAP
;
; Compute the exposure map for this field in each band
;
; INPUTS:
;  field    The SMASH field name, e.g. "Field100".
;  chstr    The chip structure
;  =outputdir  The output directory for the exposure map.
;
; OUTPUTS:
;  The exposure map is saved to OUTPUTDIR with the name
;  FIELD_combined_expmap.fits.gz'.
;
; USAGE:
;  IDL>smashred_compute_exmap,'Field100',chstr
;
; By D.Nidever Sep. 2016
;-

pro smashred_compute_expmap,field,chstr,redo=redo,outputdir=outputdir

; Not enough inputs
if n_elements(field) eq 0 or n_elements(chstr) eq 0 then begin
  print,'Syntax - smashred_compute_expmap,field,chstr,redo=redo,outputdir=outputdir'
  return
endif

reduxdir = SMASHRED_ROOTDIR()+'cp/red/photred/'
;reduxdir = '/data/smash/cp/red/photred/'
if n_elements(outputdir) eq 0 then outputdir=reduxdir+'catalogs/fina/'
outfile = outputdir+field+'_combined_expmap.fits'

; Output file already exists
if (file_test(outfile) eq 1 or file_test(outfile+'.gz') eq 1) and not keyword_set(redo) then begin
  if file_test(outfile) eq 1 then print,outfile+' EXISTS and /redo NOT set' else $
     print,outfile+'.gz EXISTS and /redo NOT set'
  return
endif

print,'Computing Exposure Map for Field=',field

; Figure out field center and RA/DEC ranges
if range(chstr.vertices_ra) gt 180 then begin
  vertices_ra = chstr.vertices_ra
  over = where(vertices_ra gt 180,nover,comp=under,ncomp=nunder)
  if nover gt 0 then vertices_ra[over]-=360
  rar = minmax(vertices_ra)
  cenra = mean(rar)
endif else begin
  rar = minmax(chstr.vertices_ra)
  cenra = mean(rar)
endelse
decr = minmax(chstr.vertices_dec)
cendec = mean(decr)
; Set up the tangent plane projection
step = 0.25/3600.0d0
delta_dec = range(decr)
delta_ra = range(rar)*cos(cendec/!radeg)
nx = ceil(delta_ra*1.05/step)
ny = ceil(delta_dec*1.05/step)
xref = nx/2
yref = ny/2
im = intarr(nx,ny)
MKHDR,head,im
SXADDPAR,head,'CDELT1',step
SXADDPAR,head,'CRPIX1',xref+1L
SXADDPAR,head,'CRVAL1',cenra
SXADDPAR,head,'CTYPE1','RA---TAN'
SXADDPAR,head,'CDELT2',step
SXADDPAR,head,'CRPIX2',yref+1L
SXADDPAR,head,'CRVAL2',cendec
SXADDPAR,head,'CTYPE2','DEC--TAN'

print,'RA range = [',strtrim(rar[0],2),',',strtrim(rar[1],2),'] deg'
print,'DEC range = [',strtrim(decr[0],2),',',strtrim(decr[1],2),'] deg'
print,'Central RA = ',strtrim(cenra,2)
print,'Central DEC = ',strtrim(cendec,2)
print,'NX = ',strtrim(nx,2)
print,'NY = ',strtrim(ny,2)

filters = ['u','g','r','i','z']
nfilters = n_elements(filters)

; Filter loop
FOR f=0,nfilters-1 do begin

  ifilter = filters[f]
  print,strtrim(f+1,2),' FILTER=',ifilter

  ; Get the FITS files for this filter
  ind = where(chstr.filter eq ifilter,nind)

  ; Loop through the files
  For i=0,nind-1 do begin
    chstr1 = chstr[ind[i]]
    CALDAT,chstr1.mjd+2400000.5d0-0.5,month,day,year,hour,minute,second
    night = string(year,format='(i04)')+string(month,format='(i02)')+string(day,format='(i02)')
    file = reduxdir+night+'/'+strtrim(chstr1.field,2)+'/'+strtrim(chstr1.file,2)
    print,' ',strtrim(f+1,2),' ',strtrim(i+1,2),'/',strtrim(nind,2),' ',file
    FITS_READ,file,im1,head1,message=error,/no_abort
    if error ne '' then stop,'PROBLEM loading '+file
    sz1 = size(im1)
    exptime = chstr1.exptime

    ; Perimeter pixels
    xperi = [lindgen(sz1[1]), lonarr(sz1[2]-2)+(sz1[1]-1), (sz1[1]-1)-lindgen(sz1[1]), lonarr(sz1[2]-2)]
    yperi = [lonarr(sz1[1]-1), lindgen(sz1[2]), lonarr(sz1[1]-2)+(sz1[2]-1), (sz1[2]-1)-lindgen(sz1[2]-1)]
    head_xyad,head1,xperi,yperi,rr_peri,dd_peri,/deg
    head_adxy,head,rr_peri,dd_peri,xx2_peri,yy2_peri,/deg

    ; Region of final image covered by this one
    xmn = floor(min(xx2_peri))
    xmx = ceil(max(xx2_peri))
    ymn = floor(min(yy2_peri))
    ymx = ceil(max(yy2_peri))
    nxout = xmx-xmn+1
    nyout = ymx-ymn+1
    xout = lindgen(nxout)+xmn
    yout = lindgen(nyout)+ymn
    expim1 = intarr(nxout,nyout)

    ; Get the region covered by our image
    cutind = polyfillv(xx2_peri-xmn,yy2_peri-ymn,nxout,nyout)
    expim1[cutind] = exptime

    ; Get the bad pixels
    gdpix1 = where(im1 lt 59000L,ngdpix1,comp=bdpix1,ncomp=nbdpix1)
    dum = array_indices(im1,bdpix1)
    xx_bd = reform(dum[0,*])
    yy_bd = reform(dum[1,*])
    ; Bad pixels, convert from X/Y to RA/DEC and then Xfinal/yfinal
    head_xyad,head1,xx_bd,yy_bd,rr_bd,dd_bd,/deg
    head_adxy,head,rr_bd,dd_bd,xx2_bd,yy2_bd,/deg
    xout_bd = xx2_bd - xmn  ; indices in our subregion image
    yout_bd = yy2_bd - ymn  ; indices in our subregion image
    badmask = intarr(nxout,nyout)
    ; spread it around a bit
    badmask[floor(xout_bd),floor(yout_bd)] = 1
    badmask[floor(xout_bd),ceil(yout_bd)] = 1
    badmask[ceil(xout_bd),floor(yout_bd)] = 1
    badmask[ceil(xout_bd),ceil(yout_bd)] = 1

    ; Set bad pixels to zero
    expim1 *= (1-badmask)

    ; Now add to final image
    im[xmn:xmx,ymn:ymx] += expim1

    print,ngdpix1,nbdpix1
  Endfor  ; chip loop

  ; Update the header
  sxaddpar,head,'filter',ifilter

  ; Write to file
  if f eq 0 then begin
    MWRFITS,im,outfile,head,/create
  endif else begin
    ; Fix the header for extension
    if strmid(head[0],0,6) eq 'SIMPLE' then head[0]="XTENSION= 'IMAGE   '           / IMAGE extension"
    MWRFITS,im,outfile,head
  endelse

ENDFOR

; Compress the file
;print,'Compressing output file'
;spawn,['gzip','-f',outfile],out,errout,/noshell

end
