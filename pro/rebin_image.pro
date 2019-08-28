pro rebin_image,file

file = '/dl1/users/dnidever/smash/cp/red/photred/20160101/F1/F1-00507800_01.fits.fz'

;; Make the final WCS and header
;; USE MAGELLANIC STREAM COORDINATES
;; -8 < MLON < +7
;; -8 < MLAT < +7.5
;; at 5" resolution that is 10800 x 11160 pixels
nx = 10800L
ny = 11300L
xref = 5900L
yref = 7500L
;cenmlon = -0.5
;cenmlat = -0.25
step = 5.0d0 / 3600.0d0  ; 5"
cenra = 81.90d0
cendec = -69.87d0
MKHDR,head,fltarr(5,5)
SXADDPAR,head,'NAXIS1',nx
SXADDPAR,head,'CDELT1',step
SXADDPAR,head,'CRPIX1',xref+1L
SXADDPAR,head,'CRVAL1',cenra
SXADDPAR,head,'CTYPE1','RA---TAN'
SXADDPAR,head,'NAXIS2',ny
SXADDPAR,head,'CDELT2',step
SXADDPAR,head,'CRPIX2',yref+1L
SXADDPAR,head,'CRVAL2',cendec
SXADDPAR,head,'CTYPE2','DEC--TAN'
;x1 = scale_vector(findgen(100),0,nx-1)
;y1 = scale_vector(findgen(100),0,ny-1)
;xx = x1#replicate(1,100)
;yy = replicate(1,100)#y1
;head_xyad,head,xx,yy,ra,dec,/deg
;glactc,ra,dec,2000.0,glon,glat,1,/deg
;gal2mag,glon,glat,mlon,mlat
;print,minmax(mlon)
;print,minmax(mlat)

dir = file_dirname(file)+'/'
base = file_basename(file,'.fits.fz')

bin = round(step*3600/0.262)
print,'Bin size = ',strtrim(bin,2),' pixels'
satlim = 60000L

;; Load the data
fits_read,file,im,head
sz = size(im)
nx = sz[1]
ny = sz[2]
mask = (im gt 59000L)

;; Load the subtracted image if it exists
sfile = dir+base+'s.fits.fz'
if file_test(sfile) eq 1 then begin
  fits_read,sfile,sim,shead
  backgim_large = sim
endif else backgim_large = im

;; Subtract the background

;-- Compute background image --

; Computing sky level and sigma
photred_sky,im,skymode,skysig1,highbad=satlim*0.95,/silent
if skysig1 lt 0.0 then skysig1 = mad(im[gdpix])
if skysig1 lt 0.0 then skysig1 = mad(im)

; First pass, no clipping (except for saturated pixels)
;backgim_large = im
; Set saturated pixels to NaN so they won't be used in the smoothing
bdpix = where(backgim_large gt satlim*0.95,nbdpix)
if nbdpix gt 0 then backgim_large[bdpix] = !values.f_nan
sm = (400 < (nx/2.0) ) < (ny/2.0)
backgim_large = smooth(backgim_large,[sm,sm],/edge_truncate,/nan,missing=skymode)

; Second pass, use clipping, and first estimate of background
backgim1 = im
; Setting hi/low pixels to NaN, they won't be used for the smoothing
;bd = where(abs(im-skymode) gt 2.0*skysig1,nbd)
bd1 = where(abs(backgim1-backgim_large) gt 3.0*skysig1,nbd1)
if nbd1 gt 0 then (backgim1)(bd1) = !values.f_nan 
sm = (400 < (nx/2.0) ) < (ny/2.0)
backgim1 = smooth(backgim1,[sm,sm],/edge_truncate,/nan,missing=skymode)

; Third pass, use better estimate of background
backgim2 = im
; Setting hi/low pixels to NaN, they won't be used for the smoothing
;bd = where(abs(im-skymode) gt 2.0*skysig1,nbd)
bd2 = where(abs(backgim2-backgim1) gt 3.0*skysig1,nbd2)
if nbd2 gt 0 then (backgim2)(bd2) = !values.f_nan 
sm = (400 < (nx/2.0) ) < (ny/2.0)
backgim2 = smooth(backgim2,[sm,sm],/edge_truncate,/nan,missing=skymode)

subim = im-backgim2

;; Rebin
nx2 = nx/bin
ny2 = ny/bin
im2 = REBIN(subim[0:nx2*bin-1,0:ny2*bin-1],nx2,ny2)

;; Mask

stop

;; Interpolate on final grid
xx = findgen(sz1[1])#replicate(1,sz1[2])
yy = replicate(1,sz1[1])#findgen(sz1[2])
xyad,head1b,xx,yy,ra,dec
glactc,ra,dec,2000.0,glon1b,glat1b,1,/deg
TRIANGULATE, glon1b, glat1b, tr, b
limits = [min(gl), min(gb), max(gl), max(gb)]
steps = [sxpar(head,'CDELT1'), sxpar(head,'CDELT2')]
im = TRIGRID(glon1b,glat1b,reform(cube1b[*,*,i]), tr, steps,limits,missing=0.0)
bd = where(finite(im) eq 0,nbd)
if nbd gt 0 then im[bd]=0

stop

end
