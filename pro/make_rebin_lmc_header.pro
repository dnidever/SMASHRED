pro make_rebin_lmc_header,tilehead

;; Make the final WCS and header
;; USE MAGELLANIC STREAM COORDINATES
;; -8 < MLON < +7
;; -10 < MLAT < +7.5
;; at 5" resolution that is 10800 x 11160 pixels
nx = 10800L
ny = 13300L
xref = 5900L
yref = 9100L
;cenmlon = -0.5
;cenmlat = -0.25
step = 5.0d0 / 3600.0d0  ; 5"
cenra = 81.90d0
cendec = -69.87d0
MKHDR,tilehead,fltarr(5,5)
SXADDPAR,tilehead,'NAXIS1',nx
SXADDPAR,tilehead,'CDELT1',step
SXADDPAR,tilehead,'CRPIX1',xref+1L
SXADDPAR,tilehead,'CRVAL1',cenra
SXADDPAR,tilehead,'CTYPE1','RA---TAN'
SXADDPAR,tilehead,'NAXIS2',ny
SXADDPAR,tilehead,'CDELT2',step
SXADDPAR,tilehead,'CRPIX2',yref+1L
SXADDPAR,tilehead,'CRVAL2',cendec
SXADDPAR,tilehead,'CTYPE2','DEC--TAN'
;x1 = scale_vector(findgen(100),0,nx-1)
;y1 = scale_vector(findgen(100),0,ny-1)
;xx = x1#replicate(1,100)
;yy = replicate(1,100)#y1
;head_xyad,tilehead,xx,yy,ra,dec,/deg
;glactc,ra,dec,2000.0,glon,glat,1,/deg
;gal2mag,glon,glat,mlon,mlat
;print,minmax(mlon)
;print,minmax(mlat)

;stop

end
