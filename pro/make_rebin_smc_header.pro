pro make_rebin_smc_header,tilehead

;; Make the final WCS and header
;; USE MAGELLANIC STREAM COORDINATES
;; -8 < MLON < +7
;; -10 < MLAT < +7.5
;; at 5" resolution that is 10800 x 11160 pixels

;; about 7 deg in dec and ~8 deg in ra
nx = 6200L
ny = 5300L
xref = 2500L
yref = 2900L
;cenmlon = -0.5
;cenmlat = -0.25
step = 5.0d0 / 3600.0d0  ; 5"
cenra = 13.183333d0
cendec = -72.828333d0
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

;head_adxy,tilehead,chstr.ra,chstr.dec,x,y,/deg
;plotc,x,y,ps=3,/xflip,xr=[-1000,8000],yr=[-1000,6000],xs=1,ys=1,charsize=1.5   
;oplot,[0,nx-1,nx-1,0,0],[0,0,ny-1,ny-1,0],co=250 

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
