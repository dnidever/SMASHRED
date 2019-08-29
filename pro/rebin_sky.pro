pro rebin_sky,im,sky

nim = n_elements(im)
if nim eq 0 then begin
  print,'Syntax - rebin_sky,im,sky'
  return
endif

satlim = 60000L
sz = size(im)
nx = sz[1]
ny = sz[2]

  ;-- Compute background image --

  ; Computing sky level and sigma
  photred_sky,im,skymode,skysig1,highbad=satlim*0.95,/silent
  if skysig1 lt 0.0 then skysig1 = mad(im[gdpix])
  if skysig1 lt 0.0 then skysig1 = mad(im)

  ; First pass, no clipping (except for saturated pixels)
  backgim_large = im
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

  ;; Maybe fit a simple linear model to the background image

  sky = backgim2

;stop

end
