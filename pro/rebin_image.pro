pro rebin_image,files,redo=redo,lmc=lmc

;file = '/dl1/users/dnidever/smash/cp/red/photred/20160101/F1/F1-00507800_01.fits.fz'
;file = '/dl1/users/dnidever/smash/cp/red/photred/20160218/F8/F8-00518838_01.fits.fz'

nfiles = n_elements(files)
if nfiles eq 0 then begin
  print,'Syntax - rebin_image,files,redo=redo'
  return
endif

outdir = '/dl1/users/dnidever/smash/cp/red/photred/rebin/'

if keyword_set(lmc) then begin
  MAKE_REBIN_LMC_HEADER,tilehead
endif else begin
  MAKE_REBIN_SMC_HEADER,tilehead
endelse
nx = sxpar(tilehead,'NAXIS1')
ny = sxpar(tilehead,'NAXIS2')
step = sxpar(tilehead,'CDELT1')

bin = round(step*3600/0.262)
print,'Bin size = ',strtrim(bin,2),' pixels'
satlim = 60000L

;; Load the list of image to mask
maskstr = importascii('/home/dnidever/projects/SMASHRED/data/image_mask.txt',/header,/silent)
add_tag,maskstr,'base','',maskstr
maskstr.base = file_basename(maskstr.filename,'.fits.fz')

;; File loop
for i=0,nfiles-1 do begin

  file = files[i]
  print,strtrim(i+1,1),'/',strtrim(nfiles,2),' ',file

  if file_test(file) eq 0 then begin
    print,file,' NOT FOUND'
    goto,BOMB
  endif
  head = headfits(file,exten=1)
  dir = file_dirname(file)+'/'
  base = file_basename(file,'.fits.fz')
  chip = long(first_el(strsplit(base,'_',/extract),/last))
  filter = sxpar(head,'filter')
  filt = strmid(filter,0,1)
  outfile = outdir+filt+'/'+base+'.fits'

  if file_test(outfile+'.gz') eq 1 and not keyword_set(redo) then begin
    print,file,' already EXISTS and /redp NOT set'
    goto,BOMB
  endif

  ;; Load the data
  fits_read,file,im,head
  sz = size(im)
  nx = sz[1]
  ny = sz[2]

  ;; Mask bad parts
  MATCH,maskstr.base,base,ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then begin
    maskstr1 = maskstr[ind1]
    xx = findgen(nx)#replicate(1,ny)
    yy = replicate(1,nx)#findgen(ny)
    bd = where(abs(yy-(xx*maskstr1.slope+maskstr1.offset)) le maskstr1.width*0.5,nbd)
    if nbd gt 0 then im[bd]=60000.
    print,'Masking ',strtrim(nbd,2),' bad pixels'
  endif

  ;; Good pixel mask
  gmask = (im lt 59000L)

  dateobs = sxpar(head,'DATE-OBS')
  mjd = PHOTRED_GETMJD('','ctio',dateobs=dateobs)

  ;; Use Gaia head if possible
  gheadfile = dir+base+'.gaiawcs.head'
  if file_test(gheadfile) eq 1 then begin
    head0 = head
    READLINE,gheadfile,head
  endif

  ;; Load the subtracted image if it exists
  sfile = dir+base+'s.fits.fz'
  if file_test(sfile) eq 1 then begin
    fits_read,sfile,sim,shead
    backgim_large = sim
  endif else backgim_large = im

  ;; Try to fix background for Chip 31
  ;if (chip eq 31) and (mjd gt 56660) then begin
  if (chip eq 31) and (mjd gt 56400) then begin
    xoff = 1024
    if nx eq 2032 then xoff=1015
    med1 = median(im[0:xoff,*])
    med2 = median(im[xoff+1:*,*])
    rebin_sky,im[xoff+1:*,*],sky2
    im[0:xoff,*] -= med1
    im[xoff+1:*,*] -= sky2
  endif

  ;; Subtract the background

  ;; Computing sky level and sigma
  gdpix = where(im lt satlim*0.90 and im ne 0.0,ngdpix)
  photred_sky,im,skymode,skysig1,highbad=satlim*0.95,/silent
  if skysig1 lt 0.0 then skysig1 = mad(im[gdpix])
  if skysig1 lt 0.0 then skysig1 = mad(im)

  ;-- Compute background image --
  ;rebin_sky,im,backgim
  ;subim = im-backgim
  ;subim = im-median(backgim)
  subim = im
  subim[gdpix] = im[gdpix]-skymode

  ; Mask bad part of Chip 31
  ;if (chip eq 31) and (mjd gt 56660) then begin
  ;  ; X: 1-1024 okay
  ;  ; X: 1025-2049 bad
  ;  print,'Masking bad half of chip 1'
  ;  im[1025:*,*] = 0
  ;  gmask[1025:*,*] = 0 
  ;endif

  ;; Rebin
  nx2 = nx/bin
  ny2 = ny/bin
  im2 = REBIN(subim[0:nx2*bin-1,0:ny2*bin-1]*gmask[0:nx2*bin-1,0:ny2*bin-1],nx2,ny2)
  gpix2 = REBIN(float(gmask[0:nx2*bin-1,0:ny2*bin-1]),nx2,ny2)
  b = where(gpix2 eq 0,nb)
  if nb gt 0 then gpix2[b]=1
  ;; Correct rebined/average image for masked pixels
  ;; im2 = Sum(image)/Nbin^2
  ;; gpix2 = Sum(good pixels)/Nbin^2
  ;; we want Sum(image)/Sum(good pixels)
  im2 /= gpix2

  ;; Interpolate on final grid
  x1 = findgen(nx2)*bin+bin/2
  y1 = findgen(ny2)*bin+bin/2
  xx = x1 # replicate(1,ny2)
  yy = replicate(1,nx2) # y1
  HEAD_XYAD,head,xx,yy,ra,dec,/deg
  HEAD_ADXY,tilehead,ra,dec,newx,newy,/deg
  TRIANGULATE, newx, newy, tr, b
  limits = [floor(min(newx)), floor(min(newy)), ceil(max(newx)), ceil(max(newy))]
  steps = [1.0,1.0]
  fim = TRIGRID(newx,newy,im2, tr, steps,limits,missing=0.0)
  fim = float(fim)
  bd = where(finite(fim) eq 0,nbd)
  if nbd gt 0 then fim[bd]=0


  MKHDR,newhead,fim
  sxaddpar,newhead,'XLO',limits[0]
  sxaddpar,newhead,'YLO',limits[1]
  sxaddpar,newhead,'XHI',limits[2]
  sxaddpar,newhead,'YHI',limits[3]

  ;; Save
  print,'Writing to ',outfile
  MWRFITS,fim,outfile,newhead,/create
  spawn,['gzip','-f',outfile],/noshell

  BOMB:
endfor

;stop

end
