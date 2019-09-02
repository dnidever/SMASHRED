pro rebin_smc,filter,redo=redo

tmpdir = '/dl1/users/dnidever/smash/cp/red/photred/rebin/tmp/'
basedir = '/dl1/users/dnidever/smash/cp/red/photred/'

if n_elements(filter) eq 0 then begin
  print,'Syntax - rebin_smc,filter,redo=redo'
  return
endif

chstr = mrdfits('/dl1/users/dnidever/smash/cp/red/photred/catalogs/final/v6/smash_chips.fits.gz',1)
chstr.smash_field = strtrim(chstr.smash_field,2)

nmulti = 20

;; MAYBE REMOVE BAD EXPOSURES USING THE BAD EXPOSURE LIST
badexp = importascii('/home/dnidever/projects/SMASHRED/obslog/smash_badexposures.txt',/header)

;; Pick chips inside this limit
;filter = 'g'
glactc,chstr.ra,chstr.dec,2000.0,glon,glat,1,/deg
gal2mag,glon,glat,mlon,mlat
g = where(mlon ge -20.5 and mlon le -10.5 and mlat ge -16.6 and mlat le -6.0 and chstr.filter eq filter,ng)
;          chstr.smash_field ne 'Field19' and chstr.smash_field ne 'Field26' and chstr.smash_field ne 'Field25' and $
;          chstr.smash_field ne 'Field65' and chstr.smash_field ne 'Field60',ng)
;; some are getting cut off at the bottom!
chstr1 = chstr[g]
chstr1.night = strtrim(chstr1.night,2)
chstr1.field = strtrim(chstr1.field,2)
chstr1.file = strtrim(chstr1.file,2)

;; Get the FITS file names
files = basedir+chstr1.night+'/'+chstr1.field+'/'+chstr1.file+'.fz'
;; The long exposures point to the deep/FIELD/ directories
gdeep = where(chstr1.exptime gt 200,ngdeep)
files[gdeep] = basedir+'deep/'+chstr1[gdeep].smash_field+'/'+chstr1[gdeep].field+'/'+chstr1[gdeep].file+'.fz'

;; Put them in groups
npergroup = 20
ngroups = ceil(float(ng)/npergroup)
cmd = strarr(ngroups)
for i=0,ngroups-1 do cmd[i] = 'rebin_image,["'+strjoin(files[i*npergroup:((i+1)*npergroup-1)<(ng-1)],'","')+'"],lmc=0'
if keyword_set(redo) then cmd+=',/redo'
dirs = strarr(ngroups)+tmpdir

;stop

;PBS_DAEMON,cmd,dirs,jobs=jobs,prefix='rebin',/hyperthread,/idle,wait=1,nmulti=nmulti

;stop

MAKE_REBIN_SMC_HEADER,tilehead
nx = sxpar(tilehead,'NAXIS1')
ny = sxpar(tilehead,'NAXIS2')
step = sxpar(tilehead,'CDELT1')

head_adxy,tilehead,chstr1.ra,chstr1.dec,xch,ych,/deg

stop

fim = fltarr(nx,ny)
num = bytarr(nx,ny)

;; Now stack them all up
outfiles = basedir+'rebin/'+filter+'/'+chstr1.file+'.gz'
for i=0,ng-1 do begin
  if i mod 100 eq 0 then print,i
  ;if xch[i] lt 6000 or xch[i] gt 7000 or ych[i] lt 9000 or ych[i] gt 10000 then goto,BOMB
  ;if xch[i] lt 6000 or ych[i] gt 4000 then goto,BOMB
  if file_test(outfiles[i]) eq 1 then begin
    fits_read,outfiles[i],im1,head1
    ;; I messed up the names the first time
    ;xlo = sxpar(head1,'xlo')
    ;ylo = sxpar(head1,'xhi')
    ;xhi = sxpar(head1,'ylo')
    ;yhi = sxpar(head1,'yhi')
    xlo = sxpar(head1,'xlo')
    xhi = sxpar(head1,'xhi')
    ylo = sxpar(head1,'ylo')
    yhi = sxpar(head1,'yhi')
    nx1 = sxpar(head1,'naxis1')
    ny1 = sxpar(head1,'naxis2')
    if xlo lt 0 or xhi gt (nx-1) or ylo lt 0 or yhi gt (ny-1) then begin
      print,outfiles[i],' is out of bounds ',xlo,xhi,ylo,yhi
      goto,BOMB
    endif

    mask1 = (im1 ne 0.0)
    ;; Correct image for exposure time
    ;;  correct o 60s
    im1 *= 60/chstr1[i].exptime
    ;; Correct image for airmass
    ;;   calmag = instmag - t.zpterm - t.amterm*am
    ;;   positive magcorr makes it fainter
    ;;   fluxcorr will be greater than 1 so divide by it
    magcorr = -chstr1[i].zpterm - chstr1[i].amterm*chstr1[i].airmass
    fluxcorr = 10^(0.4*magcorr)
    im1 /= fluxcorr

    ;print,strtrim(i+1,2),' ',file_basename(outfiles[i]),' ',chstr1[i].exptime,' ',fluxcorr
    ;if xch[i] gt 6400 and xch[i] lt 6800 and ych[i] gt 9460 and ych[i] lt 9760 then stop,'check 30dor image'

    ;; Add it to the big arrays
    fim[xlo:xhi,ylo:yhi] += im1
    num[xlo:xhi,ylo:yhi] += mask1

    ;stop
  endif else print,outfiles[i],' NOT FOUND'
  ;stop
  BOMB:
endfor

;; Take the average
fim /= (num>1)

setdisp
displayc,fim,/log,min=1,/xflip

;stop

outfile = basedir+'rebin/smash_smc_'+filter+'.fits'
print,'Writing final image to ',outfile
MWRFITS,fim,outfile,tilehead,/create
spawn,['gzip','-f',outfile],/noshell

stop

end
