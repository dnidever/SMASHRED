pro smashred_calibrate_healpix_all,version=version

;; Run smashred_calibrate_healpix on all the main body HEALPix

if n_elements(version) eq 0 then version='v6'
reduxdir = smashred_rootdir()+'cp/red/photred/'
tmpdir = reduxdir+'catalogs/final/'+version+'/tmp/'

nside = 64
radeg = 180.0d0 / !dpi
nmulti = 10

;; Load the chip structures for the LMC and SMC main body fields
lchstr = mrdfits(reduxdir+'smash_chips_lmc.fits.gz',1)
schstr = mrdfits(reduxdir+'smash_chips_smc.fits.gz',1)

;; Get HEALPix pixel numbers for the LMC and SMC chips
off = 0.5  ; 1.0
;; Wiggle the positions a bit to deal with the edges
ANG2PIX_RING,nside,(90-lchstr.dec)/radeg,lchstr.ra/radeg,lpix1
ANG2PIX_RING,nside,(90-lchstr.dec-off)/radeg,lchstr.ra/radeg,lpix2
ANG2PIX_RING,nside,(90-lchstr.dec+off)/radeg,lchstr.ra/radeg,lpix3
ANG2PIX_RING,nside,(90-lchstr.dec)/radeg,(lchstr.ra-off/cos(lchstr.dec/!radeg))/radeg,lpix4
ANG2PIX_RING,nside,(90-lchstr.dec)/radeg,(lchstr.ra+off/cos(lchstr.dec/!radeg))/radeg,lpix5
lpix = [lpix1,lpix2,lpix3,lpix4,lpix5]
lpix = lpix[uniq(lpix,sort(lpix))]
PIX2ANG_RING,nside,lpix,ltheta,lphi
lpra = lphi*radeg
lpdec = 90-ltheta*radeg
glactc,lpra,lpdec,2000.0,lpgl,lpgb,1,/deg
gal2mag,lpgl,lpgb,lpml,lpmb

ANG2PIX_RING,nside,(90-schstr.dec)/radeg,schstr.ra/radeg,spix1
ANG2PIX_RING,nside,(90-schstr.dec-off)/radeg,schstr.ra/radeg,spix2
ANG2PIX_RING,nside,(90-schstr.dec+off)/radeg,schstr.ra/radeg,spix3
ANG2PIX_RING,nside,(90-schstr.dec)/radeg,(schstr.ra-off/cos(schstr.dec/!radeg))/radeg,spix4
ANG2PIX_RING,nside,(90-schstr.dec)/radeg,(schstr.ra+off/cos(schstr.dec/!radeg))/radeg,spix5
spix = [spix1,spix2,spix3,spix4,spix5]
spix = spix[uniq(spix,sort(spix))]
PIX2ANG_RING,nside,spix,stheta,sphi
spra = sphi*radeg
spdec = 90-stheta*radeg
glactc,spra,spdec,2000.0,spgl,spgb,1,/deg
gal2mag,spgl,spgb,spml,spmb

pix = [lpix,spix]
npix = n_elements(pix)

;plot,mlon,mlat,ps=3,xr=[10,-22],yr=[-25,10],xs=1,ys=1
;oplot,ml[glmc],mb[glmc],ps=3,co=200
;oplot,ml[gsmc],mb[gsmc],ps=3,co=150
;oplot,pmlon,pmlat,ps=1,co=80
;oplot,lpml,lpmb,ps=1,co=250,sym=2
;oplot,spml,spmb,ps=1,co=200,sym=2

print,strtrim(npix,2),' HEALpix to run'

;; Create the commands
allcmd = 'smashred_calibrate_healpix,'+strtrim(pix,2)+',version="v6"'
alldirs = strarr(npix)+tmpdir
nallcmd = n_elements(allcmd)

;; Randomize
rnd = sort(randomu(0,nallcmd))
pix = pix[rnd]
allcmd = allcmd[rnd]
alldirs = alldirs[rnd]

;; Check for previous successes
outfiles = reduxdir+'catalogs/final/'+version+'/'+strtrim(pix,2)+'_combined_allobj.fits.gz'
test = file_test(outfiles)
bd = where(test eq 1,nbd,comp=gd,ncomp=ngd)
if nbd gt 0 then begin
  print,'Removing ',strtrim(nbd,2),' previous successes'
  allpix = pix
  pix = allpix[gd]
  allcmd = allcmd[gd]
  alldirs = alldirs[gd]
  nallcmd = n_elements(allcmd)
endif

;; Parcel out the jobs
spawn,'hostname',out,errout,/noshell
hostname = (strsplit(out[0],' ',/extract))[0]
thishost = first_el(strsplit(hostname,'.',/extract))
hosts = ['hulk','thing','gp09']
nhosts = n_elements(hosts)
torun = lindgen(nallcmd)
nperhost = nallcmd/nhosts
for i=0,nhosts-1 do $
  if stregex(thishost,hosts[i],/boolean) eq 1 then torun=torun[i*nperhost:(i+1)*nperhost-1]
ntorun = n_elements(torun)
cmd = allcmd[torun]
dirs = alldirs[torun]
print,'Running ',strtrim(n_elements(torun),2),' on ',thishost


stop

PBS_DAEMON,cmd,dirs,jobs=jobs,prefix='smashcalhpix',/hyperthread,/idle,nmulti=nmulti,wait=10


stop

end
