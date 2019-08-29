pro rebin_lmc

tmpdir = '/dl1/users/dnidever/smash/cp/red/photred/rebin/tmp/'
basedir = '/dl1/users/dnidever/smash/cp/red/photred/'

chstr = mrdfits('/dl1/users/dnidever/smash/cp/red/photred/catalogs/final/v6/smash_chips.fits.gz',1)
chstr.smash_field = strtrim(chstr.smash_field,2)

nmulti = 20

;; Pick chips inside this limit
glactc,chstr.ra,chstr.dec,2000.0,glon,glat,1,/deg
gal2mag,glon,glat,mlon,mlat
g = where(mlon ge -8.1 and mlon le 7.2 and mlat ge -10.0 and mlat le 7.6 and chstr.filter eq 'g' and $
          chstr.smash_field ne 'Field19' and chstr.smash_field ne 'Field26' and chstr.smash_field ne 'Field25' and $
          chstr.smash_field ne 'Field65' and chstr.smash_field ne 'Field60',ng)
;; some are getting cut off at the bottom!
chstr1 = chstr[g]
chstr1.night = strtrim(chstr1.night,2)
chstr1.field = strtrim(chstr1.field,2)
chstr1.file = strtrim(chstr1.file,2)

files = basedir+chstr1.night+'/'+chstr1.field+'/'+chstr1.file+'.fz'

cmd = 'rebin_image,"'+files+'"'
dirs = strarr(n_elements(cmd))+tmpdir

PBS_DAEMON,cmd,dirs,jobs=jobs,prefix='rebin',/hyperthread,/idle,wait=1,nmulti=nmulti

stop

end
