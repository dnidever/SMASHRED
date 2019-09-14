pro fits2db_conv_source_dr2_all

version = 'v6'
reduxdir = smashred_rootdir()+'cp/red/photred/'
tmpdir = reduxdir+'catalogs/final/'+version+'/tmp/'

lchstr = mrdfits(reduxdir+'smash_chips_lmc.fits.gz',1)
schstr = mrdfits(reduxdir+'smash_chips_smc.fits.gz',1)
mchstr = [lchstr,schstr]
mchstr.smash_field = strtrim(mchstr.smash_field,2)
ui = uniq(mchstr.smash_field,sort(mchstr.smash_field))
inner_fields = mchstr[ui].smash_field
;; 81 fields

chstr = mrdfits(reduxdir+'smash_chips.fits.gz',1)
chstr.smash_field = strtrim(chstr.smash_field,2)
ui = uniq(chstr.smash_field,sort(chstr.smash_field))
all_fields = chstr[ui].smash_field
periphery_fields = all_fields
match,periphery_fields,inner_fields,ind1,ind2,/sort
remove,ind1,periphery_fields
;; 116 periphery fields

;; Healpix 
hfiles = file_search(reduxdir+'catalogs/final/'+version+'/4*_combined_allsrc.fits.gz',count=nhfiles)
hpix = file_basename(hfiles,'_combined_allsrc.fits.gz')

;; Commands
fields = [periphery_fields,hpix]
allcmd = "fits2db_conv_source_dr2,'"+fields+"'"
nallcmd = n_elements(allcmd)
alldirs = strarr(nallcmd)+tmpdir

;; Randomize
rnd = sort(randomu(0,nallcmd))
fields = fields[rnd]
allcmd = allcmd[rnd]
alldirs = alldirs[rnd]

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


nmulti = 5

stop

pbs_daemon,cmd,dirs,jobs=jobs,/idle,/hyperthread,nmulti=nmulti,wait=10,prefix='fits2db'

stop

end
