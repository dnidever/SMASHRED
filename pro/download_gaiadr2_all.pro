pro download_gaiadr2_all

outdir = '/data/smash/cp/red/photred/gaiadr2/'

smashred_getredinfo,info,/silent
bd = where(strtrim(info.field,2) eq '',nbd)
if nbd gt 0 then remove,bd,info
ui = uniq(info.field,sort(info.field))
fields = info[ui].field
nfields = n_elements(fields)

;;; Do LAF fields field
;smash = importascii('~/projects/SMASHRED/data/smash_fields_final.txt',/header)
;g = where(smash.mslon gt 30,ng)
;fields = 'Field'+strtrim(smash[g].num,2)
;nfields = ng

cmd = 'download_gaiadr2,"'+fields+'"'
dirs = '/dl1/users/dnidever/smash/cp/red/photred/gaiadr2/tmp/'
PBS_DAEMON,cmd,dirs,/hyperthread,/idle,prefix='dwngaia',nmulti=10,wait=30

stop

end
