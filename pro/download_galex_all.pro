pro download_galex_all

;outdir = '/data/smash/cp/red/photred/gaia/'
outdir = '/data/smash/cp/red/photred/galex/'

smashred_getredinfo,info,/silent
bd = where(strtrim(info.field,2) eq '',nbd)
if nbd gt 0 then remove,bd,info
ui = uniq(info.field,sort(info.field))
fields = info[ui].field
nfields = n_elements(fields)

cmd = 'download_galex,"'+fields+'"'
dirs = '/data/smash/cp/red/photred/galex/tmp/'
PBS_DAEMON,cmd,dirs,/hyperthread,/idle,prefix='dwngaia',nmulti=20,wait=30

stop

end
