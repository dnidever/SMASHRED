pro download_gaia_all

outdir = '/data/smash/cp/red/photred/gaia/'

smashred_getredinfo,info,/silent
bd = where(strtrim(info.field,2) eq '',nbd)
if nbd gt 0 then remove,bd,info
ui = uniq(info.field,sort(info.field))
fields = info[ui].field
nfields = n_elements(fields)

cmd = 'download_gaia,"'+fields+'"'
dirs = '/data/smash/cp/red/photred/gaia/tmp/'
PBS_DAEMON,cmd,dirs,/hyperthread,/idle,prefix='dwngaia',nmulti=20,wait=30

stop

end
