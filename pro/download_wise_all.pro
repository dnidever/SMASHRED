pro download_wise_all

outdir = '/data/smash/cp/red/photred/wise/'

smashred_getredinfo,info,/silent
bd = where(strtrim(info.field,2) eq '',nbd)
if nbd gt 0 then remove,bd,info
ui = uniq(info.field,sort(info.field))
fields = info[ui].field
nfields = n_elements(fields)

nmulti = 3 ;5

cmd = 'download_wise,"'+fields+'",/compress'
dirs = '/data/smash/cp/red/photred/wise/tmp/'
PBS_DAEMON,cmd,dirs,/hyperthread,/idle,prefix='dwnwise',nmulti=nmulti,wait=5

stop

end
