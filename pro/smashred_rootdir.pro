;+
;
; SMASHRED_ROOTDIR
;
; Returns the root directory for the SMASH data
; depending on the host.
;
; By D.Nidever  Jan 2071
;-
function smashred_rootdir

spawn,'hostname',out,errout,/noshell
host = (strsplit(out[0],' ',/extract))[0]

if stregex(host,'bambam',/boolean) eq 1 then return,'/data/smash/'
if stregex(host,'dldb1',/boolean) eq 1 or stregex(host,'zeus1',/boolean) eq 1 or stregex(host,'datalab.noao.edu',/boolean) eq 1 then return,'/datalab/users/dnidever/smash/'
if stregex(host,'hulk',/boolean) eq 1 or stregex(host,'thing',/boolean) eq 1 then return,'/dl1/users/dnidever/smash/'

; Not sure about this host
print,host,' UNKNOWN'
return,'-1'

end
