pro qacheck_spatialtransvar2_results,version

; Analysis thes results of running qacheck_spatialtransvar2.pro
; on all the fields and find bad exposures

if n_elements(version) eq 0 then version = 'v2'
dir = '/data/smash/cp/red/photred/catalogs/final/'+version+'/qaspatialtransvar/'

files = file_search(dir+'Field*fits',count=nfiles)
filters = ['u','g','r','i','z']
nfilters = n_elements(filters)
for i=0,nfiles-1 do begin
  str1=mrdfits(files[i],1)
  add_tag,str1,'field','',str1
  add_tag,str1,'filtnum',0,str1
  add_tag,str1,'relmed',0.0,str1
  add_tag,str1,'bad',0B,str1
  field=first_el(strsplit(file_basename(files[i],'.fits'),'_',/extract))
  str1.field=field
  push,all,str1
endfor
nall = n_elements(all)

g=where(all.filter eq 'u') & all[g].filtnum=1
g=where(all.filter eq 'g') & all[g].filtnum=2
g=where(all.filter eq 'r') & all[g].filtnum=3
g=where(all.filter eq 'i') & all[g].filtnum=4
g=where(all.filter eq 'z') & all[g].filtnum=5

; Get relative photometric zeropoints
for i=0,nfilters-1 do begin
  ind = where(all.filter eq filters[i],nind)
  med = median(all[ind].med)
  all[ind].relmed = all[ind].med-med
endfor

bd = where(all.relmed gt 0.5 and all.diff_sig gt 0.04,nbd)
bad = strarr(nall)
bad[bd] = '  BAD!!!'
printline,all.field+'  '+all.expnum+'  '+all.filter+'  '+stringize(all.exptime,ndec=1)+'  '+$
          stringize(all.relmed,ndec=3)+'  '+stringize(all.diff_sig,ndec=3)+bad

stop

end
