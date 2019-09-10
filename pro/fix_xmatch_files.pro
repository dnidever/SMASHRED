pro fix_xmatch_files

;; Need to change GAIA_SOURCE from LONG to LONG64

dir = '/dl1/users/dnidever/smash/cp/red/photred/catalogs/final/v6/'
files = file_search(dir+'*_combined_allobj_xmatch2.fits.gz',count=nfiles)

for i=0,nfiles-1 do begin
  print,strtrim(i+1,2),' ',files[i]
  hd = headfits(files[i],exten=1)
  if strtrim(sxpar(hd,'TFORM4'),2) ne 'K' then begin
    str = mrdfits(files[i],1)
    nstr = n_elements(str)
    tags = tag_names(str)
    ntags = n_tags(str)
    typ = lonarr(ntags)
    for j=0,ntags-1 do typ[j]=size(str[0].(j),/type)
    gd = where(tags eq 'GAIA_SOURCE',ngd)
    typ[gd] = 14
    schema = {id:''}
    for j=1,ntags-1 do begin
      blank = fix('',type=typ[j])
      if typ[j] eq 1 then blank=0B
      schema = create_struct(schema,tags[j],blank)
    endfor
    new = replicate(schema,nstr)
    struct_assign,str,new
    outfile = repstr(files[i],'.fits.gz','.fits')
    MWRFITS,new,outfile,/create
    spawn,['gzip','-f',outfile],/noshell
  endif else print,files[i],' ALREADY CORRECT'
  ;stop
endfor

stop

end
