pro make_cluster_cats

; Make cluster catalogs
dir = '/dl1/users/dnidever/smash/cp/red/photred/'
catdir = dir+'catalogs/final/v6/'
clusterdir = dir+'clusters/'

; Load the cluster catalog
cat = importascii(clusterdir+'smash_clusters.txt',/header)
ncat = n_elements(cat)
smash = importascii('~/projects/SMASHRED/data/smash_fields_final.txt',/header)


for i=0,ncat-1 do begin
  print,strtrim(i+1,2),'/',strtrim(ncat,2),' ',cat[i].name
  outfile = catdir+cat[i].name+'_cat.fits'
  if file_test(outfile) eq 0 or keyword_set(redo) then begin
    ; Get SMASH field name
    dist = sphdist(smash.radeg,smash.dedeg,cat[i].ra,cat[i].dec,/deg)
    ind = first_el(minloc(dist))
    field = 'Field'+strtrim(smash[ind[0]].num,2)
    ; Load the data
    objfile = catdir+field+'_combined_allobj.fits.gz'
    print,'  Loading ',objfile
    obj = mrdfits(objfile,1)
    ; Select the cluster stars
    dist = sphdist(cat[i].ra,cat[i].dec,obj.ra,obj.dec,/deg)*60.
    gd = where(dist lt 2.5*cat[i].size,ngd)  ;3.0
    cobj = obj[gd]
    ; Save the output catalog
    print,'  Writing to ',outfile
    MWRFITS,cobj,outfile,/create   
    stop
  endif else print,outfile,' EXISTS and /redo NOT set'
endfor

stop

end
