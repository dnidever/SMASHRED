pro run_stargalaxy_separation_single,field,version,deep=deep,redo=redo

;; Make "star" catalogs.

if n_elements(field) eq 0 then begin
  print,'Syntax - run_stargalaxy_separation_single,field,version,deep=deep,redo=redo'
  return
endif

rootdir = smashred_rootdir()
if n_elements(version) eq 0 then version='v5'
catdir = rootdir+'cp/red/photred/catalogs/final/'+version+'/'
if keyword_set(deep) then tag='_deep' else tag=''


; Does the output file already exist
outfile = catdir+'stars/'+field+'_allobj'+tag+'_stars.fits'
if (file_test(outfile) eq 1 or file_test(outfile+'.gz') eq 1) and not keyword_set(redo) then begin
  print,outfile,' exists and /redo not set'
  return
endif

print,'Performing star/galaxy separation for field ',field
if keyword_set(deep) then print,'Using "deep-only" catalog'

;; Load the allobj file
file = catdir+field+'_combined_allobj'+tag+'.fits.gz'
if file_test(file) eq 0 then begin
  print,file+' NOT FOUND'
  return
endif
obj0 = mrdfits(file,1,/silent)
tags = tag_names(obj0)

; Maybe only look at depthflag=2 or 3
  
;; Deredden
obj = obj0
;print,'NO DEREDDENING'
;print,'  Dereddening'
;magext = [4.239, 3.303, 2.285, 1.698, 1.263]
;obj.u -= magext[0]*obj.ebv
;obj.g -= magext[1]*obj.ebv
;obj.r -= magext[2]*obj.ebv
;obj.i -= magext[3]*obj.ebv
;obj.z -= magext[4]*obj.ebv

;; Select the stars
gdstars = where(abs(obj.sharp) lt 1 and obj.chi lt 2 and obj.prob gt 0.2 and $
                obj.ndet gt 5 and obj.depthflag gt 1,ngdstars)
obj = obj[gdstars]

;; Morphology cuts
smashred_morphcuts,obj,ind1,nsig=3
obj1 = obj[ind1]
  
;; Color-color cuts
smashred_2cdcuts,obj1,ind2
obj2 = obj1[ind2]
; this cuts out some stars with g>24 on Field56

; Write out pruned catalogs
print,'Writing out pruned catalogs'
if file_test(catdir+'stars/',/directory) eq 0 then file_mkdir,catdir+'stars/'
file_delete,[outfile,outfile+'.gz'],/allow
MWRFITS,obj2,outfile,/create
spawn,['gzip',outfile],out,errout,/noshell   ; compress

;stop

end
