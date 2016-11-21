pro smashred_matchcats,field,version=version,redo=redo

; Crossmatch SMASH catalog with GAIA, 2MASS and ALLWISE

; Not enough inputs
if n_elements(field) eq 0 then begin
  print,'Syntax - smashred_matchcats,field,version=version,redo=redo'
  return
endif

if n_elements(version) eq 0 then version='v3'
dir = '/data/smash/cp/red/photred/catalogs/final/'+version+'/'

print,'Running SMASHRED_MATCHCATS on ',field,' version=',version

if file_test(dir+field+'_combined_allobj.fits.gz') eq 0 then begin
  print,'ALLOBJ file not found for ',field,' and version ',version
  return
endif

outfile = dir+field+'combined_allobj_xmatch.fits'
if (file_test(outfile) eq 1 or file_test(outfile+'.gz') eq 1) and not keyword_set(redo) then begin
  print,'Output file exists and /redo NOT set'
  return
endif

; Restore the SMASH catalog
oldobj = MRDFITS(dir+field+'_combined_allobj.fits.gz',1)
nallobj = n_elements(oldobj)

; Copy to new schema
nan = !values.f_nan
allobj_schema = {id:'',ra:0.0d0,dec:0.0d0,raerr:99.99,decerr:99.99,rascatter:99.99,decscatter:99.99,ndet:0,depthflag:0B,$
                 u:99.99,uerr:9.99,uscatter:99.99,ndetu:0,g:99.99,gerr:9.99,$
                 gscatter:99.99,ndetg:0,r:99.99,rerr:9.99,rscatter:99.9,ndetr:0,i:99.99,ierr:9.99,iscatter:99.99,ndeti:0,$
                 z:99.99,zerr:9.99,zscatter:99.99,ndetz:0,chi:nan,sharp:nan,flag:-1,prob:nan,ebv:99.99,$
                 gaia_match:0B,gaia_matchdist:99.99,gaia_source:-1L,gaia_gmag:99.99,gaia_gmagerr:99.99,tmass_match:0B,tmass_matchdist:99.99,$
                 tmass_id:'',tmass_jmag:99.99,tmass_jmagerr:99.99,tmass_hmag:99.99,tmass_hmagerr:99.99,tmass_kmag:99.99,tmass_kmagerr:99.99,tmass_qflg:'',$
                 wise_match:0B,wise_matchdist:99.99,wise_id:'',wise_w1mag:99.99,wise_w1magerr:99.99,wise_w2mag:99.99,wise_w2magerr:99.99,$
                 wise_w3mag:99.99,wise_w3magerr:99.99,wise_w4mag:99.99,wise_w4magerr:99.99,wise_qph:''}
allobj = REPLICATE(allobj_schema,nallobj)
STRUCT_ASSIGN,oldobj,allobj,/nozero
undefine,oldobj

; Restore the GAIA file
print,'Matching with GAIA'
gaiadir = '/data/smash/cp/red/photred/gaia/'
if file_test(gaiadir+field+'_gaia.fits.gz') eq 0 then begin
  print,'No GAIA file found for field ',field
  return
endif
gaia = MRDFITS(gaiadir+field+'_gaia.fits.gz',1)
SRCMATCH,allobj.ra,allobj.dec,gaia.ra_icrs,gaia.de_icrs,0.5,ind1,ind2,/sph,count=nmatch
print,strtrim(nmatch,2),' GAIA matches'
if nmatch gt 0 then begin
  allobj[ind1].gaia_source = gaia[ind2].source
  allobj[ind1].gaia_gmag = gaia[ind2]._gmag_
  gmagerr = 2.5*alog10(1.0+gaia[ind2].e__fg_/gaia[ind2]._fg_)
  allobj[ind1].gaia_gmagerr = gmagerr
  dist = sphdist(allobj[ind1].ra,allobj[ind1].dec,gaia[ind2].ra_icrs,gaia[ind2].de_icrs,/deg)*3600.0d0
  allobj[ind1].gaia_matchdist = dist
endif

; Restore the 2MASS file
print,'Matching with 2MASS'
tmassdir = '/data/smash/cp/red/photred/tmass/'
if file_test(tmassdir+field+'_tmass.fits.gz') eq 0 then begin
  print,'No 2MASS file found for field ',field
  return
endif
tmass = MRDFITS(tmassdir+field+'_tmass.fits.gz',1)
SRCMATCH,allobj.ra,allobj.dec,tmass.raj2000,tmass.dej2000,0.5,ind1,ind2,/sph,count=nmatch
print,strtrim(nmatch,2),' 2MASS matches'
if nmatch gt 0 then begin
  allobj[ind1].tmass_id = tmass[ind2]._2mass
  allobj[ind1].tmass_jmag = tmass[ind2].jmag
  allobj[ind1].tmass_jmagerr = tmass[ind2].e_jmag
  allobj[ind1].tmass_hmag = tmass[ind2].hmag
  allobj[ind1].tmass_hmagerr = tmass[ind2].e_hmag
  allobj[ind1].tmass_kmag = tmass[ind2].kmag
  allobj[ind1].tmass_kmagerr = tmass[ind2].e_kmag
  allobj[ind1].tmass_qflg = tmass[ind2].qflg
  dist = sphdist(allobj[ind1].ra,allobj[ind1].dec,tmass[ind2].raj2000,tmass[ind2].dej2000,/deg)*3600.0d0
  allobj[ind1].tmass_matchdist = dist
endif

; Restore the WISE file
print,'Matching with WISE'
wisedir = '/data/smash/cp/red/photred/wise/'
if file_test(wisedir+field+'_wise.fits.gz') eq 0 then begin
  print,'No WISE file found for field ',field
  return
endif
wise = MRDFITS(wisedir+field+'_wise.fits.gz',1)
SRCMATCH,allobj.ra,allobj.dec,wise.raj2000,wise.dej2000,0.5,ind1,ind2,/sph,count=nmatch
print,strtrim(nmatch,2),' WISE matches'
if nmatch gt 0 then begin
  allobj[ind1].wise_id = wise[ind2].allwise
  allobj[ind1].wise_w1mag = wise[ind2].w1mag
  allobj[ind1].wise_w1magerr = wise[ind2].e_w1mag
  allobj[ind1].wise_w2mag = wise[ind2].w2mag
  allobj[ind1].wise_w2magerr = wise[ind2].e_w2mag
  allobj[ind1].wise_w3mag = wise[ind2].w3mag
  allobj[ind1].wise_w3magerr = wise[ind2].e_w3mag
  allobj[ind1].wise_w4mag = wise[ind2].w4mag
  allobj[ind1].wise_w4magerr = wise[ind2].e_w4mag
  allobj[ind1].wise_qph = wise[ind2].qph
  dist = sphdist(allobj[ind1].ra,allobj[ind1].dec,wise[ind2].raj2000,wise[ind2].dej2000,/deg)*3600.0d0
  allobj[ind1].wise_matchdist = dist
endif

; Change any NANs to 99.99
tags = tag_names(allobj)
cols = ['gaia_gmag','gaia_gmagerr','tmass_jmag','tmass_jmagerr','tmass_hmag','tmass_hmagerr','tmass_kmag','tmass_kmagerr',$
        'wise_w1mag','wise_w1magerr','wise_w2mag','wise_w2magerr','wise_w3mag','wise_w3magerr','wise_w4mag','wise_w4magerr']
for i=0,n_elements(cols)-1 do begin
  ind = where(tags eq strupcase(cols[i]),nind)
  bd = where(finite(allobj.(ind)) eq 0,nbd)
  if nbd gt 0 then begin
    print,'Fixing ',strtrim(nbd,2),' NANs in ',strupcase(cols[i])
    allobj[bd].(ind) = 99.99
  endif
endfor

; Now save the file
MWRFITS,allobj,outfile,/create
spawn,['gzip',outfile],/noshell

;stop

end
