pro smashred_matchcats_gaiadr2,field,version=version,redo=redo

; Crossmatch SMASH catalog with GAIA DR2, 2MASS and ALLWISE

; Not enough inputs
if n_elements(field) eq 0 then begin
  print,'Syntax - smashred_matchcats_gaiadr2,field,version=version,redo=redo'
  return
endif

if n_elements(version) eq 0 then version='v6'
dir = smashred_rootdir()+'cp/red/photred/catalogs/final/'+version+'/'
;dir = '/data/smash/cp/red/photred/catalogs/final/'+version+'/'

print,'Running SMASHRED_MATCHCATS on ',field,' version=',version

if file_test(dir+field+'_combined_allobj.fits.gz') eq 0 then begin
  print,'ALLOBJ file not found for ',field,' and version ',version
  return
endif

outfile = dir+field+'_combined_allobj_xmatch2.fits'
if (file_test(outfile) eq 1 or file_test(outfile+'.gz') eq 1) and not keyword_set(redo) then begin
  print,'Output file exists and /redo NOT set'
  return
endif

; Restore the SMASH catalog
allobj = MRDFITS(dir+field+'_combined_allobj.fits.gz',1)
nallobj = n_elements(allobj)

; Copy to new schema
nan = !values.f_nan
xtra_schema = {id:'',gaia_match:0B,gaia_matchdist:9999.0,gaia_source:-1L,gaia_ra:9999.0d0,gaia_dec:9999.0d0,$
               gaia_raerr:9999.0,gaia_decerr:9999.0,gaia_parallax:9999.0,gaia_parallax_error:9999.0,$
               gaia_pmra:9999.0,gaia_pmra_error:9999.9,gaia_pmdec:9999.0,gaia_pmdec_error:9999.0,$
               gaia_gmag:9999.0,gaia_gerr:9999.0,gaia_bpmag:9999.0,gaia_bperr:9999.0,gaia_rpmag:9999.0,gaia_rperr:9999.0,$
               tmass_match:0B,tmass_matchdist:9999.0,tmass_id:'',tmass_ra:9999.0d0,tmass_dec:9999.0d0,$
               tmass_raerr:9999.0,tmass_decerr:9999.0,tmass_jmag:9999.0,tmass_jerr:9999.0,$
               tmass_hmag:9999.0,tmass_herr:9999.0,tmass_kmag:9999.0,tmass_kerr:9999.0,tmass_qflg:'',$
               wise_match:0B,wise_matchdist:9999.0,wise_id:'',wise_ra:9999.0d0,wise_dec:9999.0d0,$
               wise_raerr:9999.0,wise_decerr:9999.0,wise_w1mag:9999.0,wise_w1err:9999.0,$
               wise_w2mag:9999.0,wise_w2err:9999.0,wise_w3mag:9999.0,wise_w3err:9999.0,$
               wise_w4mag:9999.0,wise_w4err:9999.0,wise_qph:''}
xtra = REPLICATE(xtra_schema,nallobj)
xtra.id = allobj.id

; Restore the GAIADR2 file
print,'Matching with GAIADR2'
gaiadir = smashred_rootdir()+'cp/red/photred/gaiadr2/'
;gaiadir = '/data/smash/cp/red/photred/gaia/'
if file_test(gaiadir+field+'_gaiadr2.fits.gz') eq 0 then begin
  print,'No GAIA file found for field ',field
  return
endif
gaia = MRDFITS(gaiadir+field+'_gaiadr2.fits.gz',1)
SRCMATCH,allobj.ra,allobj.dec,gaia.ra,gaia.dec,0.5,ind1,ind2,/sph,count=nmatch
print,strtrim(nmatch,2),' GAIADR2 matches'
if nmatch gt 0 then begin
  xtra[ind1].gaia_match = 1
  xtra[ind1].gaia_source = gaia[ind2].source_id
  xtra[ind1].gaia_ra = gaia[ind2].ra
  xtra[ind1].gaia_dec = gaia[ind2].dec
  xtra[ind1].gaia_raerr = gaia[ind2].ra_error
  xtra[ind1].gaia_decerr = gaia[ind2].dec_error
  xtra[ind1].gaia_parallax = gaia[ind2].parallax
  xtra[ind1].gaia_parallax_error = gaia[ind2].parallax_error
  xtra[ind1].gaia_pmra = gaia[ind2].pmra
  xtra[ind1].gaia_pmra_error = gaia[ind2].pmra_error
  xtra[ind1].gaia_pmdec = gaia[ind2].pmdec
  xtra[ind1].gaia_pmdec_error = gaia[ind2].pmdec_error
  xtra[ind1].gaia_gmag = gaia[ind2].phot_g_mean_mag
  xtra[ind1].gaia_gerr = 2.5*alog10(1.0+gaia[ind2].phot_g_mean_flux_error/gaia[ind2].phot_g_mean_flux)
  xtra[ind1].gaia_bpmag = gaia[ind2].phot_bp_mean_mag
  xtra[ind1].gaia_bperr = 2.5*alog10(1.0+gaia[ind2].phot_bp_mean_flux_error/gaia[ind2].phot_bp_mean_flux)
  xtra[ind1].gaia_rpmag = gaia[ind2].phot_rp_mean_mag
  xtra[ind1].gaia_rperr = 2.5*alog10(1.0+gaia[ind2].phot_rp_mean_flux_error/gaia[ind2].phot_rp_mean_flux)
  dist = sphdist(allobj[ind1].ra,allobj[ind1].dec,gaia[ind2].ra,gaia[ind2].dec,/deg)*3600.0d0
  xtra[ind1].gaia_matchdist = dist
endif

; Restore the 2MASS file
print,'Matching with 2MASS'
tmassdir = smashred_rootdir()+'cp/red/photred/tmass/'
;tmassdir = '/data/smash/cp/red/photred/tmass/'
if file_test(tmassdir+field+'_tmass.fits.gz') eq 0 then begin
  print,'No 2MASS file found for field ',field
  return
endif
tmass = MRDFITS(tmassdir+field+'_tmass.fits.gz',1)
SRCMATCH,allobj.ra,allobj.dec,tmass.raj2000,tmass.dej2000,0.5,ind1,ind2,/sph,count=nmatch
print,strtrim(nmatch,2),' 2MASS matches'
if nmatch gt 0 then begin
  xtra[ind1].tmass_match = 1
  xtra[ind1].tmass_id = tmass[ind2]._2mass
  xtra[ind1].tmass_ra = tmass[ind2].raj2000
  xtra[ind1].tmass_dec = tmass[ind2].dej2000
  xtra[ind1].tmass_raerr = tmass[ind2].errmaj
  xtra[ind1].tmass_decerr = tmass[ind2].errmin
  xtra[ind1].tmass_jmag = tmass[ind2].jmag
  xtra[ind1].tmass_jerr = tmass[ind2].e_jmag
  xtra[ind1].tmass_hmag = tmass[ind2].hmag
  xtra[ind1].tmass_herr = tmass[ind2].e_hmag
  xtra[ind1].tmass_kmag = tmass[ind2].kmag
  xtra[ind1].tmass_kerr = tmass[ind2].e_kmag
  xtra[ind1].tmass_qflg = tmass[ind2].qflg
  dist = sphdist(allobj[ind1].ra,allobj[ind1].dec,tmass[ind2].raj2000,tmass[ind2].dej2000,/deg)*3600.0d0
  xtra[ind1].tmass_matchdist = dist
endif

; Restore the WISE file
print,'Matching with WISE'
wisedir = smashred_rootdir()+'cp/red/photred/wise/'
;wisedir = '/data/smash/cp/red/photred/wise/'
if file_test(wisedir+field+'_wise.fits.gz') eq 0 then begin
  print,'No WISE file found for field ',field
  return
endif
wise = MRDFITS(wisedir+field+'_wise.fits.gz',1)
SRCMATCH,allobj.ra,allobj.dec,wise.ra,wise.dec,0.5,ind1,ind2,/sph,count=nmatch
print,strtrim(nmatch,2),' WISE matches'
if nmatch gt 0 then begin
  xtra[ind1].wise_match = 1
  ;; from queryvizier
  if tag_exist(wise,'RAJ2000') then begin
    xtra[ind1].wise_id = wise[ind2].allwise
    xtra[ind1].wise_ra = wise[ind2].raj2000
    xtra[ind1].wise_dec = wise[ind2].dej2000
    xtra[ind1].wise_raerr = wise[ind2].eemaj
    xtra[ind1].wise_decerr = wise[ind2].eemin
    xtra[ind1].wise_w1mag = wise[ind2].w1mag
    xtra[ind1].wise_w1err = wise[ind2].e_w1mag
    xtra[ind1].wise_w2mag = wise[ind2].w2mag
    xtra[ind1].wise_w2err = wise[ind2].e_w2mag
    xtra[ind1].wise_w3mag = wise[ind2].w3mag
    xtra[ind1].wise_w3err = wise[ind2].e_w3mag
    xtra[ind1].wise_w4mag = wise[ind2].w4mag
    xtra[ind1].wise_w4err = wise[ind2].e_w4mag
    xtra[ind1].wise_qph = wise[ind2].qph
    dist = sphdist(allobj[ind1].ra,allobj[ind1].dec,wise[ind2].raj2000,wise[ind2].dej2000,/deg)*3600.0d0
    xtra[ind1].wise_matchdist = dist
  ;; some datalab
  endif else begin
    xtra[ind1].wise_id = wise[ind2].designation
    xtra[ind1].wise_ra = wise[ind2].ra
    xtra[ind1].wise_dec = wise[ind2].dec
    xtra[ind1].wise_raerr = wise[ind2].sigra
    xtra[ind1].wise_decerr = wise[ind2].sigdec
    xtra[ind1].wise_w1mag = wise[ind2].w1mpro
    xtra[ind1].wise_w1err = wise[ind2].w1sigmpro
    xtra[ind1].wise_w2mag = wise[ind2].w2mpro
    xtra[ind1].wise_w2err = wise[ind2].w2sigmpro
    xtra[ind1].wise_w3mag = wise[ind2].w3mpro
    xtra[ind1].wise_w3err = wise[ind2].w3sigmpro
    xtra[ind1].wise_w4mag = wise[ind2].w4mpro
    xtra[ind1].wise_w4err = wise[ind2].w4sigmpro
    xtra[ind1].wise_qph = wise[ind2].ph_qual
    dist = sphdist(allobj[ind1].ra,allobj[ind1].dec,wise[ind2].ra,wise[ind2].dec,/deg)*3600.0d0
    xtra[ind1].wise_matchdist = dist
  endelse
endif

; Change any NANs to 9999.0
tags = tag_names(xtra)
cols = ['gaia_gmag','gaia_gerr','tmass_jmag','tmass_jerr','tmass_hmag','tmass_herr','tmass_kmag','tmass_kerr',$
        'wise_w1mag','wise_w1err','wise_w2mag','wise_w2err','wise_w3mag','wise_w3err','wise_w4mag','wise_w4err']
for i=0,n_elements(cols)-1 do begin
  ind = where(tags eq strupcase(cols[i]),nind)
  bd = where(finite(xtra.(ind)) eq 0,nbd)
  if nbd gt 0 then begin
    ;print,'Fixing ',strtrim(nbd,2),' NANs in ',strupcase(cols[i])
    xtra[bd].(ind) = 9999.0
  endif
endfor

; Only keep ones with some matches
gd = where(xtra.gaia_match eq 1 or xtra.tmass_match eq 1 or xtra.wise_match eq 1,ngd)
xtra = xtra[gd]

; Now save the file
print,'Writing results to ',outfile
if file_test(outfile) eq 1 then file_delete,outfile
MWRFITS,xtra,outfile,/create
print,'Compressing'
if file_test(outfile+'.gz') eq 1 then file_delete,outfile+'.gz'
spawn,['gzip',outfile],/noshell

;stop

end
