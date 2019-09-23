pro make_field_table

  ;; Make the SMASH field table for DR2

  smash = importascii('~/projects/SMASHRED/data/smash_fields_final.txt',/header)
  calib = mrdfits('/dl1/users/dnidever/smash/cp/red/photred/check_calibrated_v6.fits',1)
  calib.field = strtrim(calib.field,2)
  ;dr1 = mrdfits('~/smash/reduction/smash_dr1_field_table.fits',1)
  chstr = mrdfits('/dl1/users/dnidever/smash/cp/red/photred/smash_chips.fits.gz',1)
  chstr.smash_field = strtrim(chstr.smash_field,2)
  
  schema = {name:'',fieldid:0L,ra:0.0d0,dec:0.0d0,glon:0.0d0,glat:0.0d0,mslon:0.0d0,mslat:0.0d0,$
            ring256:0L,nchips:0L,nsrc:0L,nobj:0L,nexp:0L,nexp_u:0L,nexp_g:0L,nexp_r:0L,nexp_i:0L,$
            nexp_z:0L,ucalib:0L,gcalib:0L,rcalib:0L,icalib:0L,zcalib:0L}
  tags = tag_names(schema)
  cat = replicate(schema,197)
  struct_assign,calib,cat
  cat.name = calib.field
  cat.fieldid = long(strmid(cat.name,5))
  glactc,cat.ra,cat.dec,2000.0,glon,glat,1,/deg
  cat.glon = glon
  cat.glat = glat
  gal2mag,glon,glat,mlon,mlat
  cat.mslon = mlon
  cat.mslat = mlat
  cat.ucalib = calib.nucalib
  cat.gcalib = calib.ngcalib
  cat.rcalib = calib.nrcalib
  cat.icalib = calib.nicalib
  cat.zcalib = calib.nzcalib

  chindex = create_index(chstr.smash_field)
  MATCH,cat.name,chindex.value,ind1,ind2,/sort,count=nmatch
  if nmatch ne 197 then stop,'not all match'
  cat[ind1].nchips = chindex.num[ind2]
  filters = ['u','g','r','i','z']
  for i=0,nmatch-1 do begin
    ind = chindex.index[chindex.lo[ind2[i]]:chindex.hi[ind2[i]]]
    chstr1 = chstr[ind]
    cat[ind1[i]].nsrc = long(total(chstr1.nsrc))
    for j=0,4 do begin
      gd = where(chstr1.filter eq filters[j],ngd)
      ui = uniq(chstr1[gd].expnum,sort(chstr1[gd].expnum))
      nexp = n_elements(ui)
      colind = where(tags eq 'NEXP_'+strupcase(filters[j]),ncolind)
      cat[ind1[i]].(colind[0]) = nexp
    endfor
    ui = uniq(chstr1.expnum,sort(chstr1.expnum))
    cat[ind1[i]].nexp = n_elements(ui)
  endfor

  ang2pix_ring,256,(90-cat.dec)/!radeg,cat.ra/!radeg,pix
  cat.ring256 = pix

  ;; Get number of objects from FIELDXX_combined_allobj.fits.gz files
  for i=0,n_elements(cat)-1 do begin
    file = '/dl1/users/dnidever/smash/cp/red/photred/catalogs/final/v6/'+cat[i].name+'_combined_allobj.fits.gz'
    hd = headfits(file,exten=1)
    cat[i].nobj = long(sxpar(hd,'NAXIS2'))
  endfor

  ;MWRFITS,cat,'/dl1/users/dnidever/smash/cp/red/photred/catalogs/dr2/smash_dr2_field_table.fits',/create
  
  stop

  end
