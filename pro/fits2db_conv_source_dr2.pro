pro fits2db_conv_source_dr2,ifield,redo=redo

  ;; Convert source table for database

  version = 'v6'
  dir = '/dl1/users/dnidever/smash/cp/red/photred/catalogs/final/'+version+'/'
  outdir = dir+'db/'
  ;outdir = '/data0/dnidever/smash/dr2/'
  if file_test(outdir,/directory) eq 0 then file_mkdir,outdir

  fieldid = ifield
  if strmid(ifield,0,5) eq 'Field' then fieldid = fix(strmid(strtrim(ifield,2),5))
  print,ifield

  outfile = outdir+ifield+'_source.fits'
  if file_test(outfile+'.gz') eq 1 and not keyword_set(redo) then begin
    print,outfile+'.gz exists and /redo not set'
    return
  endif

  ;SOURCE table
  ;id -> origid
  ;
  ;FID is the object ID
  ;ID original ID in ALF/ALF files
  ;sourceid???

  ; ALLSRC structure schema
  ;  cmbindx       index into ALLOBJ, added later on in smashred_crossmatch.pro
  ;  fid           final ID, added later on in smashred_crossmatch.pro
  ;  id            original ID in ALS/ALF file
  ;  idref         ID in PHOT/AST file
  ;It looks like there's no UNIQUE SOURCEID column.  ID is the id/number in that chip.alf file.
  ;Would need something like EXPNUM_CHIP.NUMBER
  ;There's also no way to get the FILTER, EXPTIME, MJD information for the source table.
  ;DR1 had:
  ;-fieldid
  ;-expnum
  ;-chip
  ;-mjd
  ;-filter
  ;should also maybe add some kind of chipID, actually does the CHIP table have a unique ID?

  
  allsrc = mrdfits(dir+ifield+'_combined_allsrc.fits.gz',1)
  nallsrc = n_elements(allsrc)
  chips = mrdfits(dir+ifield+'_combined_chips.fits.gz',1)
  nchips = n_elements(chips)

  schema_allsrc = {sourceid:'',id:'',origid:0L,refid:0L,fieldid:0L,expnum:'',chip:0,chipid:'',$
                   mjd:0.0d0,filter:'',x:0.0,y:0.0,xref:0.0,yref:0.0,forced:0b,mag:0.0,err:0.0,$
                   cmag:0.0,cerr:0.0,chi:0.0,sharp:0.0,flag:0,prob:0.0,ra:0.0d0,dec:0.0d0,$
                   raerr:0.0,decerr:0.0,raindiv:0.0d0,decindiv:0.0d0,raref:0.0d0,decref:0.0d0}
  newsrc = replicate(schema_allsrc,nallsrc)
  STRUCT_ASSIGN,allsrc,newsrc
  newsrc.chipid = strtrim(chips[allsrc.chipindx].expnum,2)+'_'+strtrim(string(chips[allsrc.chipindx].chip,format='(i02)'),2)
  newsrc.sourceid = newsrc.chipid+'.'+strtrim(allsrc.id,2)
  if strmid(ifield,0,5) eq 'Field' then newsrc.id = strmid(strtrim(allsrc.fid,2),5) else $
    newsrc.id=strtrim(allsrc.fid,2)  ; strip 'Field' portion
  newsrc.fieldid = fieldid
  newsrc.origid = allsrc.id         ; original allstar/allframe ID
  newsrc.refid = allsrc.idref
  newsrc.expnum = chips[allsrc.chipindx].expnum
  newsrc.chip = chips[allsrc.chipindx].chip
  chips_mjd = dblarr(nchips)
  for j=0,nchips-1 do chips_mjd[j]=date2jd(chips[j].utdate+'T'+chips[j].uttime,/mjd)  
  newsrc.mjd = chips_mjd[allsrc.chipindx]
  newsrc.filter = chips[allsrc.chipindx].filter

  print,'Writing to ',outfile
  MWRFITS,newsrc,outfile,/create
  spawn,['gzip','-f',outfile],/noshell
  undefine,allsrc    ; free up memory
  
  ;stop

  end
