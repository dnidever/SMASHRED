;+
;
; SMASHRED_CALIBRATE_HEALPIX
;
; This program calibrates all of the photometry for a single HEALPix region
; using transformation equations and ubercal techniques.
;
; INPUTS:
;  pix        The HEALPix pixel number
;  =nside      The HEALPix nside.
;  =version   The version of the final catalogs.  This is only used
;               if OUTPUTDIR is not input.
;  =sumfiles  Input list of summary files to use.  Otherwise smashred_getredinfo
;               uses all summary files for FIELD in dir+'/20??????/*summary.fits'
;  =transfile The file with the photometric transformation equations.
;  =reduxdir  The base directory for the PHOTRED reductions.  The
;               default is "/data/smash/cp/red/photred/"
;  =outputdir The output directory for the catalogs.
;  /usegaia   Use GAIA photometry to set the photometric zeropoint.
;               This is the default
;  /compress  Gzip compress the output FITS files.  This is the default.
;  /redo      Redo a field that was already calibrated.
;  /silent    Don't print anything to the screen.
;
; OUTPUTS:
;  Five binary FITS files are created in the output directory:
;   FIELD_combined_exposures.fits - information on each unique exposure
;   FIELD_combined_chips.fits - information on each unique chip
;   FIELD_combined_allsrc.fits - information on each unique source detection
;   FIELD_combined_allobj.fits - information on each unique object including
;                                     average magnitudes.
;   FIELD_combined_allobj_bright.fits - bright stars in allobj used for
;                                       cross-matching between fields.
;   =error    The error message, if one occurred.
;
; USAGE:
;  IDL>smashred_calibrate_healpix,1001,64
;
; By D.Nidever August 2019, copied from smashred_calibrate_field.pro
;-

pro smashred_calibrate_healpix,pix,nside=nside,version=version,sumfiles=sumfiles,transfile=transfile,reduxdir=reduxdir,outputdir=outputdir,$
                             usegaia=usegaia,redo=redo,compress=compress,silent=silent,error=error

undefine,error
t0 = systime(1)
radeg = 180.0d0 / !dpi

; Not enough inputs
if n_elements(pix) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_calibrate_healpix,pix,nside=nside,version=version,sumfiles=sumfiles,reduxdir=reduxdir,outputdir=outputdir,redo=redo,'
  print,'                                  usegaia=usegaia,compress=compress,error=error'
  return
endif

; Defaults
if n_elements(nside) eq 0 then nside=64
if n_elements(version) eq 0 then version='v6'
if n_elements(reduxdir) eq 0 then reduxdir=SMASHRED_ROOTDIR()+'cp/red/photred/'
;if n_elements(reduxdir) eq 0 then reduxdir='/data/smash/cp/red/photred/'
if file_test(reduxdir,/directory) eq 0 then begin
  error = reduxdir+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif
if n_elements(outputdir) eq 0 then begin
  outputdir = reduxdir+'catalogs/final/'
  if n_elements(version) gt 0 then outputdir+=version+'/'
endif
if file_test(outputdir,/directory) eq 0 then begin
  if not keyword_set(silent) then print,outputdir+' does NOT exist.  Creating it.'
  FILE_MKDIR,outputdir
endif
if n_elements(transfile) eq 0 then transfile=SMASHRED_ROOTDIR()+'cp/red/photred/stdred/smashred_transphot_eqns.fits'
;if n_elements(transfile) eq 0 then transfile='/data/smash/cp/red/photred/stdred/smashred_transphot_eqns.fits'
; Temporary directory
tmpdir = outputdir+'/tmp/'
if file_test(tmpdir,/directory) eq 0 then FILE_MKDIR,tmpdir
; Compression
if n_elements(compress) eq 0 then compress=1
; GAIA photometry
if n_elements(usegaia) eq 0 then usegaia=0

; Check if the final files already exists
outfile = outputdir+strtrim(pix,2)+'_combined'
testfiles = outfile+['_exposures','_chips','_allsrc','_allobj']+'.fits'
ntestfiles = n_elements(testfiles)
if (total(file_test(testfiles)) eq ntestfiles or total(file_test(testfiles+'.gz')) eq ntestfiles) and $
  not keyword_set(redo) then begin
  print,'Final output files already exist for ',strtrim(pix,2),' and /redo NOT set.'
  return
endif


; Get the boundary coordinates
;   healpy.boundaries but not sure how to do it in IDL
;   pix2vec_ring/nest can optionally return vertices but only 4
;     maybe subsample myself between the vectors
; Expand the boundary to include a "buffer" zone
;  to deal with edge cases
;PIX2VEC_RING,nside,pix,vec,vertex

; Use python code to get the boundary
;  this takes ~2s mostly from import statements
tempfile = MKTEMP('bnd')
file_delete,tempfile+'.fits',/allow
step = 100
;; for some reason on thing I have to import matplotlib other healpy
;;    throws an error on import.  no idea why.
pylines = 'python -c "import matplotlib; from healpy import boundaries; from astropy.io import fits;'+$
          ' v=boundaries('+strtrim(nside,2)+','+strtrim(pix,2)+',step='+strtrim(step,2)+');'+$
          " fits.writeto('"+tempfile+".fits'"+',v)"'
spawn,pylines,out,errout
vecbound = MRDFITS(tempfile+'.fits',0,/silent)
file_delete,[tempfile,tempfile+'.fits'],/allow
VEC2ANG,vecbound,theta,phi
rabound = phi*radeg
decbound = 90-theta*radeg

; Expand the boundary by the buffer size
PIX2ANG_RING,nside,pix,centheta,cenphi
cenra = cenphi*radeg
cendec = 90-centheta*radeg
; reproject onto tangent plane
ROTSPHCEN,rabound,decbound,cenra,cendec,lonbound,latbound,/gnomic
; expand by a fraction, it's not an extact boundary but good enough
buffsize = 10.0/3600. ; in deg
radbound = sqrt(lonbound^2+latbound^2)
frac = 1.0 + 1.5*max(buffsize/radbound)
lonbuff = lonbound*frac
latbuff = latbound*frac
ROTSPHCEN,lonbuff,latbuff,cenra,cendec,rabuff,decbuff,/gnomic,/reverse
if range(rabuff) gt 100 then begin  ; deal with RA=0 wraparound
  bd = where(rabuff gt 180,nbd)
  if nbd gt 0 then rabuff[bd] -=360.0
endif
buffer = {cenra:cenra,cendec:cendec,rar:minmax(rabuff),decr:minmax(decbuff),ra:rabuff,dec:decbuff,lon:lonbuff,lat:latbuff,lr:minmax(lonbuff),br:minmax(latbuff)}


; Get reduction info
print,'Getting reduction information'
if n_elements(sumfiles) gt 0 then begin
  SMASHRED_GETREDINFO,info,sumfiles=sumfiles,/silent
endif else begin
  SMASHRED_GETREDINFO,allinfo,/silent
  nallinfo = n_elements(allinfo)
  NEIGHBOURS_RING,nside,pix,neipix,nneipix
  keep = lonarr(nallinfo)
  ;; Find sumfiles that overlap this pixel or neighbors
  for i=0,nallinfo-1 do begin
    fstr1 = *allinfo[i].fstr
    for j=0,n_elements(fstr1)-1 do begin  ;; exposure loop
      ; Calculate the healpix
      x = [0.473, 0.940, 1.100,  1.100,  0.940,  0.473, -0.473, -0.940, -1.100, -1.100, -0.940, -0.473]
      y = [1.000, 0.504, 0.167, -0.167, -0.504, -1.000, -1.000, -0.504, -0.167,  0.167,  0.504,  1.000]
      r = x/cos(fstr1[j].dec/radeg)+fstr1[j].ra
      d = y+fstr1[j].dec
      ROTSPHCEN,r,d,buffer.cenra,buffer.cendec,lon,lat,/gnomic
      olap = dopolygonsoverlap(lon,lat,buffer.lon,buffer.lat)
      if olap eq 1 then keep[i] = 1
      ;ANG2PIX_RING,nside,(90-d)/radeg,r/radeg,ipring      
      ;MATCH,ipring,neipix,ind1,ind2,/sort,count=nmatch
      ;if nmatch gt 0 then keep[i]=1
    endfor
  endfor
  gdinfo = where(keep eq 1,ngdinfo)
  if ngdinfo eq 0 then begin
    error = 'No reduced exposures for '+strtrim(pix,2)
    if not keyword_set(silent) then print,error
    return
  endif
  info = allinfo[gdinfo]
endelse
ninfo = n_elements(info)


; MAYBE GET THE CATALOGS FROM THE PHOTRED DIRECTOREIS DIRECTLY!!!!???
; HOW DO WE KNOW IF THE CATALOGS/PHOT FILES ARE CALIBRATED OR INST??
;  Want instrumental photometry.  Maybe only do instrumental mags with PHOTRED
;  Maybe just us the AST files!!!??  We are going to calibrate the photometry
;  anyway so why not also correct for exptime and apcor as well??!!
;  Another option would be to always include the instrumental phot in
;  the PHOT files and use those instead.

; The various steps for this calibration:
; 1.) Load all of the chip data from .phot files
; 2.) Crossmatch all of the sources and build ALLSRC and ALLOBJ structures
; 3.) Calibrate:  smashred_photcalib.pro
;     a) measure photometric offsets between chip pairs (using calibrated photometry)
;     b) determine relative magnitude offsets per chip using ubercal and
;     c) set zeropoint per band with "photometric" data or 0.9m data
;     d) calculate average photometry per object per filter
;     e) calibrate instrumental photometry with trans eqns. and ubercal mag offsets
;
; HOW DOES THE 0.9M DATA FIT IN????  Can't be done just as an afterburner,
; since this will affect the mean mags which will in turn affect the color
; terms.  So needs to be a 
;
; 4.) Output the data and compress

; Load the transformation equations for all chips, nights and filters
; Load the list of photometric/non-photometric nights
; Load 0.9m data

print,'Running SMASHRED_CALIBRATE_HEALPIX on ',strtrim(pix,2),' with NSIDE=',strtrim(nside,2)
print,strtrim(ninfo,2),' PHOTRED catalogs found'
print,info.file


; Check for temporarily saved crossmatch file
crossmatchfile = tmpdir+strtrim(pix,2)+'_crossmatch.dat'
if file_test(crossmatchfile) eq 1 and not keyword_set(redo) then goto,crossmatch  ; skip loading


; Loop through the catalogs
print,'--------------------------------------------------'
print,'--- STEP 1. Load all of the PHOTRED photometry ---'
print,'=================================================='
ninfo = n_elements(info)
undefine,allfstr,allchstr,allsrc
count_allfstr = 0L
count_allchstr = 0L
count_allsrc = 0L
For c=0,ninfo-1 do begin
  info1 = info[c]
  print,strtrim(c+1,2),'/',strtrim(ninfo,2),' ',info1.file
  undefine,chstr1,allsrc1
  SMASHRED_LOAD_CATPHOT_HEALPIX,pix,info1,chstr1,allsrc1,/useast,reduxdir=reduxdir,redo=redo,outputdir=outputdir,buffer=buffer

  ; Some good data to add
  if n_elements(allsrc1) gt 0 then begin

    ;; Only keep chips that overlap
    norigchstr1 = n_elements(chstr1)
    gdch = where(chstr1.overlap eq 1 and chstr1.nsrc gt 0,ngdch)
    chstr1 = chstr1[gdch]
    ;; Updating the allsrc CHIPINDX with a lookup table
    newchipindx = lonarr(norigchstr1)-1
    newchipindx[gdch] = lindgen(ngdch)
    allsrc1.chipindx = newchipindx[allsrc1.chipindx]

    ;; Remove some exposures
    ui = uniq(chstr1.expnum,sort(chstr1.expnum))
    uexp = chstr1[ui].expnum
    nuexp = n_elements(uexp)
    if nuexp ne n_elements(*info1.fstr) then begin
      print,'Removing some exposures that are not used'
      fstr = *info1.fstr
      MATCH,fstr.expnum,uexp,ind1,ind2,/sort,count=nmatch
      torem = lindgen(n_elements(fstr))
      remove,ind1,torem  ; remove the ones that matched
      REMOVE,torem,fstr
      info1.nexp = n_elements(fstr)
      uibands = uniq(fstr.filter,sort(fstr.filter))
      expbands = strjoin(fstr[uibands].filter)
      info1.bands = expbands
      info1.fstr = ptr_new(fstr)
    endif

    ;; offset indices
    chstr1.allsrcindx += count_allsrc
    allsrc1.chipindx += count_allchstr

    ;; Combine structures
    ;;--------------------
    ;; Create ALL structures
    if n_elements(allfstr) eq 0 then allfstr = make_structure(*info1.fstr,1000,count=nallfstr)
    if n_elements(allchstr) eq 0 then allchstr = make_structure(chstr1,100000L,count=nallchstr)
    if n_elements(allsrc) eq 0 then allsrc = make_structure(allsrc1,1000000L,count=nallsrc)
    ;; Add more elements
    if count_allfstr+n_elements(*info1.fstr) gt nallfstr then allfstr=add_elements(allcstr,n_elements(*info1.fstr)>1000,count=nallfstr)
    if count_allchstr+n_elements(chstr1) gt nallchstr then allchstr=add_elements(allchstr,n_elements(chstr1)>10000L,count=nallchstr)
    if count_allsrc+n_elements(allsrc1) gt nallsrc then allsrc=add_elements(allsrc,n_elements(allsrc1)>100000L,count=nallsrc)
    ;; Put in the new information
    allfstr[count_allfstr:count_allfstr+n_elements(*info1.fstr)-1] = *info1.fstr
    count_allfstr += n_elements(*info1.fstr)
    allchstr[count_allchstr:count_allchstr+n_elements(chstr1)-1] = chstr1
    count_allchstr += n_elements(chstr1)
    allsrc[count_allsrc:count_allsrc+n_elements(allsrc1)-1] = allsrc1
    count_allsrc += n_elements(allsrc1)

    ;push,allfstr,*info1.fstr
    ;push,allchstr,chstr1
    ;push,allsrc,allsrc1
  endif
Endfor ; catalog loop
;; Trim extra elements
if count_allfstr lt nallfstr then allfstr=allfstr[0:count_allfstr-1]
if count_allchstr lt nallchstr then allchstr=allchstr[0:count_allchstr-1]
if count_allsrc lt nallsrc then allsrc=allsrc[0:count_allsrc-1]
; Kludge!  Remove duplicate exposures in short/deep for Field130
;if field eq 'Field130' then begin
;  print,'KLUDGE! Removing duplicate exposures in short/deep for Field130'
;  bd = where((allfstr.expnum eq '00423440' and allfstr.alf_nsources lt 0) or $
;             (allfstr.expnum eq '00426607' and allfstr.alf_nsources lt 0),nbd)
;  REMOVE,bd,allfstr
;endif
; Rename
fstr = allfstr
chstr = allchstr
undefine,allfstr,allchstr

stop

print,'---------------------------------------------------------------------'
print,'--- STEP 2. Crossmatch all of the sources and build ALLSRC/ALLOBJ ---'
print,'====================================================================='
crossmatch:
SMASHRED_CROSSMATCH,strtrim(pix,2),fstr,chstr,allsrc,allobj,reduxdir=reduxdir,redo=redo,outputdir=outputdir

;; Apply the HEALPix boundary cut (without buffer)
ANG2PIX_RING,nside,(90-allobj.dec)/radeg,allobj.ra/radeg,ipring      
gdobj = where(ipring eq pix,ngdobj,comp=bdobj,ncomp=nbdobj)
print,strtrim(ngdobj,2),' fall inside the healpix boundary'
; allobj.srcindx    points to allsrc
; allobj.srcfindx   points to allsrc
; allsrc.fid        is the object name   doesn't matter
; allsrc.cmbindx    points to allobj
; allsrc.chipindx   points to allchstr
; chstr.nsrc        number of sources in this chip
; chstr.allsrcindx  chip sources start in allsrc
;NEED TO UPDATE ALL OF THE ARRAYS IF WE REMOVE OBJECTS
; what if whole exposures are removed?
; need an array that translates old srcindex to new srcindex
; check nsc_instcal_combine.pro at end where it does this same thing

;; Some objects to remove
if nbdobj gt 0 then begin
  print,'Updating structures'
  ;; Get index array of ALLSRC elements to remove
  nbdsrc = long(total(allobj[bdobj].ndet))
  bdsrc = lonarr(nbdsrc)
  cnt = 0L
  for i=0,nbdobj-1 do begin
    ndet = allobj[bdobj[i]].ndet
    ind = allobj[bdobj[i]].srcindx[0:ndet-1]
    bdsrc[cnt:cnt+ndet-1] = ind
    cnt += ndet
  endfor
  gdsrc = lindgen(n_elements(allsrc))
  REMOVE,bdsrc,gdsrc
  ngdsrc = n_elements(gdsrc)
  ;; Update NSRC in CHSTR
  for i=0,nbdsrc-1 do chstr[allsrc[bdsrc[i]].chipindx].nsrc--
  ;; Chips we need to remove
  bdch = where(chstr.nsrc le 0,nbdch,comp=gdch,ncomp=ngdch)
  if nbdch gt 0 then begin
    ;; Indices into the new CHSTR structure using old index
    newchstrindx = lonarr(n_elements(chstr))-1
    newchstrindx[gdch] = lindgen(ngdch)
    ;; Update CHIPINDX in ALLSRC, only for rows we are keeping
    allsrc[gdsrc].chipindx = newchstrindx[allsrc[gdsrc].chipindx]
    ;; Update NCHIPS in FSTR
    for k=0,nbdch-1 do begin
      MATCH,chstr[bdch[k]].expnum,fstr.expnum,ind1,ind2,/sort,count=nmatch
      fstr[ind2].nchips--
    endfor
  endif
  ;; Update CMBINDX in ALLSRC
  ;;   Indices into the new ALLOBJ using old index
  newallobjindx = lonarr(n_elements(allobj))-1
  newallobjindx[gdobj] = lindgen(ngdobj)
  allsrc[gdsrc].cmbindx = newallobjindx[allsrc[gdsrc].cmbindx]
  ;; Update ALLSRCINDX in CHSTR
  sich = sort(chstr.allsrcindx)  ; chstr are not sorted
  newallsrcindx = [0,long(total(chstr[sich].nsrc,/cum))]
  newallsrcindx = newallsrcindx[0:n_elements(chstr)-1]
  chstr[sich].allsrcindx = newallsrcindx
  ;; Update SRCINDX/SRCFINDX in ALLOBJ
  ;;   srcindx
  srcindx = allobj.srcindx
  newallsrcindx = lonarr(n_elements(allsrc))-1
  newallsrcindx[gdsrc] = lindgen(ngdsrc)
  g = where(srcindx ge 0,ng)
  srcindx[g] = newallsrcindx[srcindx[g]]
  allobj.srcindx = srcindx
  ;;  srcfindx
  srcfindx = allobj.srcfindx
  g = where(srcfindx ge 0,ng)
  srcfindx[g] = newallsrcindx[srcfindx[g]]
  allobj.srcfindx = srcfindx
  ;; Exposures to remove
  bdfstr = where(fstr.nchips le 0,nbdfstr,comp=gdfstr,ncomp=ngdfstr)
  ;; Shifting SRCFINDX rows in ALLOBJ
  if nbdfstr gt 0 then begin
    ;; shift them over
    oldsrcfindx = allobj.srcfindx
    allobj.srcfindx = -1
    for k=0,ngdfstr-1 do allobj.srcfindx[k]=oldsrcfindx[gdfstr[k],*]
    undefine,oldsrcfindx
  endif
  ;; Remove the rows from all structures
  REMOVE,bdobj,allobj
  REMOVE,bdsrc,allsrc
  if nbdch gt 0 then REMOVE,bdch,chstr
  if nbdfstr gt 0 then REMOVE,bdfstr,fstr
endif

print,'-----------------------------------------------'
print,'--- STEP 3. Calibrate all of the photometry ---'
print,'==============================================='
info[0].field = strtrim(pix,2)
SMASHRED_PHOTCALIB,info,fstr,chstr,allsrc,allobj,transfile=transfile,usegaia=usegaia,reduxdir=reduxdir,$
                   outputdir=outputdir


; Compute average morphology and coordinate values
print,'Calculating average morphology and coordinate parameters'
SMASHRED_AVERAGEMORPHCOORD,fstr,chstr,allsrc,allobj

; Compute exposure map
;print,'Computing the exposure maps'
;SMASHRED_COMPUTE_EXPMAP,strtrim(pix,2),chstr,redo=redo,outputdir=outputdir

; Set non-detections based on the exposure map
;print,'Setting non-detections based on the exposure maps'
;SMASHRED_SET_NONDETECTIONS,strtrim(pix,2),allobj,dir=outputdir

; Calculate extinction
print,'Getting SFD E(B-V) extinctions'
GLACTC,allobj.ra,allobj.dec,2000.0,lon,lat,1,/deg
ebv = DUST_GETVAL(lon,lat,/interp,/noloop)
; add EBV, G0 and I0 to the catalogs
allobj.ebv = ebv

; Write out the final files
print,'Writing combined file to ',outfile
MWRFITS,fstr,outfile+'_exposures.fits',/create
MWRFITS,chstr,outfile+'_chips.fits',/create
MWRFITS,allsrc,outfile+'_allsrc.fits',/create
MWRFITS,allobj,outfile+'_allobj.fits',/create
; Compress
if keyword_set(compress) then begin
  print,'Compressing output files'
  spawn,['gzip','-f',outfile+'_exposures.fits'],out,errout,/noshell
  spawn,['gzip','-f',outfile+'_chips.fits'],out,errout,/noshell
  spawn,['gzip','-f',outfile+'_allsrc.fits'],out,errout,/noshell
  spawn,['gzip','-f',outfile+'_allobj.fits'],out,errout,/noshell
  ;spawn,['gzip','-f',outfile+'_expmap.fits'],out,errout,/noshell
endif

; Make bright allobj catalog
SMASHRED_MAKE_BRIGHTCAT,strtrim(pix,2),redo=redo,dir=outputdir

; Make DEEP allobj catalog
;   redo the average photometry, morphology parameters and coordinates
SMASHRED_AVERAGEPHOT,fstr,chstr,allsrc,allobj,/usecalib,/deeponly
SMASHRED_AVERAGEMORPHCOORD,fstr,chstr,allsrc,allobj,/deeponly
deepoutfile = outfile+'_allobj_deep.fits'
print,'Writing deep to ',deepoutfile
file_delete,[deepoutfile,deepoutfile+'.gz'],/allow
MWRFITS,allobj,deepoutfile,/create
spawn,['gzip',deepoutfile],out,errout,/noshell

;; Crossmatching with other catalogs
;SMASHRED_MATCHCATS_GAIADR2
;; Do smashred_matchcats_gaiadr2.pro


; Print processing time
dt = systime(1)-t0
print,'' & print,'Processing time = ',strtrim(string(dt/60.0,format='(F20.3)'),2),' min.'

;stop

end
