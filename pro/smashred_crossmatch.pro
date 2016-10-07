;+
;
; SMASHRED_CROSSMATCH
;
; This program crossmatches sources for a field and creates
; the ALLSTR and ALLOBJ structures.
;
; INPUTS:
;  field      The name of the field.
;  fstr       The structure with information for each exposure.
;  chstr      The structure with information for each chip.
;  allsrc     The structure with information for each source detection.
;  =dcr       The matching radius in arcsec.  The default is 0.5 arcsec.
;  =reduxdir  The reduction directory, the default is "/data/smash/cp/red/photred/"
;  =outputdir  The output directory, the default is reduxdir+"catalogs/final/"
;  /redo      Reget the photometry even if the temporary output file                                                   
;                already exists.  
;  /silent    Don't print anything to the screen.
;
; OUTPUTS:
;  allobj  The structure with information for each unique object.
;  The CHSTR structure is also updated with the ALLSRCINDX column.
;  =error  The error message if one occurred.
;
; USAGE:
;  IDL>smashred_crossmatch,fstr,chstr,allstr,allobj
; 
; By D.Nidever  March 2016
;-

pro smashred_crossmatch,field,fstr,chstr,allsrc,allobj,dcr=dcr,error=error,reduxdir=reduxdir,$
                        redo=redo,silent=silent,outputdir=outputdir

undefine,allobj

; Not enough inputs
if n_elements(field) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_crossmatch,field'
  return
endif

; Defaults
if n_elements(reduxdir) eq 0 then reduxdir='/data/smash/cp/red/photred/'
if file_test(reduxdir,/directory) eq 0 then begin
  error = reduxdir+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif
if n_elements(outputdir) eq 0 then outputdir=reduxdir+'catalogs/final/'
if file_test(outputdir,/directory) eq 0 then begin
  if not keyword_set(silent) then print,outputdir+' does NOT exist.  Creating it.'
  FILE_MKDIR,outputdir
endif
; Temporary directory
tmpdir = outputdir+'/tmp/'
if file_test(tmpdir,/directory) eq 0 then FILE_MKDIR,tmpdir
; Default matching radius
if n_elements(dcr) eq 0 then dcr=0.5  ; 0.7


; Output filename
outfile = tmpdir+field+'_crossmatch.dat'


; Crossmatch the sources
If file_test(outfile) eq 0 or keyword_set(redo) then begin

  ; Not enough inputs if running for the first time
  if n_elements(fstr) eq 0 or n_elements(chstr) eq 0 or n_elements(allsrc) eq 0 then begin
    error = 'Not enough inputs if running for the first time'
    print,'Need FSTR, CHSTR, ALLSRC if running for the first time.'
    return
  endif
  nfstr = n_elements(fstr)
  nchstr = n_elements(chstr)


  ; Unique exposures from CHSTR
  uiexp = uniq(chstr.expnum,sort(chstr.expnum))
  uexp = chstr[uiexp].expnum
  nuexp = n_elements(uexp)
  ; Unique exposures from FSTR
  uifexp = uniq(fstr.expnum,sort(fstr.expnum))
  ufexp = fstr[uifexp].expnum
  nufexp = n_elements(ufexp)
  ; unqiue FSTR and CHSTR exposures don't match
  if nuexp ne nfstr then begin
    error = 'Exposures do NOT match in FSTR and CHSTR'
    if not keyword_set(silent) then print,error
    return
  endif
  ; duplicate exposures in FSTR
  if nufexp ne nfstr then begin
    error = 'Duplicate exposures in FSTR'
    if not keyword_set(silent) then print,error
    return
  endif


  nan = !values.f_nan
  dnan = !values.d_nan

  ; Make the ALLOBJ structure schema
  allobj_schema = {id:'',ra:0.0d0,dec:0.0d0,rascatter:99.99,decscatter:99.99,ndet:0,depthflag:0B,srcindx:lonarr(nuexp)-1,$
                   srcfindx:lonarr(nuexp)-1,u:99.99,uerr:9.99,uscatter:99.99,ndetu:0,g:99.99,gerr:9.99,gscatter:99.99,ndetg:0,$
                   r:99.99,rerr:9.99,rscatter:99.9,ndetr:0,i:99.99,ierr:9.99,iscatter:99.99,ndeti:0,z:99.99,zerr:9.99,$
                   zscatter:99.99,ndetz:0,chi:nan,sharp:nan,flag:-1,prob:nan,ebv:99.99}
  ur_allsrc_indx = 0LL         ; next one starts from HERE
  nallsrc = n_elements(allsrc)  ; current number of total Allsrc elements, NOT all filled
  ; SRCINDX has NDET indices at the front of the array
  ; SRCFINDX has them in the element that matches the frame they were
  ; detected in
  allobjtags = tag_names(allobj_schema)
  lallobjtags = strlowcase(allobjtags)
  ; Loop through the exposures
  for i=0,nfstr-1 do begin
    ; CHSTR indices for this exposure
    expind = where(chstr.expnum eq fstr[i].expnum,nexpind)
    if fstr[i].exptime lt 100 then depthbit=1 else depthbit=2  ; short or long

    print,strtrim(i+1,2),'/',strtrim(nfstr,2),' adding exposure ',fstr[i].expnum

    ; Get all chip ALLSRC information for this exposure
    undefine,expnew,expnewallsrcindx
    for j=0,nexpind-1 do begin
      ; Get the source data from CHSTR and ALLSRC
      temp_allsrcindx = lindgen(chstr[expind[j]].nsrc)+chstr[expind[j]].allsrcindx
      temp = allsrc[temp_allsrcindx]
      push,expnew,temp         ; add to exposure "new" structure
      push,expnewallsrcindx,temp_allsrcindx  ; allsrc index for this exposure
    endfor

    ;------------------------------
    ; PUT IN ALLOBJ MERGED CATALOG
    ;------------------------------
    ; Copy to new structure type
    If i eq 0 then begin
      allobj = replicate(allobj_schema,n_elements(expnew))
      allobj.id = field+'.'+strtrim(lindgen(n_elements(expnew))+1,2)
      allobj.ra = expnew.ra
      allobj.dec = expnew.dec

      ; Put ALLSRC index in ALLOBJ
      ;allobj.srcindx[0] = lindgen(n_elements(expnew))
      ;allobj.srcfindx[i] = lindgen(n_elements(expnew))
      allobj.srcindx[0] = expnewallsrcindx
      allobj.srcfindx[i] = expnewallsrcindx
      allobj.ndet = 1
      allobj.depthflag OR= depthbit             ; OR combine to depthflag, 1-short, 2-long, 3-short+long
      ; Put ID, CMBINDX in ALLSRC
      ;allsrc[0:cur_allsrc_indx-1].fid = allobj.id
      ;allsrc[0:cur_allsrc_indx-1].cmbindx = lindgen(n_elements(allobj))
      allsrc[expnewallsrcindx].fid = allobj.id
      allsrc[expnewallsrcindx].cmbindx = lindgen(n_elements(allobj))

    ; 2nd and later exposures, check for repeats/overlap
    Endif else begin

      ; Match sources
t0 = systime(1)
      SRCMATCH,allobj.ra,allobj.dec,expnew.ra,expnew.dec,dcr,ind1,ind2,count=nmatch,/sph,/usehist  ; use faster histogram_nd method
;print,'dt=',systime(1)-t0,' sec.  matching time'
      print,' ',strtrim(nmatch,2),' matched sources'
      ; Some matches, add data to existing record for these sources
      if nmatch gt 0 then begin
        for k=0LL,nmatch-1 do allobj[ind1[k]].srcindx[allobj[ind1[k]].ndet] = expnewallsrcindx[ind2[k]]     ; put SRCINDX in ALLOBJ
        allobj[ind1].srcfindx[i] = expnewallsrcindx[ind2]    ; put SRCFINDX in ALLOBJ
        allobj[ind1].ndet++
        allobj[ind1].depthflag OR= depthbit                  ; OR combine to depthflag  
        allsrc[expnewallsrcindx[ind2]].fid = allobj[ind1].id  ; put ID, CMBINDX in ALLSRC
        allsrc[expnewallsrcindx[ind2]].cmbindx = ind1
        ; Remove stars
        if nmatch lt n_elements(expnew) then remove,ind2,expnew,expnewallsrcindx else undefine,expnew,expnewallsrcindx
      endif

      ; Some left, add records for these sources
      if n_elements(expnew) gt 0 then begin
        print,' ',strtrim(n_elements(expnew),2),' sources left to add'
        newallobj = replicate(allobj_schema,n_elements(expnew))
        newallobj.id = field+'.'+strtrim(lindgen(n_elements(expnew))+1+n_elements(allobj),2)
        newallobj.ra = expnew.ra
        newallobj.dec = expnew.dec
        newallobj.srcindx[0] = expnewallsrcindx      ; put SRCINDX in ALLOBJ
        newallobj.srcfindx[i] = expnewallsrcindx     ; put SRCFINDX in ALLOBJ
        newallobj.ndet = 1
        newallobj.depthflag OR= depthbit             ; OR combine to depthflag
        allsrc[expnewallsrcindx].fid = newallobj.id  ; put ID, CMBINDX in ALLSRC
        allsrc[expnewallsrcindx].cmbindx = lindgen(n_elements(expnew))+n_elements(allobj)
        ; concatenating two large structures causes lots to be zerod out
;t0 = systime(1)
        nold = n_elements(allobj)
        nnew = n_elements(newallobj)
        new = replicate(allobj_schema,nold+nnew)
        new[0:nold-1] = allobj
        new[nold:*] = newallobj
        allobj = new
        undefine,new
;print,'dt=',systime(1)-t0
      endif

    Endelse  ; 2nd or later

bd = where(allobj.ra eq 0.0,nbd)
if nbd gt 0 then stop,'ALLOBJ Zerod out elements problem!!!'

    ;stop

  Endfor  ; exposure loop

  ; Saving the results to output file
  print,'Saving information temporarily to ',outfile
  SAVE,fstr,chstr,allsrc,allobj,file=outfile


; Load the existing temporary file                                                                                                                                                  
Endif else begin
  print,'Restoring previously saved ',outfile
  RESTORE,outfile
Endelse

;stop

end
