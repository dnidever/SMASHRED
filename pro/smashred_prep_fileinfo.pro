pro smashred_prep_fileinfo,fstr,files=files,silent=silent,image=image

; Get "calibrated" images ready for PHOTRED

undefine,fstr

nfiles = n_elements(files)
if nfiles eq 0 then begin
  undefine,files
  ;files1 = file_search('c4d_*_ooi_*.fits*',count=nfiles1)
  searchstring = 'c4d_*_*.fits*'
  if keyword_set(image) then searchstring = 'c4d_*_*_ooi*.fits*'
  files1 = file_search(searchstring,count=nfiles1)
  if nfiles1 gt 0 then push,files,files1
  files2 = file_search('tu*.fits*',count=nfiles2)
  if nfiles2 gt 0 then push,files,files2
  nfiles = n_elements(files)
  ; No archive files, get chip 1 split files, 00317720_01.fits
  if nfiles eq 0 then begin
    print,'NO ARCHIVE FILES.  Getting Chip 1 SPLIT files'
    undefine,files
    files1 = file_search('????????_01.fits',count=nfiles1)
    if nfiles1 gt 0 then push,files,files1
    files2 = file_search('????????c_01.fits',count=nfiles2)  ; combined files
    if nfiles2 gt 0 then push,files,files2
    nfiles = n_elements(files)
  endif
endif
print,strtrim(nfiles,2),' images'
if nfiles eq 0 then return

; Get information on all the images
fstr = replicate({file:'',base:'',archivename:'',object:'',ra:'',dec:'',dateobs:'',exptime:0.0,filter:'',desext:'',desext2:'',$
                  dtacqnam:'',dsidname:'',expnum:'',writedate:'',bitpix:0L,shname:'',fname:'',newfile:'',deep:0,combfile:'',ncombine:0},nfiles)
for i=0,nfiles-1 do begin
  fil = file_basename(files[i])
  head = headfits(files[i],exten=0)
  ; Archive file
  if strmid(fil,0,3) eq 'c4d' or strmid(fil,0,2) eq 'tu' then begin
    head1 = headfits(files[i],exten=1)
    head2 = headfits(files[i],exten=2)
    base = strmid(file_basename(files[i]),0,17)

    ; Get header information
    orig_filename = sxpar(head,'DTACQNAM',/silent)
    fstr[i].dtacqnam = orig_filename
    ; /data_local/images/DTS/2013B-0440/DECam_00317720.fits.fz
    expnum = strmid(file_basename(orig_filename),6,8)  

    fstr[i].file = files[i]
    fstr[i].archivename = files[i]  ; the original name of the download file
    fstr[i].base = base
    fstr[i].desext = strtrim(sxpar(head1,'DES_EXT',/silent),2)
    fstr[i].desext2 = strtrim(sxpar(head2,'DES_EXT',/silent),2)
    fstr[i].dsidname = sxpar(head,'DTNSANAM',/silent)
    if strmid(fstr[i].dsidname,0,3) ne 'c4d' then begin  ; construct name
      dtutc = sxpar(head,'DTUTC',/silent)
      ; c4d_140106_052819, 2014-01-06T05:28:19
      fstr[i].dsidname = 'c4d_'+strmid(dtutc,2,2)+strmid(dtutc,5,2)+strmid(dtutc,8,2)+'_'+$
                         strmid(dtutc,11,2)+strmid(dtutc,14,2)+strmid(dtutc,17,2)+'_ori.fits'
    endif
    fstr[i].expnum = expnum
    fstr[i].writedate = sxpar(head1,'date',/silent)  
    if strmid(fstr[i].file,1,2,/reverse) eq 'fz' then fstr[i].bitpix=sxpar(head1,'zbitpix',/silent) else $
      fstr[i].bitpix = sxpar(head1,'bitpix',/silent)

  ; Chip file
  endif else begin
    fstr[i].file = files[i]
    dum = strsplit(files[i],'_',/extract)
    fstr[i].expnum = dum[0]

    ; Renamed chip file
    if strmid(fstr[i].expnum,0,1) eq 'F' and stregex(fstr[i].expnum,'F[0-9]+-*',/boolean) then begin
      fstr[i].newfile = fstr[i].expnum
    endif
  endelse
  fstr[i].object = strtrim(sxpar(head,'object',/silent),2)
  fstr[i].ra = sxpar(head,'ra',/silent)
  fstr[i].dec = sxpar(head,'dec',/silent)
  fstr[i].exptime = sxpar(head,'exptime',/silent)
  fstr[i].dateobs = sxpar(head,'DATE-OBS',/silent)
  filter = sxpar(head,'filter',count=nfilter,/silent)
  if nfilter eq 0 then filter=sxpar(head1,'filter',/silent)
  fstr[i].filter = strmid(strtrim(filter,2),0,1)

  fstr[i].ncombine = sxpar(head,'NCOMBINE',/silent)

  if not keyword_set(silent) then $
    print,strtrim(i+1,2),' ',file_basename(files[i]),'  ',fstr[i].expnum,'  ',fstr[i].object,'  ',fstr[i].dateobs,'  ',fstr[i].filter,'  ',$
        strtrim(fstr[i].exptime,2),'  ',fstr[i].ra,'  ',fstr[i].dec
endfor

end
