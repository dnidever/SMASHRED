pro smashred_load_catphot,info,fstr,chstr,useast=useast

; Load all of the photometry for a PHOTRED catalog

fstr = *info.fstr
fbaseold = file_basename(info.file,'_summary.fits')  ; the observed field name
fbase = field  ; use the CORRECT field name with "sh"
if info.sh eq 1 then fbase+='sh'
print,strtrim(c+1,2),'/',strtrim(ninfo,2),'  ',field
chstr = mrdfits(info.file,2,/silent)  ; load the chip structure
nchstr = n_elements(chstr)
night = info.night

print,'Loading data for ',field

; Load in the data
;   don't have the "reference" frame information, just search
;   for PHOT files
print,'---Loading the data for ',field,'---'

add_tag,chstr,'refexpnum','',chstr
add_tag,chstr,'vertices_ra',dblarr(4),chstr
add_tag,chstr,'vertices_dec',dblarr(4),chstr
add_tag,chstr,'ndata',-1L,chstr
add_tag,chstr,'data',ptr_new(),chstr
; Use PHOT files
if not keyword_set(useast) then begin
  ext = '.phot'
  photfiles = file_search(reduxdir+night+'/'+chstr[0].field+'/'+chstr[0].field+'-*_??'+ext,count=nphotfiles)
  print,strtrim(nphotfiles,2),' PHOT files for ',fbase,' in ',reduxdir+night+'/'+chstr[0].field+'/'
; Use AST files
endif else begin
  ext = '.ast'
  photfiles = file_search(reduxdir+night+'/'+chstr[0].field+'/'+chstr[0].field+'-*_??.'+ext,count=nphotfiles)
  print,strtrim(nphotfiles,2),' AST files for ',fbase,' in ',reduxdir+night+'/'+chstr[0].field+'/'
endelse

; Loop through individual chip AST/PHOT files
for i=0,nphotfiles-1 do begin
  
  print,strtrim(i+1,2),' ',photfiles[i]
  phot = importascii(photfiles[i],/header,/silent)
  tags = tag_names(phot)
  phbase = file_basename(photfiles[i],ext)
  ichip = long(first_el(strsplit(phbase,'_',/extract),/last))
  dum = first_el(strsplit(phbase,'-',/extract),/last)
  refexpnum = first_el(strsplit(dum,'_',/extract))

  chind = where(chstr.chip eq ichip,nchind)
  ; Loop through all of the exposures
  for j=0,nchind-1 do begin
    mind = where(tags eq strtrim(chstr[chind[j]].calib_magname,2),nmind)
    gd = where(phot.(mind[0]) lt 50,ngd)
    ; cmag/cerr are for calibrated photometry that will be added later
    temp = replicate({id:-1L,x:0.0,y:0.0,mag:0.0,err:0.0,cmag:-1.0,cerr:-1.0,chi:0.0,sharp:0.0,flag:-1,prob:-1.0,ra:0.0d0,dec:0.0d0},ngd)
    struct_assign,phot[gd],temp,/nozero
    temp.mag = phot[gd].(mind[0])
    temp.err = phot[gd].(mind[0]+1)   ; assume error is the next column

    chstr[chind[j]].refexpnum = refexpnum
    chstr[chind[j]].ndata = ngd
    chstr[chind[j]].data = ptr_new(temp)  ; save the data

    ; Get astrometric vertices from header
    fitsfile = reduxdir+night+'/'+chstr[0].field+'/'+strtrim(chstr[chind[j]].base,2)+'.fits'
    head = headfits(fitsfile)
    nx = sxpar(head,'NAXIS1')
    ny = sxpar(head,'NAXIS2')
    head_xyad,head,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/degree
    chstr[chind[j]].vertices_ra = vra
    chstr[chind[j]].vertices_dec = vdec

    print,j+1,chstr[chind[j]].expnum,chstr[chind[j]].chip,ngd,format='(I4,A12,I5,I10)'
  endfor
endfor

;stop

end
