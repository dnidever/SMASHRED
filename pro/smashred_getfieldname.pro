;+
;
; SMASHRED_GETFIELDNAME
;
; Get the closest matching SMASH field name given RA/DEC coordinates.
; The angular distance is also returned.
;
; INPUTS:
;  str        The input structure.  Must have RA and DEC fields.
;  /verbose   Print results to the screen.
;
; OUTPUTS:
;  info     The structure with information on the field, coordinates
;             and distance (in deg) from input coordinates.
;
; USAGE:
;  IDL>smashred_getfieldname,str,info
;
; By D.Nidever Dec. 2016
;-

pro smashred_getfieldname,str,info,verbose=verbose

undefine,info

nstr = n_elements(str)
if nstr eq 0 then begin
  print,'Syntax - smashred_getfieldname,str,info,verbose=verbose'
  return
endif

; Check for RA and DEC fields
if tag_exist(str,'RA') eq 0 or tag_exist(str,'DEC') eq 0 then begin
  print,'STR must have RA and DEC fields'
  return
endif

; Load the SMASH fields file
rootdir = SMASHRED_ROOTDIR()
smash = importascii(rootdir+'cp/red/photred/catalogs/pro/smash_fields_final.txt',/header,/silent)
;smash = importascii('/data/smash/cp/red/photred/catalogs/pro/smash_fields_final.txt',/header,/silent)
add_tag,smash,'field','',smash
smash.field = 'Field'+strtrim(smash.num,2)

; Add pilot fields
schema = smash[0]
struct_assign,{dum:''},schema
add = replicate(schema,3)
; FieldB RA=80.5746, DEC=-55.4250
add[0].field = 'FieldB'
add[0].radeg = 80.5746
add[0].dedeg = -55.4250
; FieldE RA=33.337246, DEC=-45.265196
add[1].field = 'FieldE'
add[1].radeg = 33.337246
add[1].dedeg = -45.265196
; FieldF RA=1.1352460, DEC=0.39827800
add[2].field = 'FieldF'
add[2].radeg = 1.1352460
add[2].dedeg = 0.39827800
smash = [smash,add]

; Create info structure
info = replicate({input_ra:0.0d0,input_dec:0.0d0,field:'',ra:0.0d0,dec:0.0d0,dist:0.0},nstr)
info.input_ra = str.ra
info.input_dec = str.dec
info.dist = -1

; Load the summary file information
if keyword_set(verbose) then print,'   NUM   IN_RA     IN_DEC    SMASH       RA       DEC       DIST  '
for i=0L,nstr-1 do begin

  ; Get the "standard" SMASH field number
  dist = sphdist(str[i].ra,str[i].dec,smash.radeg,smash.dedeg,/deg)
  bestind = first_el(minloc(dist))
  ;if dist[bestind] lt 1.2 then begin
  if dist[bestind] lt 2.2 then begin
    ;smashfield = 'Field'+strtrim(smash[bestind].num,2)
    info[i].field = smash[bestind].field
    info[i].ra = smash[bestind].radeg
    info[i].dec = smash[bestind].dedeg
    info[i].dist = min(dist)
    ;print,'This is SMASH field: ',smashfield
  endif else begin
    ;print,'NO SMASH field MATCH. Using observed field name: ',ofield
    ;smashfield = ofield
  endelse

  if keyword_set(verbose) then $
    print,i+1,info[i].input_ra,info[i].input_dec,info[i].field,info[i].ra,info[i].dec,info[i].dist,format='(I5,2F10.4,A10,2F10.4,F10.3)'
endfor

end
