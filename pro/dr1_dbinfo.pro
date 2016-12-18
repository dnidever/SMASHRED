pro dr1_dbinfo

; Make an ASCII file with information on the DR1 database tables.

dir = '/data/smash/cp/red/photred/catalogs/final/v4/db/'

tables = ['field','exposure','chip','source','object','xmatch']
ntables = n_elements(tables)

undefine,lines,allinfo

push,lines,'Description of SMASH DR1 Database tables'
push,lines,''

; Loop through the tables
for i=0,ntables-1 do begin
  ;str = mrdfits(dir+'Field100_'+tables[i]+'.fits',1)
  head = headfits(dir+'Field100_'+tables[i]+'.fits',exten=1)
  ncols = sxpar(head,'tfields')
  info = replicate({name:'',form:'',type:'',unit:'',ucd:'',description:''},ncols)
  push,lines,strupcase(tables[i])+' TABLE'
  push,lines,'----------------------------------------------------------------------------------------------------------------------------------------------'
  push,lines,'NUM  NAME             TYPE  UNITS          UCD                           DESCRIPTION'
  push,lines,'=============================================================================================================================================='
  for j=0,ncols-1 do begin
    info[j].name = strtrim(sxpar(head,'TTYPE'+strtrim(j+1,2)),2)
    info[j].form = strtrim(sxpar(head,'TFORM'+strtrim(j+1,2)),2)
    info[j].unit = strtrim(sxpar(head,'TUNIT'+strtrim(j+1,2)),2)
    info[j].ucd = strtrim(sxpar(head,'TUCD'+strtrim(j+1,2)),2)
    info[j].description = strtrim(sxpar(head,'TCOMM'+strtrim(j+1,2)),2)
    type = info[j].form
    if strlen(type) gt 1 then type='A' 
    line = string(j+1,format='(I-5)')+string(info[j].name,format='(A-19)')+string(type,format='(A-4)')+$
           string(info[j].unit,format='(A-15)')+string(info[j].ucd,format='(A-30)')+$
           string(info[j].description,format='(A-70)')
    push,lines,line
  endfor
  push,lines,'----------------------------------------------------------------------------------------------------------------------------------------------'
  push,lines,''
  push,allinfo,info

  ;stop

endfor

printline,lines

writeline,'dr1_dbinfo.txt',lines

stop

end
