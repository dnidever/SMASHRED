pro read_votable,file,str

undefine,str

if n_elements(file) eq 0 then begin
  print,'Syntax - read_votable,file,str'
  return
endif

if file_test(file) eq 0 then begin
  print,file,' NOT FOUND'
  return
endif

readline,file,lines
lines0 = lines
lines = strtrim(lines,2)
bd = where(lines eq '',nbd)
if nbd gt 0 then remove,bd,lines

lotable = where(stregex(lines,'^<TABLE ',/boolean) eq 1,nlotable)
hitable = where(stregex(lines,'^</TABLE>',/boolean) eq 1,nhitable)
tlines = lines[lotable:hitable]

; Read in header lines
loheader = first_el(where(stregex(tlines,'^<FIELD',/boolean) eq 1))
hiheader = where(stregex(tlines,'^<DATA',/boolean) eq 1)-1
header = tlines[loheader:hiheader]
headind = where(stregex(header,'^<FIELD ',/boolean) eq 1,nheader)

namearr = strarr(nheader)
typearr = lonarr(nheader)
datatypearr = strarr(nheader)
for i=0,nheader-1 do begin

  if i ne nheader-1 then begin
    ihead = header[headind[i]:headind[i+1]-1]
  endif else begin
    ihead = header[headind[i]:*]
  endelse
  hd = strjoin(ihead,' ')
  arr = strsplit(hd,' ',/extract)
  datatype_ind = where(stregex(arr,'datatype=',/boolean),ndatatype_ind)
  name_ind = where(stregex(arr,'name="',/boolean),nname_ind)

  datatype = first_el(strsplit(arr[datatype_ind[0]],'=',/extract),/last)
  datatype = repstr(datatype,'"','')
  datatypearr[i] = datatype
  case datatype of
  'char': typearr[i]=7
  'float': typearr[i]=5 ;4
  'int': typearr[i]=2
  'long': typearr[i]=3
  'byte': typearr[i]=1
  'double': typearr[i]=5
  else: typearr[i]=-1
  endcase

  name = first_el(strsplit(arr[name_ind[0]],'=',/extract),/last)
  name = repstr(name,'"','')
  namearr[i] = name

  if typearr[i] eq 7 then zero='' else zero=fix(0.0,type=typearr[i])
  if i eq 0 then dum=create_struct(name,zero) else dum=create_struct(dum,name,zero)
  
endfor

; Read in data
lodata = where(stregex(tlines,'^<TABLEDATA>',/boolean) eq 1,nlodata)
hidata = where(stregex(tlines,'^</TABLEDATA>',/boolean) eq 1,nhidata)
datalines = tlines[lodata+1:hidata-1]

dataind = where(stregex(datalines,'^<TR>',/boolean) eq 1,ndataind)

str = replicate(dum,ndataind)

for i=0,ndataind-1 do begin

  if i ne ndataind-1 then begin
    idata = datalines[dataind[i]:dataind[i+1]-1]
  endif else begin
    idata = datalines[dataind[i]:*]
  endelse
  idata = idata[1:*]
  idata = idata[0:n_elements(idata)-2]

  idata = repstr(idata,'<TD>','')
  idata = repstr(idata,'</TD>','')

  for j=0,nheader-1 do str[i].(j)=idata[j]

endfor

;stop

end
