pro qacheck_astids

; Check the problem duplicate IDs in the AST files

;gaiastr = mrdfits('/data/smash/cp/red/photred/gaia/qacheck_gaiawcs.fits',1)
;goto,findoutliers

reduxdir = '/data/smash/cp/red/photred/'
dirs = file_search(reduxdir+'20??????',/test_directory,count=ndirs)
undefine,gaiastr
for i=0,ndirs-1 do begin
  files = file_search(dirs[i]+'/F*/F*-*_01.ast',count=nfiles)
  print,strtrim(i+1,2),'/',strtrim(ndirs,2),' ',dirs[i],' ',strtrim(nfiles,2)
  str1 = replicate({dir:'',file:'',okay:0},nfiles)
  for j=0,nfiles-1 do begin
    ast = importascii(files[j],/header,/silent)
    ui = uniq(ast.id,sort(ast.id))
    nui = n_elements(ui)
    if nui ne n_elements(ast) then str1[j].okay=0 else str1[j].okay=1

    str1[j].dir = dirs[i]
    str1[j].file = repstr(files[j],'.gaiawcs.head','.fits')
    print,'  ',file_basename(files[j]),' ',str1[j].okay
  endfor
  push,str,str1
endfor

stop

end
