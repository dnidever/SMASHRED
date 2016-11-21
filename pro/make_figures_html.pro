pro make_figures_html

dir = '/data/smash/cp/red/photred/'
figdir = '/d0/choi/Figures/DEcam_stamp/'
outdir = '/home/dnidever/public_html/smash/qahtml/'
;outdir = '/data/smash/cp/red/photred/qahtml/'

nightdirs = file_search(dir+'20??????',/test_directory,count=nnights)

nightstr = replicate({name:'',dir:'',nexp:0L},nnights)
for i=0,nnights-1 do begin

  night = file_basename(nightdirs[i])
  sumfiles = file_search(nightdirs[i]+'/*_summary.fits',count=nsumfiles)
  undefine,fstr,fstr1
  for j=0,nsumfiles-1 do begin
    fstr1 = mrdfits(sumfiles[j],1,/silent)
    chstr1 = mrdfits(sumfiles[j],2,/silent)
    add_tag,fstr1,'field','',fstr1
    fstr1.field = chstr1[0].field
    push,fstr,fstr1
  endfor
  nfstr = n_elements(fstr)
  nightstr[i].name = night
  nightstr[i].dir = nightdirs[i]
  nightstr[i].nexp = nfstr
  print,night,' ',strtrim(nfstr,2),' exposures'

  ; Make the HTML page
  htmlfile = outdir+night+'_qapsf.html'

  ; Write the HTML file
  openw,unit,/get_lun,htmlfile

  printf,unit,'<HTML>'
  printf,unit,'<BODY>'

  ; Make an index
  ;if keyword_set(index) then begin
  ;  printf,unit,'<h1>Index</h1><br>'
  ;  for i=0,n-1 do begin
  ;    printf,unit,'<a href="#'+names[i]+'">'+names[i]+'</a>'
  ;  endfor
  ;  printf,unit,'<p>'
  ;endif

  printf,unit,'<h1>Night '+night+' - '+strtrim(nfstr,2)+' Exposures</h1>'
  printf,unit,'<hr>'
  printf,unit,'<h4>Clicking on exposure name will open a new tab with the three images for that exposure.<br>'
  printf,unit,'Clicking on an image will open a larger version of that image in a new tab.</h4>'

  printf,unit,'<table border=1>'
  printf,unit,'<tr><th>Num</th><th>Name</th><th>Image</th><th>PSF</th><th>Difference</th></tr>'

  for j=0,nfstr-1 do begin
    expnum = fstr[j].field+'-'+fstr[j].expnum
    printf,unit,'<tr>'
    printf,unit,'<td>'+strtrim(j+1,2)+'</td>'
    printf,unit,'<td><a href="'+expnum+'_qapsf.html" target="main">'+expnum+'</a></td>'
    printf,unit,'<td><a href="'+expnum+'_img.png" target="main"><img src="'+expnum+'_img.png" width=400></a></td>'
    printf,unit,'<td><a href="'+expnum+'_psf.png" target="main"><img src="'+expnum+'_psf.png" width=400></a></td>'
    printf,unit,'<td><a href="'+expnum+'_sub.png" target="main"><img src="'+expnum+'_sub.png" width=400></a></td>'
    printf,unit,'</tr>'
  endfor
  printf,unit,'</table>'
  printf,unit,'</BODY>'
  printf,unit,'</HTML>'
  close,unit
  free_lun,unit
  print,'HTML file ',htmlfile,' created'

  ; Make individual html pages per exposure
  for j=0,nfstr-1 do begin
    expnum = fstr[j].field+'-'+fstr[j].expnum
    htmlfile = outdir+expnum+'_qapsf.html'
    openw,unit,/get_lun,htmlfile
    printf,unit,'<HTML>'
    printf,unit,'<BODY>'
    printf,unit,'<h1>Exposure '+expnum+'</h1>'
    printf,unit,'<hr>'
    printf,unit,'<h4>Clicking on an image will open a larger version of that image in a new tab.</h4>'
    printf,unit,'<table border=1>'
    printf,unit,'<tr><th>Image</th><th>PSF</th><th>Difference</th></tr>'
    printf,unit,'<tr>'
    printf,unit,'<td><a href="'+expnum+'_img.png" target="main"><img src="'+expnum+'_img.png" width=400></a></td>'
    printf,unit,'<td><a href="'+expnum+'_psf.png" target="main"><img src="'+expnum+'_psf.png" width=400></a></td>'
    printf,unit,'<td><a href="'+expnum+'_sub.png" target="main"><img src="'+expnum+'_sub.png" width=400></a></td>'
    printf,unit,'</tr>'
    printf,unit,'</table>'
    printf,unit,'</BODY>'
    printf,unit,'</HTML>'
    close,unit
    free_lun,unit
    print,'HTML file ',htmlfile,' created'
  endfor

endfor

; Make the overall index page
htmlfile = outdir+'index.html'
openw,unit,/get_lun,htmlfile
printf,unit,'<HTML>'
printf,unit,'<BODY>'
printf,unit,'<h1>SMASH PSF Quality Assurance Figures</h1>'
printf,unit,'<p>'
printf,unit,'<h2>'+strtrim(nnights,2),' nights</h2>'
printf,unit,'<hr>'
printf,unit,'<h4>Clicking on a link will open up a new tab for that night.</h4>'
printf,unit,'<table border=1>'
printf,unit,'<tr><th>Num</th><th>Night</th><th>Nexposures</th></tr>'
for i=0,nnights-1 do begin
  printf,unit,'<tr>'
  printf,unit,'<td>'+strtrim(i+1,2)+'</td>'
  printf,unit,'<td><a href="'+nightstr[i].name+'_qapsf.html" target="main">'+nightstr[i].name+'</a></td>'
  printf,unit,'<td align=center>'+strtrim(nightstr[i].nexp,2)+'</td>'
  printf,unit,'</tr>'
endfor
printf,unit,'</table>'
printf,unit,'</BODY>'
printf,unit,'</HTML>'
close,unit
free_lun,unit
print,'HTML file ',htmlfile,' created'

stop

end
