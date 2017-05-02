pro laf_thumbnails

; Make thumbnails for the faint, blue LAF sources

rootdir = smashred_rootdir()+'cp/red/photred/'
dir = rootdir+'catalogs/final/v5/'
field = 'Field170'
psdir = '/datalab/users/dnidever/smash/cp/red/photred/lafthumbs/'
;psdir = '/home/dnidever/papers/smash_laf/lafthumbs/'

; Load the catalogs
chstr = mrdfits(dir+field+'_combined_chips.fits.gz',1)
chstr.file = strtrim(chstr.file,2)
expstr = mrdfits(dir+field+'_combined_exposures.fits.gz',1)
add_tag,expstr,'expind',0L,expstr
expstr.expind = lindgen(n_elements(expstr))
allsrc = mrdfits(dir+field+'_combined_allsrc.fits.gz',1)
allsrc.fid = strtrim(allsrc.fid,2)
stars = mrdfits(dir+'stars1/'+field+'_allobj_stars.fits.gz',1)
stars.id = strtrim(stars.id,2)

; Deredden
;magext = [4.239, 3.303, 2.285, 1.698, 1.263]
magext = [5.155, 3.793, 2.751, 2.086, 1.479]  ; from SFD98
stars.u -= magext[0]*stars.ebv
stars.g -= magext[1]*stars.ebv
stars.r -= magext[2]*stars.ebv
stars.i -= magext[3]*stars.ebv
stars.z -= magext[4]*stars.ebv

; Get best deep exposure in each band
bands = ['u','g','r','i','z']
nbands = n_elements(bands)
for i=0,nbands-1 do begin
  ind = where(expstr.filter eq bands[i] and expstr.exptime gt 200 and expstr.photometric eq 1,nind)
  if nind eq 0 then ind = where(expstr.filter eq bands[i] and expstr.exptime gt 200,nind)
  bestind = first_el(minloc(expstr[ind].fwhm))
  push,bestexp,expstr[ind[bestind]]
endfor

; Pick the sources
nstars = n_elements(stars)
srcindex = lon64arr(nstars,5)
for i=0,nbands-1 do srcindex[*,i]=stars.srcfindx[bestexp[i].expind]
numdet = total(srcindex ge 0,2)
gd = where(stars.g-stars.i ge -0.3 and stars.g-stars.i le 0.3 and stars.g ge 22.5 and stars.g le 24.0 and numdet eq nbands,ngd) 
; require that it's detected in all the "best" exposures

;goto,html

; Loop through the sources
for i=0,ngd-1 do begin
  objid = stars[gd[i]].id
  print,strtrim(i+1,2),' ',objid

  for j=0,nbands-1 do begin
    srcind = stars[gd[i]].srcfindx[bestexp[j].expind]
    allsrc1 = allsrc[srcind]
    x = allsrc1.x-1
    y = allsrc1.y-1
    chipindx = allsrc1.chipindx    
    chstr1 = chstr[chipindx]
    chipfile = rootdir+'deep/'+field+'/'+strtrim(chstr1.field,2)+'/'+strtrim(chstr1.file,2)

    psfile = psdir+objid+'_'+bands[j]
    if file_test(psfile+'.png') eq 0 or keyword_set(redo) then begin
      dx = 20
      subim = fltarr(2*dx+1,2*dx+1)
      fits_read,chipfile,im,head
      sz = size(im)
      nx = sz[1] & ny=sz[2]
      x0 = round(x-dx)
      y0 = round(y-dx)
      xlo = round(x-dx) > 0
      xhi = round(x+dx) < (nx-1)
      ylo = round(y-dx) > 0
      yhi = round(y+dx) < (ny-1)
      subim[xlo-x0:xhi-x0,ylo-y0:yhi-y0] = im[xlo:xhi,ylo:yhi]
      xx = lindgen(2*dx+1)+x0
      yy = lindgen(2*dx+1)+y0

      !p.font = 0
      ps_open,psfile,/color,thick=4,/encap
      device,/inches,xsize=9.5,ysize=9.0
      loadcol,3
      black = fsc_color('black',0)
      displayc,subim,xx,yy,/z,xtit='X',ytit='Y',tit=objid+' '+bands[j]+'='+stringize(allsrc1.cmag,ndec=2)+' ('+chstr1.file+')',charsize=1.1
      loadct,39,/silent
      oplot,[x],[y],ps=1,sym=3,co=250
      ps_close
      ps2png,psfile+'.eps',/eps
    endif

    ; Make PSF-subtracted thumbnail as well
    psfile2 = psdir+objid+'_'+bands[j]+'_psfsub'
    if file_test(psfile2+'.png') eq 0 or keyword_set(redo) then begin
      dx = 20
      subim = fltarr(2*dx+1,2*dx+1)
      subchipfile = repstr(chipfile,'.fits','s.fits.fz')
      fits_read,subchipfile,sim,shead
      sz = size(sim)
      nx = sz[1] & ny=sz[2]
      x0 = round(x-dx)
      y0 = round(y-dx)
      xlo = round(x-dx) > 0
      xhi = round(x+dx) < (nx-1)
      ylo = round(y-dx) > 0
      yhi = round(y+dx) < (ny-1)
      subim[xlo-x0:xhi-x0,ylo-y0:yhi-y0] = sim[xlo:xhi,ylo:yhi]
      xx = lindgen(2*dx+1)+x0
      yy = lindgen(2*dx+1)+y0

      !p.font = 0
      ps_open,psfile2,/color,thick=4,/encap
      device,/inches,xsize=9.5,ysize=9.0
      loadcol,3
      black = fsc_color('black',0)
      displayc,subim,xx,yy,/z,xtit='X',ytit='Y',tit=objid+' '+bands[j]+'='+stringize(allsrc1.cmag,ndec=2)+' ('+file_basename(subchipfile)+')',charsize=1.1
      loadct,39,/silent
      oplot,[x],[y],ps=1,sym=3,co=250
      ps_close
      ps2png,psfile2+'.eps',/eps
    endif

    BOMB:
  endfor  ; band loop
endfor  ; object loop

; Make the HTML page
HTML:
htmlfile = psdir+'laf_thumbnails.html'

; Write the HTML file
openw,unit,/get_lun,htmlfile

printf,unit,'<HTML>'
printf,unit,'<BODY>'

printf,unit,'<h1>Faint, blue LAF objects in field '+field+'</h1>'
printf,unit,'<hr>'
printf,unit,'<table border=1>'
printf,unit,'<tr><th>Num</th><th>Object</th><th>u</th><th>g</th><th>r</th><th>i</th><th>z</th></tr>'

for i=0,ngd-1 do begin
  objid = stars[gd[i]].id
  printf,unit,'<tr>'
  printf,unit,'<td>'+strtrim(i+1,2)+'</td>'
  printf,unit,'<td>'+objid+'</td>'
  printf,unit,'<td><a href="'+objid+'_u.png" target="main"><img src="'+objid+'_u.png" height=300></a></td>'
  printf,unit,'<td><a href="'+objid+'_g.png" target="main"><img src="'+objid+'_g.png" height=300></a></td>'
  printf,unit,'<td><a href="'+objid+'_r.png" target="main"><img src="'+objid+'_r.png" height=300></a></td>'
  printf,unit,'<td><a href="'+objid+'_i.png" target="main"><img src="'+objid+'_i.png" height=300></a></td>'
  printf,unit,'<td><a href="'+objid+'_z.png" target="main"><img src="'+objid+'_z.png" height=300></a></td>'
  printf,unit,'</tr>'
  printf,unit,'<tr>'
  printf,unit,'<td></td><td></td>'
  printf,unit,'<td><a href="'+objid+'_u_psfsub.png" target="main"><img src="'+objid+'_u_psfsub.png" height=300></a></td>'
  printf,unit,'<td><a href="'+objid+'_g_psfsub.png" target="main"><img src="'+objid+'_g_psfsub.png" height=300></a></td>'
  printf,unit,'<td><a href="'+objid+'_r_psfsub.png" target="main"><img src="'+objid+'_r_psfsub.png" height=300></a></td>'
  printf,unit,'<td><a href="'+objid+'_i_psfsub.png" target="main"><img src="'+objid+'_i_psfsub.png" height=300></a></td>'
  printf,unit,'<td><a href="'+objid+'_z_psfsub.png" target="main"><img src="'+objid+'_z_psfsub.png" height=300></a></td>'
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
