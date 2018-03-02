pro make_cluster_isoplots,redo=redo

; Make plots comparing isochrones to the clusters
dir = '/dl1/users/dnidever/smash/cp/red/photred/'
catdir = dir+'catalogs/final/v6/'
clusterdir = dir+'clusters/'
plotdir = clusterdir+'plots/'
if file_test(plotdir,/directory) eq 0 then file_mkdir,plotdir

setdisp
!p.font = 0

; Load the cluster catalog
cat = importascii(clusterdir+'smash_clusters.txt',/header)
ncat = n_elements(cat)
print,strtrim(ncat,2),' clusters'
smash = importascii('~/projects/SMASHRED/data/smash_fields_final.txt',/header)

; Load the Basti isochrones
print,'Loading the Basti isochrones'
isosolar = MRDFITS('~/isochrone/basti/new/solar/basti_solar_iso.fits',1)
isosolar_info = MRDFITS('~/isochrone/basti/new/solar/basti_solar_iso.fits',2)
isoaenh = MRDFITS('~/isochrone/basti/new/aenh/basti_aenh_iso.fits',1)
isoaenh_info = MRDFITS('~/isochrone/basti/new/aenh/basti_aenh_iso.fits',2)

; Loop over the clusters
for i=0,ncat-1 do begin
  print,strtrim(i+1,2),'/',strtrim(ncat,2),' ',cat[i].name
  outfile = plotdir+cat[i].name+'_isoplots.pdf'
  if file_test(outfile) eq 0 or keyword_set(redo) then begin
    ; Load the data
    catfile = clusterdir+cat[i].name+'_cat.fits'
    if file_test(catfile) eq 0 then begin
      print,catfile,' NOT FOUND'
      goto,BOMB
    endif
    obj = MRDFITS(catfile,1)
    ; Get the appropriate Basti isochrone
    if cat[i].afe ge 0.2 then begin
      isotype = 'aenh'
      ; Get closest metallicity
      uimh = uniq(isoaenh_info.mh,sort(isoaenh_info.mh))
      umh = isoaenh_info[uimh].mh
      mhind = first_el(minloc(abs(umh-cat[i].feh)))
      isomh = umh[mhind]
      gdmh = where(abs(isoaenh_info.mh-isomh) lt 0.01,ngdmh)
      ; Get closest age
      uiage = uniq(isoaenh_info[gdmh].age,sort(isoaenh_info[gdmh].age))
      uage = isoaenh_info[gdmh[uiage]].age
      ageind = first_el(minloc(abs(uage-cat[i].age)))
      isoage = uage[ageind]
      ; Final isochrone
      isoind = where(abs(isoaenh_info.mh-isomh) lt 0.01 and abs(isoaenh_info.age-isoage) lt 0.001,nisoind)
      isoinfo = isoaenh_info[isoind[0]]
      iso = isoaenh[isoinfo.index:isoinfo.index+isoinfo.np-1]
    endif else begin
      isotype = 'solar'
      ; Get closest metallicity
      uimh = uniq(isosolar_info.mh,sort(isosolar_info.mh))
      umh = isosolar_info[uimh].mh
      mhind = first_el(minloc(abs(umh-cat[i].feh)))
      isomh = umh[mhind]
      gdmh = where(abs(isosolar_info.mh-isomh) lt 0.01,ngdmh)
      ; Get closest age
      uiage = uniq(isosolar_info[gdmh].age,sort(isosolar_info[gdmh].age))
      uage = isosolar_info[gdmh[uiage]].age
      ageind = first_el(minloc(abs(uage-cat[i].age)))
      isoage = uage[ageind]
      ; Final isochrone
      isoind = where(abs(isosolar_info.mh-isomh) lt 0.01 and abs(isosolar_info.age-isoage) lt 0.001,nisoind)
      isoinfo = isosolar_info[isoind[0]]
      iso = isosolar[isoinfo.index:isoinfo.index+isoinfo.np-1]
    endelse
    print,'Catalog values:    [Fe/H]=',stringize(cat[i].feh,ndec=2),' age=',stringize(cat[i].age,ndec=2)
    print,'Isochrones values:  [M/H]=',stringize(isomh,ndec=2),' age=',stringize(isoage,ndec=2)

    undefine,psfiles

    ;DECam_u, R_u = 3.9631
    ;DECam_g, R_g = 3.1863
    ;DECam_r, R_r = 2.1401
    ;DECam_i, R_i = 1.5690
    ;DECam_z, R_z = 1.1957
    ;DECam_Y, R_Y = 1.0476
    ;For example:
    ;A_u = R_u * EBV_SFD98
    ext = {u:3.9631,g:3.1863,r:2.1401,i:1.5690,z:1.1957,y:1.0476}

    ; Plots
    ;=========
    ; g v. g-i
    ;---------
    file = plotdir+cat[i].name+'_ggi_iso'
    ps_open,file,/color,thick=4,/encap
    device,/inches,xsize=8.5,ysize=9.5
    loadcol,3
    black = fsc_color('black',0)
    gd = where(obj.g lt 50 and obj.i lt 50,ngd)
    hess,obj[gd].g-obj[gd].i,obj[gd].g,dx=0.02,dy=0.05,xr=[-1,3.5],yr=[25,12],/log,xtit='g-i',ytit='g',tit=cat[i].name+' (g vs. g-i)',charsize=1.5
    loadct,39,/silent
    isocol = iso.g+cat[i].ebv*ext.g - (iso.i+cat[i].ebv*ext.i)
    isomag = iso.g+cat[i].ebv*ext.g + cat[i].distmod
    oplot,isocol,isomag,co=250
    legend_old,['[M/H]='+stringize(isomh,ndec=2)+' age='+stringize(isoage,ndec=2)+' Gyr  distmod='+stringize(cat[i].distmod,ndec=2)+$
                '  E(B-V)='+stringize(cat[i].ebv,ndec=2)],textcolor=[250],charsize=1.4,/bottom,/left,box=0
    ps_close
    ps2png,file+'.eps',/eps
    push,psfiles,file

    ; r v. g-i
    ;---------
    file = plotdir+cat[i].name+'_rgi_iso'
    ps_open,file,/color,thick=4,/encap
    device,/inches,xsize=8.5,ysize=9.5
    loadcol,3
    black = fsc_color('black',0)
    gd = where(obj.g lt 50 and obj.r lt 50,ngd)
    hess,obj[gd].g-obj[gd].i,obj[gd].r,dx=0.02,dy=0.05,xr=[-1,3.5],yr=[25,11],/log,xtit='g-i',ytit='r',tit=cat[i].name+' (r vs. g-i)',charsize=1.5
    loadct,39,/silent
    isocol = iso.g+cat[i].ebv*ext.g - (iso.i+cat[i].ebv*ext.i)
    isomag = iso.r+cat[i].ebv*ext.r + cat[i].distmod
    oplot,isocol,isomag,co=250
    legend_old,['[M/H]='+stringize(isomh,ndec=2)+' age='+stringize(isoage,ndec=2)+' Gyr  distmod='+stringize(cat[i].distmod,ndec=2)+$
                '  E(B-V)='+stringize(cat[i].ebv,ndec=2)],textcolor=[250],charsize=1.4,/bottom,/left,box=0
    ps_close
    ps2png,file+'.eps',/eps
    push,psfiles,file

    ; z v. r-i
    ;---------
    file = plotdir+cat[i].name+'_zri_iso'
    ps_open,file,/color,thick=4,/encap
    device,/inches,xsize=8.5,ysize=9.5
    loadcol,3
    black = fsc_color('black',0)
    gd = where(obj.r lt 50 and obj.z lt 50,ngd)
    hess,obj[gd].r-obj[gd].i,obj[gd].z,dx=0.02,dy=0.05,xr=[-1,1.5],yr=[25,11],/log,xtit='r-i',ytit='z',tit=cat[i].name+' (z vs. r-i)',charsize=1.5
    loadct,39,/silent
    isocol = iso.r+cat[i].ebv*ext.r - (iso.i+cat[i].ebv*ext.i)
    isomag = iso.z+cat[i].ebv*ext.z + cat[i].distmod
    oplot,isocol,isomag,co=250
    legend_old,['[M/H]='+stringize(isomh,ndec=2)+' age='+stringize(isoage,ndec=2)+' Gyr  distmod='+stringize(cat[i].distmod,ndec=2)+$
                '  E(B-V)='+stringize(cat[i].ebv,ndec=2)],textcolor=[250],charsize=1.4,/bottom,/left,box=0
    ps_close
    ps2png,file+'.eps',/eps
    push,psfiles,file

    ; g-z v. g-i
    ;-----------
    file = plotdir+cat[i].name+'_gzgi_iso'
    ps_open,file,/color,thick=4,/encap
    device,/inches,xsize=9.5,ysize=9.5
    loadcol,3
    black = fsc_color('black',0)
    gd = where(obj.g lt 50 and obj.i lt 50 and obj.z lt 50,ngd)
    hess,obj[gd].g-obj[gd].i,obj[gd].g-obj[gd].z,dx=0.02,dy=0.05,xr=[-1,3.5],yr=[-1,3.5],/log,xtit='g-i',ytit='g-z',tit=cat[i].name+' (g-z vs. g-i)',charsize=1.5
    loadct,39,/silent
    isocol1 = iso.g+cat[i].ebv*ext.g - (iso.i+cat[i].ebv*ext.i)
    isocol2 = iso.g+cat[i].ebv*ext.g - (iso.z+cat[i].ebv*ext.z)
    isomag = iso.g+cat[i].ebv*ext.g + cat[i].distmod
    gdiso = where(isomag ge min(obj[gd].g),ngdiso)
    oplot,isocol1[gdiso],isocol2[gdiso],co=250
    legend_old,['[M/H]='+stringize(isomh,ndec=2)+' age='+stringize(isoage,ndec=2)+' Gyr  distmod='+stringize(cat[i].distmod,ndec=2)+$
                '  E(B-V)='+stringize(cat[i].ebv,ndec=2)],textcolor=[250],charsize=1.4,/bottom,/left,box=0
    ps_close
    ps2png,file+'.eps',/eps
    push,psfiles,file

    ; u-r v. g-i
    ;------------
    file = plotdir+cat[i].name+'_urgi_iso'
    ps_open,file,/color,thick=4,/encap
    device,/inches,xsize=9.5,ysize=9.5
    loadcol,3
    black = fsc_color('black',0)
    gd = where(obj.g lt 50 and obj.i lt 50 and obj.u lt 50 and obj.r lt 50,ngd)
    hess,obj[gd].g-obj[gd].i,obj[gd].u-obj[gd].r,dx=0.02,dy=0.05,xr=[-1,3.0],yr=[0,5],/log,xtit='g-i',ytit='u-r',tit=cat[i].name+' (u-r vs. g-i)',charsize=1.5
    loadct,39,/silent
    isocol1 = iso.g+cat[i].ebv*ext.g - (iso.i+cat[i].ebv*ext.i)
    isocol2 = iso.u+cat[i].ebv*ext.u - (iso.r+cat[i].ebv*ext.r)
    isomag = iso.g+cat[i].ebv*ext.g + cat[i].distmod
    gdiso = where(isomag ge min(obj[gd].g),ngdiso)
    oplot,isocol1[gdiso],isocol2[gdiso],co=250
    legend_old,['[M/H]='+stringize(isomh,ndec=2)+' age='+stringize(isoage,ndec=2)+' Gyr  (m-M)='+stringize(cat[i].distmod,ndec=2)+$
                '  E(B-V)='+stringize(cat[i].ebv,ndec=2)],textcolor=[250],charsize=1.4,/bottom,/left,box=0
    ps_close
    ps2png,file+'.eps',/eps
    push,psfiles,file

    ; Convert all to PDF and combine
    for j=0,n_elements(psfiles)-1 do spawn,'epstopdf '+psfiles[j]+'.eps'
    cmd = 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+outfile+' '
    cmd += strjoin(psfiles+'.pdf',' ')
    spawn,cmd

    ;stop
  endif else print,outfile,' EXISTS and /redo NOT set'
  BOMB:
endfor


stop

end
