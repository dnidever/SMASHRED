pro lmc_completeness

  ;; measure the LMC RC completeness

  dir = '/dl1/users/dnidever/smash/cp/red/photred/catalogs/final/v6/ast/'
  ;dir = '/Users/nidever/smash/reduction/catalogs/final/v5/stars1/'
  ;dir = '/Users/nidever/smash/reduction/catalogs/final/v5/ast/'

  ;smash = importascii('/Users/nidever/smash/tiling/smash_fields_final.txt',/header)
  smash = importascii('~/projects/SMASHRED/data/smash_fields_final.txt',/header)
  add_tag,smash,'pa',0.0,smash
  add_tag,smash,'rad',0.0,smash
  add_tag,smash,'rcgrec',0.0,smash
  add_tag,smash,'gmag08_50',0.0,smash
  add_tag,smash,'gmag14_50',0.0,smash
  add_tag,smash,'rcirec',0.0,smash
  add_tag,smash,'imag08_50',0.0,smash
  add_tag,smash,'imag14_50',0.0,smash
  cel2lmc,smash.radeg,smash.dedeg,lmcpa,lmcrad
  smash.pa = lmcpa
  smash.rad = lmcrad
  gd = where(lmcrad lt 12,ngd)
  lmc = smash[gd]
  
  ;; Loop through the fields
  for i=0,ngd-1 do begin
     fname = 'Field'+strtrim(lmc[i].num,2)
     ;file = dir+fname+'_complete_stars.fits.gz'
     file = dir+fname+'_complete.fits.gz'
     if file_test(file) eq 0 then begin
       print,file,' NOT FOUND'
       goto,BOMB
     endif
     ast = mrdfits(file,1,/silent)
     ;ast = mrdfits(dir+fname+'_complete.fits.gz',1)
     ;gdrecover = where(ast.recovered eq 1,ngdrecover)
     gdrecover = where(ast.recovered eq 1 and ast.ndetg gt 0 and ast.ndeti gt 0,ngdrecover)
     dx = 0.2
     dy = 0.4
     xr = [-1,3.5]
     yr = [17.0,26.0]
     ; g vs. g-i
     hess,ast.inp_g-ast.inp_i,ast.inp_g,dum,gimall,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
     hess,ast[gdrecover].inp_g-ast[gdrecover].inp_i,ast[gdrecover].inp_g,dum,gimrec,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
     grecim = float(gimrec)/(gimall>1)
     ; i vs. g-i
     hess,ast.inp_g-ast.inp_i,ast.inp_i,dum,iimall,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
     hess,ast[gdrecover].inp_g-ast[gdrecover].inp_i,ast[gdrecover].inp_i,dum,iimrec,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
     irecim = float(iimrec)/(iimall>1)

     ;hess,str.inp_g-str.inp_i,str.inp_g,dum,imall,dx=0.02,dy=0.02,xr=[-1,3.5],yr=[24,16],xarr=xarr,yarr=yarr,/noplot
     ;recind = where(str.recovered eq 1,nrecind)
     ;rec = str[recind]
     ;hess,rec.inp_g-rec.inp_i,rec.inp_g,dum,imrec,dx=0.02,dy=0.02,xr=[-1,3.5],yr=[24,16],/noplot
     ;displayc,recim,xarr,yarr,/yflip

     ;rec = total(recim[9:10,*],1)/2.
     ;plot,yarr,rec,ps=-1

     ;The paper says io=18.47 and (g-i)o=0.82 so that gives go=19.29

     ; g-band completeness at the RC
     ;------------------------------
     ;maxrec = median(recim[*,0:9])  ;max(recim)
     maxgrec = max(grecim)
     grec = reform(grecim[9,6]) / maxgrec  ; g-i=0.8, g=19.4
     ; Get 50% completeness at g-i=0.8
     grecall = reform(grecim[9,*])
     yarr2 = scale_vector(findgen(1000),min(yarr),max(yarr))
     interp,yarr,grecall,yarr2,grecall2
     gind50 = first_el(where(grecall2/maxgrec lt 0.5))
     gmag08_50 = yarr2[gind50]
     ; 50% completeness at g-i=1.4
     grecall = reform(grecim[12,*])
     yarr2 = scale_vector(findgen(1000),min(yarr),max(yarr))
     interp,yarr,grecall,yarr2,grecall2
     gind50 = first_el(where(grecall2/maxgrec lt 0.5))
     if gind50 eq -1 then gind50=n_elements(grecall2)-1
     gmag14_50 = yarr2[gind50]
     ; i-band completeness at the RC
     ;------------------------------
     ;maxrec = median(recim[*,0:9])  ;max(recim)
     maxirec = max(irecim)
     irec = reform(irecim[9,4]) / maxirec  ; g-i=0.8, i=18.6
     ; Get 50% completeness at g-i=0.8
     irecall = reform(irecim[9,*])
     yarr2 = scale_vector(findgen(1000),min(yarr),max(yarr))
     interp,yarr,irecall,yarr2,irecall2
     iind50 = first_el(where(irecall2/maxirec lt 0.5))
     imag08_50 = yarr2[iind50]
     ; 50% completeness at g-i=1.4
     irecall = reform(irecim[12,*])
     yarr2 = scale_vector(findgen(1000),min(yarr),max(yarr))
     interp,yarr,irecall,yarr2,irecall2
     iind50 = first_el(where(irecall2/maxirec lt 0.5))
     if iind50 eq -1 then iind50=n_elements(irecall2)-1
     imag14_50 = yarr2[iind50]
     print,i+1,fname,grec,gmag08_50-19.3,gmag14_50-19.3,irec,imag08_50-18.5,imag14_50-18.5,format='(I5,A10,6F10.3)'
     ;print,i+1,fname,rec,gmag08_50,gmag14_50,format='(I5,A10,3F10.3)'
     lmc[i].rcgrec = grec
     lmc[i].gmag08_50 = gmag08_50
     lmc[i].gmag14_50 = gmag14_50
     lmc[i].rcirec = irec
     lmc[i].imag08_50 = imag08_50
     lmc[i].imag14_50 = imag14_50
     BOMB:
;if lmc[i].num gt 183 then stop
     ;stop
  endfor

  x = (90+lmc.dedeg)*cos(lmc.radeg/!radeg)
  y = (90+lmc.dedeg)*sin(lmc.radeg/!radeg)
  gd = where(lmc.rcgrec gt 0,ngd)
  plotc,x[gd],y[gd],lmc[gd].rcgrec,ps=1
  cel2lmc,lmc.radeg,lmc.dedeg,lmcpa,lmcrad  

  stop

  end
