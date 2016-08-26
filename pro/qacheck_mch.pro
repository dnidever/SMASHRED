pro qacheck_mch

forward_function trans_coo_dev
; Compile MATCHSTARS.PRO
RESOLVE_ROUTINE,'matchstars',/compile_full_file

dirs = file_search('/data/smash/cp/red/photred/20??????',/test_directory,count=ndirs)

for i=0,ndirs-1 do begin
  print,strtrim(i+1,2),'  '+dirs[i]
  cd,dirs[i]
  readcol,'fields',shortname,fieldname,format='A,A',/silent
  shind = where(strmid(fieldname,1,/reverse) eq 'sh',nshind)
  for j=0,nshind-1 do begin
    shfield1 = shortname[shind[j]]
    ;print,'  ',strtrim(j+1,2),'  '+shfield1
    cd,dirs[i]+'/'+shfield1
    mchfiles = file_search(shfield1+'-*_??.mch',count=nmchfiles)
    undefine,allmchstr
    for k=0,nmchfiles-1 do begin
      loadmch,mchfiles[k],mchlist,mchtrans
      mchlistbase = reform( (strsplitter(file_basename(mchlist,'.als'),'_',/extract))[0,*] )
      base = file_basename(mchfiles[k],'.mch')
      ichip = long( (strsplit(base,'_',/extract))[1] )
      mchstr = replicate({field:'',chip:0L,mchfile:'',alsfile:'',alsbase:'',trans0:fltarr(6),trans:fltarr(6)},n_elements(mchlist))
      mchstr.field = shfield1
      mchstr.chip = ichip
      mchstr.mchfile = mchfiles[k]
      mchstr.alsfile = mchlist
      mchstr.alsbase = mchlistbase
      mchstr.trans0 = transpose( mchtrans )

      ; Make array of X/Y points in reference frame
      nx = 2048
      ny = 4096
      xx = ( findgen(100)/99*(nx-2)+1 )#replicate(1,100)
      xx = xx[*]
      yy = replicate(1,100)#( findgen(100)/99*(ny-2)+1 )
      yy = yy[*]

      ; Use ref FITS header to convert X/Y to RA/DEC
      reffile = dirs[i]+'/'+shfield1+'/'+file_basename(mchlist[0],'.als')+'.fits'
      refhead = headfits(reffile)
      head_xyad,refhead,xx,yy,ra,dec,/deg

      ; Loop through files in the MCH file
      for l=0,n_elements(mchlist)-2 do begin
        file = dirs[i]+'/'+shfield1+'/'+file_basename(mchlist[l+1],'.als')+'.fits'
        ; Load header of second file
        head = headfits(file)
        ; Use second FITS header to convert RA/DEC to X2/Y2 in second image
        head_adxy,head,ra,dec,xx2,yy2,/deg
        ; Only select stars that fall within the second image
        gdstars = where(xx2 ge 0 and xx2 le nx-1 and yy2 ge 0 and yy2 le ny-1,ngdstars)
        ; if more than ~100 points then randomly subsample
        if ngdstars gt 100 then begin
          rndind = (sort(randomu(seed,ngdstars)))[0:99]
          gdstars0 = gdstars
          gdstars = gdstars0[rndind]
          ngdstars = n_elements(gdstars)
        endif
        ; Use leftover points and X/Y and X2/Y2 to measure transformation
        xdiff = xx[gdstars] - xx2[gdstars]
        ydiff = yy[gdstars] - yy2[gdstars]
        ;xdiff = refals[gdref[ind1]].x-als[gdals[ind2]].x
        ;ydiff = refals[gdref[ind1]].y-als[gdals[ind2]].y
        xmed = median([xdiff],/even)
        ymed = median([ydiff],/even)
        trans = [xmed, ymed, 1.0, 0.0, 0.0, 1.0]
        ; Fit rotation if there are enough points
        ;fa = {x1:refals[gdref[ind1]].x,y1:refals[gdref[ind1]].y,x2:als[gdals[ind2]].x,y2:als[gdals[ind2]].y}
        fa = {x1:xx[gdstars],y1:yy[gdstars],x2:xx2[gdstars],y2:yy2[gdstars]}
        initpar = [xmed, ymed, 1.0, 0.0, 0.0, 1.0]
        fpar = MPFIT('trans_coo_dev',initpar,functargs=fa, perror=perror, niter=iter, status=status,$
                      bestnorm=chisq, dof=dof, autoderivative=1, /quiet)
        trans = fpar
        mchstr[l+1].trans = trans
      endfor
      mchstr[0].trans = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
      xydiff = mchstr.trans[0:1] - mchstr.trans0[0:1]       
      ;print,strtrim(k+1,2),' ',mchfiles[k],' ',max(abs(xydiff))
      push,allmchstr,mchstr
    endfor

    xydiff = allmchstr.trans[0:1] - allmchstr.trans0[0:1]       
    print,'  ',strtrim(j+1,2),'  '+shfield1,'  ',max(abs(xydiff))

    ;; Make sure the shifts for each exposure are consistent for all chips
    ;uiexp = uniq(allmchstr.alsbase,sort(allmchstr.alsbase))
    ;uexp = allmchstr[uiexp].alsbase
    ;nexp = n_elements(uexp)
    ;expstr = replicate({exp:'',nchips:0L,xmn:0.0,ymn:0.0,xrms:0.0,yrms:0.0,xsig:0.0,ysig:0.0,sig:0.0,nbad:0L},nexp)
    ;for k=0,nexp-1 do begin
    ;  match,uexp[k],allmchstr.alsbase,ind1,ind2,/sort,count=nmatch
    ;  xtrans = allmchstr[ind2].trans[0]
    ;  ytrans = allmchstr[ind2].trans[1]
    ;  xmed = median([xtrans])
    ;  ymed = median([ytrans])
    ;  dxtrans = xtrans-xmed
    ;  dytrans = ytrans-ymed
    ;  xrms = stddev(xtrans)
    ;  yrms = stddev(ytrans)
    ;  xsig = mad(xtrans)
    ;  ysig = mad(ytrans)
    ;  diff = sqrt( dxtrans^2 + dytrans^2 )
    ;  sig = mad(diff,/zero)
    ;  bad = where(diff/(1.0<sig>0.1) gt 10,nbad)
    ;  expstr[k].exp = uexp[k]
    ;  expstr[k].nchips = nmatch
    ;  expstr[k].xmn = xmed
    ;  expstr[k].ymn = ymed
    ;  expstr[k].xrms = xrms
    ;  expstr[k].yrms = yrms
    ;  expstr[k].xsig = xsig
    ;  expstr[k].ysig = ysig
    ;  expstr[k].sig = sig
    ;  expstr[k].nbad = nbad
    ;  print,k+1,'  '+uexp[k],xmed,ymed,xrms,yrms,xsig,ysig,sig,nbad
    ;endfor
    ;if total(expstr.nbad) gt 0 then stop,'bad chips'
    ;stop
  endfor ; field loop
  print,''
  ;stop
endfor  ; directory loop


stop

end
