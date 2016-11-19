pro qacheck_mch,night,field,verbose=verbose,allfields=allfields,thresh=thresh

forward_function trans_coo, trans_coo_dev
; Compile MATCHSTARS.PRO
RESOLVE_ROUTINE,'matchstars',/compile_full_file

if n_elements(thresh) eq 0 then thresh=5.0  ; pixel difference threshold

if n_elements(night) gt 0 then begin
  dirs = file_search('/data/smash/cp/red/photred/'+night,/test_directory,count=ndirs)
endif else begin
  dirs = file_search('/data/smash/cp/red/photred/20??????',/test_directory,count=ndirs)
endelse
print,strtrim(ndirs,2),' nights'

for i=0,ndirs-1 do begin
;for i=29,ndirs-1 do begin
  print,strtrim(i+1,2),'  '+dirs[i]
  cd,dirs[i]
  readcol,'fields',shortname,fieldname,format='A,A',/silent
  if not keyword_set(allfields) then begin
    shind = where(strmid(fieldname,1,/reverse) eq 'sh',nshind)
  endif else begin  ; running on all fields
    shind = lindgen(n_elements(fieldname))
    nshind = n_elements(fieldname)
  endelse

  ; Specific field input
  if n_elements(field) gt 0 then begin
    MATCH,shortname,field,shind,ind2,/sort,count=nshind
    ;shind = where(shortname eq field,nshind)
    if nshind eq 0 then begin
      print,field,' NOT FOUND'
      return
    endif
  endif

  ; Loop through the fields
  for j=0,nshind-1 do begin
    shfield1 = shortname[shind[j]]
    if keyword_set(verbose) then print,' ',strtrim(j+1,2),'  '+shfield1
    cd,dirs[i]+'/'+shfield1
    mchfiles = file_search(shfield1+'-*_??.mch',count=nmchfiles)
    if nmchfiles eq 0 then goto,bomb
    undefine,allmchstr
    for k=0,nmchfiles-1 do begin
      loadmch,mchfiles[k],mchlist,mchtrans
      mchlistbase = reform( (strsplitter(file_basename(mchlist,'.als'),'_',/extract))[0,*] )
      base = file_basename(mchfiles[k],'.mch')
      ichip = long( (strsplit(base,'_',/extract))[1] )
      mchstr = replicate({field:'',chip:0L,mchfile:'',alsfile:'',alsbase:'',trans0:fltarr(6),$
                          trans:fltarr(6),midxdiff:0.0,midydiff:0.0,meddiff:0.0},n_elements(mchlist))
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

        ; Compare X/Y transformations in the overlap region
        ;  using the random positions previously created
        out = trans_coo(xx2[gdstars],yy2[gdstars],mchstr[l+1].trans0)
        out0 = trans_coo(xx2[gdstars],yy2[gdstars],mchstr[l+1].trans)
        mchstr[l+1].meddiff = median(abs(out-out0))
        ; Compare X/Y transformations in the middle of the chip
        midout = trans_coo(nx/2,ny/2,mchstr[l+1].trans0)
        midout0 = trans_coo(nx/2,ny/2,mchstr[l+1].trans)
        mchstr[l+1].midxdiff = midout[0]-midout0[0]
        mchstr[l+1].midydiff = midout[1]-midout0[1]
if mchstr[l+1].meddiff gt thresh then print,mchstr[l+1].alsfile,' ',mchstr[l+1].meddiff
;if max(abs([mchstr[l+1].midxdiff,mchstr[l+1].midydiff])) gt 2 then print,mchstr[l+1].alsfile,' ',out-out0
;if max(abs([mchstr[l+1].midxdiff,mchstr[l+1].midydiff])) gt 2 then stop,'high offset'
;if file_basename(mchlist[l+1],'.als') eq '
      endfor ; files in MCH file
      mchstr[0].trans = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
      ;;xydiff = mchstr.trans[0:1] - mchstr.trans0[0:1]       
      ;;maxdiff = max(abs(xydiff))
      ;xdiff = mchstr.midxdiff
      ;ydiff = mchstr.midydiff
      ;maxdiff = max(abs([xdiff,ydiff]))
      maxdiff = max(mchstr.meddiff)
      if maxdiff gt thresh then cmt=' ***' else cmt=''
      if keyword_set(verbose) then print,'  ',strtrim(k+1,2),' ',mchfiles[k],' ',maxdiff,cmt
      push,allmchstr,mchstr
;if maxdiff gt 2 then stop,'large offset'
    endfor  ; MCH file loop

    ;;xydiff = allmchstr.trans[0:1] - allmchstr.trans0[0:1]       
    ;;maxdiff = max(abs(xydiff))
    ;xdiff = allmchstr.midxdiff
    ;ydiff = allmchstr.midydiff
    ;maxdiff = max(abs([xdiff,ydiff]))
    maxdiff = max(allmchstr.meddiff)
    print,'  ',strtrim(j+1,2),'  '+shfield1,'  ',maxdiff

    bomb:

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

;stop

end
