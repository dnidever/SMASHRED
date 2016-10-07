;+
;
; SMASHRED_MEASURE_MAGOFFSET_MAINBODY
;
; This program measures photometric offsets between pairs
; of main-body fields.
;
; INPUTS:
;  /silent      Don't print anything to the screen.
;  /verbose     Print lots of information to the screen.
;
; OUTPUTS:
;  =error      The error message if one occurred.
;
; USAGE:
;  IDL>smashred_measure_magoffset_mainbody
;
; By D.Nidever  Oct 2016
;-

pro smashred_measure_magoffset_mainbody,version,smc=smc,dir=dir,silent=silent,verbose=verbose,error=error

undefine,overlapstr,error

nchips = n_elements(chstr)
nallsrc = n_elements(allsrc)

if n_elements(dir) eq 0 then begin
  if n_elements(version) eq 0 then begin
    print,'Please enter the catalog version, e.g. "v2"'
    return
  endif
  dir = '/data/smash/cp/red/photred/catalogs/final/'+version+'/'
endif

; Load the field overlap file
fieldoverlapstr = MRDFITS(dir+'smash_fieldoverlaps.fits',1,/silent)
fieldoverlapstr.field = strtrim(fieldoverlapstr.field,2)  ; remove any trailing spaces

; Pick out the main body fields
;  LMC
if not keyword_set(smc) then begin
  galtype = 'lmc'
  cel2lmc,fieldoverlapstr.cenra,fieldoverlapstr.cendec,lmcpa,lmcrad
  gdfields = where(lmcrad lt 14 and fieldoverlapstr.noverlap gt 0,ngdfields)
;  SMC
endif else begin
  galtype = 'smc'
  cel2smc,fieldoverlapstr.cenra,fieldoverlapstr.cendec,smcpa,smcrad
  gdfields = where(smcrad lt 6 and fieldoverlapstr.noverlap gt 0,ngdfields)
endelse

filters = ['u','g','r','i','z']
nfilters = n_elements(filters)

;overlapdata = mrdfits(dir+'lmcoverlap_magoffsets.fits',1)
;noverlapdata = n_elements(overlapdata)
;;goto,skiptohere
;goto,doplots

; Initialize overlapstr structure
overlapdata_schema = {field1:'',index1:-1L,field2:'',index2:-1L,nmatch:-1L,filter:'',ngood:lonarr(5)-1,magoff:fltarr(5)+99.99,magofferr:fltarr(5)+9.99}
overlapdata = replicate(overlapdata_schema,total(fieldoverlapstr[gdfields].noverlap))
noverlapdata = n_elements(overlapdata)
count = 0

; Field outer loop
For i=0,ngdfields-1 do begin
  field1 = strtrim(fieldoverlapstr[gdfields[i]].field,2)
  noverlap = fieldoverlapstr[gdfields[i]].noverlap
  overlapindex = fieldoverlapstr[gdfields[i]].overlapindex
  calib1 = fieldoverlapstr[gdfields[i]].calib
  print,strtrim(i+1,2),'/',strtrim(ngdfields,2),' ',field1,' ',strtrim(noverlap,2)
  if total(calib1) eq 0 then begin
    print,'Not calibrated.  Skipping'
    goto,BOMB1
  endif
  brightfile1 = dir+field1+'_combined_allobj_bright.fits'
  str1 = MRDFITS(brightfile1,1,/silent)
  tags1 = tag_names(str1)

  ; Field inner loop
  For j=0,noverlap-1 do begin
    field2 = strtrim(fieldoverlapstr[overlapindex[j]].field,2)
    calib2 = fieldoverlapstr[overlapindex[j]].calib
    print,'  ',strtrim(j+1,2),' ',field2
    if total(calib2) eq 0 then begin
      print,'Not calibrated.  Skipping'
      goto,BOMB2
    endif
    brightfile2 = dir+field2+'_combined_allobj_bright.fits'
    str2 = MRDFITS(brightfile2,1,/silent)
    tags2 = tag_names(str2)

    ; Match then up
    SRCMATCH,str1.ra,str1.dec,str2.ra,str2.dec,0.5,ind1,ind2,/sph,count=nmatch
    print,'    ',strtrim(nmatch,2),' matches'
    str1b = str1[ind1]
    str2b = str2[ind2]

    ; Fill in some info
    overlapdata[count].field1 = field1
    overlapdata[count].index1 = gdfields[i]
    overlapdata[count].field2 = field2
    overlapdata[count].index2 = overlapindex[j]
    overlapdata[count].nmatch = nmatch

    ; Filter loop
    For f=0,nfilters-1 do begin

      if calib1[f] eq 1 and calib2[f] eq 1 then begin

        magind1 = where(tags1 eq strupcase(filters[f]),nmagind1)
        errind1 = where(tags1 eq strupcase(filters[f])+'ERR',nerrind1)
        magind2 = where(tags2 eq strupcase(filters[f]),nmagind2)
        errind2 = where(tags2 eq strupcase(filters[f])+'ERR',nerrind2)

        ; Measure the mag offset
        mag1 = str1b.(magind1)
        err1 = str1b.(errind1)
        mag2 = str2b.(magind2)
        err2 = str2b.(errind2)

        magdiff = mag1 - mag2
        magerr = sqrt( err1^2.0 + err2^2.0 )
        gdmag = where(mag1 lt 50. and mag2 lt 50. and magerr lt 0.07,ngdmag)
        if ngdmag lt 10 then $   ; not enough points, lower error threshold
          gdmag = where(mag1 lt 50. and mag2 lt 50. and magerr lt 0.1,ngdmag)
        if ngdmag lt 10 then $   ; not enough points, lower error threshold
          gdmag = where(mag1 lt 50. and mag2 lt 50. and magerr lt 0.2,ngdmag)
        if ngdmag lt 10 then $   ; not enough points, lower error threshold
          gdmag = where(mag1 lt 50. and mag2 lt 50. and magerr lt 0.5,ngdmag)

        ; Some sources with decent photometry
        if ngdmag gt 0 then begin
          ROBUST_MEAN,magdiff[gdmag],magoffset,magoffsig,sig=magerr[gdmag],numrej=numrej,/usemad
          magofferr = magoffsig/sqrt(ngdmag-numrej)  ; stdev of mean
        endif

        ; Plotting
        if keyword_set(pl) and ngdmag gt 0 then begin
          plotc,mag1[gdmag],magdiff[gdmag],magerr[gdmag],ps=3,xr=[14,23],yr=[-1,1],xs=1,ys=1,$
                xtit=filters[f],ytit='Residuals',tit=field1+'-'+fields2+' ('+filters[f]+')'
          al_legend,[stringize(magoffset,ndec=4)+'+/-'+stringize(magofferr,ndec=4)+' mag'],/top,/left,charsize=1.5
        endif

        ; Save the information in the OVERLAPSTR structure
        if nmatch gt 1 and ngdmag gt 1 then begin
          ;magofferr = magoffsig/sqrt(nmatch)  ; error in the mean
          if magoffsig lt 1e-8 then magofferr=10.0  ; if it's zero that normally means it's pretty bad
          ;magoffsig = magoffsig > 1e-6  ; let lower limit, can sometimes be zero

          overlapdata[count].ngood[f] = ngdmag
          overlapdata[count].magoff[f] = magoffset
          overlapdata[count].magofferr[f] = magofferr  ; error in the mean

          ; Print info
          if keyword_set(verbose) then $
            print,overlapdata[count].field1,overlapdata[count].field2,overlapdata[count].ngood[f],filters[f],$
                 overlapdata[count].magoff[f],overlapdata[count].magofferr[f],format='(A10,A10,I7,A5,F10.4,F10.5)'
        endif ; some matches

      ; Not both calibrated
      endif else begin
        if keyword_set(verbose) then $
          print,field1,field2,'----',filters[f],$
                '------','-------',format='(A10,A10,A7,A5,A10,A10)'
      endelse
    Endfor  ; filter loop

    ; Increment
    count++

    BOMB2:
  Endfor  ; field inner loop
  ;stop
  BOMB1:
Endfor  ; field outer loop

; Trim excess elements
if n_elements(overlapdata) gt count then begin
  overlapdata = overlapdata[0:count-1]
  noverlapdata = n_elements(overlapdata)
endif

;MWRFITS,overlapdata,dir+galtype+'overlap_magoffsets.fits',/create

doplots:

; Make some plots
!p.font = 0
setdisp,/silent
for f=0,nfilters-1 do begin
  file = dir+galtype+'overlap_magoffhist_'+filters[f]
  magoff = overlapdata.magoff[f]
  gd = where(magoff lt 50,ngd)
  ps_open,file,/color,thick=4,/encap
  bin = 0.005
  hist = histogram(magoff[gd],bin=bin,min=-0.1,max=0.1,locations=xhist)
  xhist += 0.5*bin
  yr = [0,max(hist)*1.2]
  plot,xhist,hist,ps=10,xs=1,ys=1,yr=yr,xtit=filters[f],ytit='N',tit=filters[f]+' magnitude offsets'
  med = median(magoff[gd])
  sig = mad(magoff[gd])
  ; Gaussian area = ht*sig*sqrt(2*pi)
  ht = bin * total(hist) / (sig*sqrt(2*!dpi))
  x = scale_vector(findgen(1000),-0.2,0.2)
  oplot,x,gfunc(x,[ht,med,sig]),co=250
  al_legend,['Median offset = '+stringize(med,ndec=3)+' mag','Sigma offset = '+stringize(sig,ndec=3)+' mag'],/top,/left,charsize=1.5
  ps_close
  ps2png,file+'.eps',/eps
  spawn,['epstopdf',file+'.eps'],/noshell
endfor
cd,current=curdir
cd,dir
spawn,'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+galtype+'overlap_magoffhist.pdf '+galtype+'overlap_magoffhist_?.pdf'
cd,curdir

;stop
skiptohere:

; Initialize the final output structure
overlapstr = {field:strtrim(fieldoverlapstr[gdfields].field,2),noverlap:lonarr(ngdfields),$
              ind0:lonarr(ngdfields)-1,ind1:lonarr(ngdfields)-1,$
              index:lonarr(noverlapdata),revindex:lonarr(noverlapdata),data:overlapdata}
; ind0 and ind1 are the index into "index" and "revindex" which in
; turn are indices into the overlapdata structure,
; e.g. data[index[ind0[i]:ind[i]]] gives the elements of the overlap
; data structure for the ith chip.

; Get reverse indices
oid1 = strtrim(overlapdata.field1,2)
oid2 = strtrim(overlapdata.field2,2)
index_cnt = 0LL
for i=0,ngdfields-1 do begin
  MATCH,oid1,strtrim(fieldoverlapstr[gdfields[i]].field,2),ind1a,ind1b,count=nmatch1,/sort
  ; reverse ones
  MATCH,oid2,strtrim(fieldoverlapstr[gdfields[i]].field,2),ind2a,ind2b,count=nmatch2,/sort
  overlapstr.noverlap[i] = nmatch1
  if nmatch1 ne nmatch2 then stop,'DIFFERENT MATCHES!!'
  if nmatch1 gt 0 then begin
    overlapstr.ind0[i] = index_cnt
    overlapstr.ind1[i] = index_cnt+nmatch1-1
    overlapstr.index[index_cnt:index_cnt+nmatch1-1] = ind1a
    overlapstr.revindex[index_cnt:index_cnt+nmatch1-1] = ind2a
    index_cnt += nmatch1
  endif
endfor

; Print out some statistics of the offsets
gd = where(overlapstr.data.magoff lt 50,ngd)
;gd = where(overlapstr.data.magoff lt 50 and overlapstr.data.overlap eq 1 and overlapstr.data.primary eq 1,ngd)
if not keyword_set(silent) then begin
  print,ngd,' mag offsets','med=',median(overlapstr.data[gd].magoff),'rms=',$
        mad(overlapstr.data[gd].magoff),format='(I8,A12,A5,F11.6,A5,F11.6)'
  print,'',' errors','med=',median(overlapstr.data[gd].magofferr),'rms=',$
        mad(overlapstr.data[gd].magofferr),format='(A8,A12,A5,F11.6,A5,F11.6)'
endif

; Save the results
print,'Saving results'
MWRFITS,overlapdata,dir+galtype+'overlap_magoffsets.fits',/create
SAVE,overlapstr,file=dir+galtype+'overlaps.dat'

; Solve the ubercal problems
SMASHRED_SOLVE_UBERCAL_MAINBODY,overlapstr,ubercalstr,verbose=verbose
MWRFITS,ubercalstr,dir+galtype+'ubercal.fits',/create

stop

end

