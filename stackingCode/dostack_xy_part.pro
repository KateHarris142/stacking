pro dostack_xy_part, image, x_arr, y_arr, img_stack,wei_stack,datcube,$
                SIZE=size, SILENT=silent, $
                DCUBE=dcube, WCUBE=wcube, $
                Nstacked=Nstacked, ID_stacked=id_stacked, $
                RESAMPLING=resampling, $
                MEDIAN=MEDIAN, ROTATE=ROTATE, $
                UNDEF=undef, USE_ALL=use_all, $
                WEIGHT=weight

;+
;  NAME: 
;     DOSTACK_XY
;
;  PURPOSE:
;     Do the stacking of a catalog of sources
;
;  EXPLANATION:
;
;  CALLING SEQUENCE: 
;     dostack_xy, image, x, y, imout, ...
;
;  INPUTS:     
;     IMAGE   : input image 
;     X_ARR   : x coordinates of the sources to stack (in pixel, IDL convention), 
;     Y_ARR   : y coordinates of the sources to stack (in pixel, IDL convention), 
;;
;  OPTIONAL INPUT KEYWORDS:
;     SIZE    : Output image will have a size of 2*size+1 x 2*size+1 in pixels.
;               Default value is 30 pix.
;     SILENT  : If it is not set, display warnings.
;     RESAMPLING : Interpolate (bilinearly) subimages in order to be
;                  precisely centered on the sources. 
;     MEDIAN  : Compute median instead of mean. Default is to compute mean.
;     ROTATE  : If set, rotate each source by 90deg * (k MOD 4).
;     UNDEF   : Only subimages without bad pixels are stacked. Values
;               of bad pixels are set with this keyword. Default value
;               is NaN.
;     HEALPIX : if set, the input image is considered as an Healpix
;               vector and coordinates are supposed to be galactic
;               coordinates (in degree).
;     WEIGHT: input map with uncertainties. Must have same size than
;                 the image map. Needed if ERR_STACK and/or WEIGHTCUBE
;                 output keyowrd are set. If set, the stacking will be
;                 weighted by this uncertainties.
;
;  OUTPUTS:   
;     IM_STACKED : 2D array of the result of the stacking. It
;                  corresponds, either to the mean or to the
;                  median. This array has 2*size+1 x 2*size+1 elements.
;
;  OPTIONAL OUTPUT KEYWORDS:
;      NSTACKED : Will contains the number of sources effectively
;                 stacked. It may differ of the number of elements of
;                 input arrays because sources that are close of the edges or out
;                 of the field of view or with bad pixels are
;                 rejected.
;      ID_STACKED : Return the ID of the stacked sources.
;      DCUBE      : Contains the whole cube of data. 3D array which size
;                 is 2*size+1 x 2*size+1 x Nstacked 
;      WCUBE      : Same as datacube but with the wheighting.
;      WEI_STACK : Stacking of weight map.
;
;  PROCEDURE: 
;      ZSCALE_RANGE
;
;  EXAMPLE:
;
;  MODIFICATION HISTORY:
;       Written by NB : 19 jan 2007
;       Add ID_STACKED keyword: NB - 23 mar 2007
;       Change background estimation (median instead of mean): NB - june 2007
;       Add HEALPIX keyword: NB - august 2007
;       Add uncertainties weighted stacking: NB - 3 october 2007
;       Replace NOBACKGROUND keyword by BACKSUBTRACT keyword:
;           NB - 30 october 2007 
;       Save all images in the datacube only if requested or if the
;           MEDIAN keyword is set (to save memory - NB - 11 Feb 2008
;       Rectify err_stack estimation and median estimation - 
;           Nicolas Bavouzet + Alexandre Beelen - 22 Feb 2008
;       Keyword renaming Alexandre Beelen & Matthieu Bethermin, Dec 2009
;
;     Copyright Insitut d'Astrophysique Spatiale
;
;     This program is free software: you can redistribute it and/or modify
;     it under the terms of the GNU General Public License as published by
;     the Free Software Foundation, either version 2 of the License, or
;     (at your option) any later version.
;
;     This program is distributed in the hope that it will be useful,
;     but WITHOUT ANY WARRANTY; without even the implied warranty of
;     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;     GNU General Public License for more details.
;
;     You should have received a copy of the GNU General Public License
;     along with this program.  If not, see <http://www.gnu.org/licenses/>.
;
;-      

  Nsour = N_ELEMENTS(x_arr)   ; number of elements in the positions

; Evaluate size arrays and check keywords
  s = SIZE(image)
  Nx = s(1)
  Ny = s(2)
  print, Nx
  print, Ny
  IF KEYWORD_SET(size) eq 0 THEN size = 30 ; in pixel, how big you want to output image
  IF KEYWORD_SET(weight) THEN BEGIN
     se = SIZE(weight)
     IF se(1) ne Nx or se(2) ne Ny THEN BEGIN    ;Weight map must be same size as image
        MESSAGE, 'Uncertainty map must have same size than data map.'
     ENDIF
  ENDIF


  IF (ARG_PRESENT(wei_stack) OR ARG_PRESENT(wcube)) AND $   ;need a weight map to have wcube
     ( NOT KEYWORD_SET(weight) )  THEN BEGIN
     MESSAGE, 'Warning : Uncertainty map is not defined.', /CONT
  ENDIF

  IF NOT KEYWORD_SET(weight) THEN $   ;if no weight set, each pixel gets weight of 1
     weight = image * 0. + 1

;print, 'min of weight250=',min(weight)
;print, 'max of weight250=',max(weight)

; ------------------------------------------------------------------------
; Padding of the map by the size of the submap plus one pixel (size
; half a pixel does not exist...)
;The final image will have an even number of pixels. We want the central pixel to be
; thcentre of the image. So have to add one as can't have half pixels.

  paddingSize=size+1       ; Default is 31
  stackSize=2*size+1       ; Default is 61
  
  ;creates an array, in this case, (62+13781, 62+2849) with all pixels containing NaN
  signal_map = FLTARR(Nx+2*paddingSize,Ny+2*paddingSize)+!VALUES.F_NAN
  
  ;creates an array for the weight map but pixels contain zeros. 
  weight_map = FLTARR(Nx+2*paddingSize,Ny+2*paddingSize)

; Puts the image data into the signal_map aray within the limits [31:13811,31:2879] 
; Does the same for weight map
  signal_map[paddingSize:paddingSize+Nx-1,paddingSize:paddingSize+Ny-1] = image 
  weight_map[paddingSize:paddingSize+Nx-1,paddingSize:paddingSize+Ny-1] = weight

; For sake of simplicity replace all the undef values by NaN if needed
  IF KEYWORD_SET(undef) THEN $
     IF FINITE(undef) THEN $
        signal_map[WHERE(signal_map EQ undef)] = !VALUES.F_NAN

; Sometimes, error_map can have Nan Values, but translated into
; weight, it means only 0, but 1./Nan^2 = Nan in IDL
; this just tests for NaN in the weight map and if found subsitutes in zero instead.
  shouldBeZero = WHERE(FINITE(weight_map) EQ 0,nZero)
  IF nZero GT 0 THEN $
     weight_map[shouldBeZero] = 0

; if we want to use ALL the pixel in the map, change the undef pixel
; values and weights
  IF KEYWORD_SET(use_all) THEN BEGIN
     toBeUsed = WHERE(FINITE(signal_map) EQ 0)
     signal_map[toBeUsed] = 0
     weight_map[toBeUsed] = 0
     
  ENDIF
;creates arrays sized 61*61 pixels (ie, centre of star is on pixel 31 with 30 pixels either
;side). Deos this for both an image array and a weight array.
  img_stack  = FLTARR(stackSize, stackSize)
  wei_stack  = FLTARR(stackSize, stackSize)

; If you want a datacube or weightcube, this sets up the 3D arrays for these 
  IF KEYWORD_SET(dcube) OR KEYWORD_SET(median) THEN $
     dcube  = FLTARR(stackSize, stackSize, Nsour)
  IF KEYWORD_SET(wcube) THEN $
     wcube = FLTARR(stackSize, stackSize, Nsour)

;---------

  Nstacked   = 0L
  id_stacked = [0]

; main loop on sources to stack
;loops around each of the images to be stacked.

  FOR k=0L, Nsour-1 DO BEGIN 
     
; check position :
; 
; NB : +0.5 because "pixel numbers refer to the center of the pixel in
; each axis" (Greisen & Calabretta, 2002, A&A) and map[X,Y] treat X
; and Y as floor(X) and floor(Y)
;so this would give values of 0.5, 1.5 etc as this would be the centre of the pixel
     
; (size+1) comes from the padding...
; ( +1 ) comes from the resampling if needed
;the 0.5 gives the pixel centre.
; xarr and yarr are the positions of the objects to be stack.
     xmin = x_arr[k] + 0.5 + paddingSize - size       ; same as adding 1.5
     xmax = x_arr[k] + 0.5 + paddingSize + size + 1   ; same as adding 0.5+2*paddingsize
     ymin = y_arr[k] + 0.5 + paddingSize - size
     ymax = y_arr[k] + 0.5 + paddingSize + size + 1
     
; Sources has to be inside the image with at least the central pixel
; valid
     
; Stupid IDL tries the second test B, if the first one A=1 in (A or B),
; so if B is undefined, it fails....
; For xmin or ymin to be less than zero, xarr would have to -1.0 or less.
; so would show the mapping from wcs to pixels coords has gone wrong
; For xmax and ymax, its asking for if xarr+0.5 GT than Nx, ie if the coord plus 0.5 is 
; greater than the total numner of pixels.
; If either of the above are true give the warning. If not, carry on.

     IF ( xmin LT 0 OR xmax GT Nx+paddingSize+size+1 OR $
          ymin LT 0 OR ymax GT Ny+paddingSize+size+1) THEN  BEGIN
        IF NOT KEYWORD_SET(SILENT) THEN BEGIN
           MESSAGE,/CONTINUE, $
                   STRING(FORMAT='("Warning : Source [", I5, "] (", F6.1,",", F6.1,") is out of field !")',$
                          k, x_arr[k], y_arr[k])
        ENDIF
        CONTINUE
     ENDIF
     
; If, in the weight map, the value in the pixel position is equal to 0, give a warning. 
; 
; Also don't understand why using different conditions for weight and signal

     IF weight_map[x_arr[k] + 0.5 + paddingSize, y_arr[k] + 0.5 + paddingSize] EQ 0 THEN  BEGIN
        IF NOT KEYWORD_SET(SILENT) THEN BEGIN
           MESSAGE,/CONTINUE, $
                   STRING(FORMAT='("Warning : Source [", I5, "] (", F6.1,",", F6.1,") is out of field !")',$
                          k, x_arr[k], y_arr[k])
        ENDIF
        CONTINUE
     ENDIF
     
     
     
; extract a sub-image
     
; *_tmp are 2*size+2 sized for the moment
; Takes a 61*61 pixel subimage from the signal_map version of the image.
; Same for the weight map.    
     img_tmp    = signal_map[xmin:xmax, ymin:ymax]
     wei_tmp = weight_map[xmin:xmax, ymin:ymax]

; A bad pixel is one with 0 in sub-image file.     
     bad_pix = WHERE(FINITE(img_tmp) EQ 0, count)

;If there are any bad pixels, do this. 
; Warning is also triggered if the source is too close to the edge of the image as 
; anything outside the image is set as 0 and so classed as a bad pixel.    
     IF count NE 0 THEN BEGIN
        IF NOT KEYWORD_SET(SILENT) THEN $
           MESSAGE,/CONTINUE, STRING(FORMAT='("Warning : Source [", I5, "] (", F6.1,",", F6.1,") too close to border !")',$
                                     k, x_arr[k], y_arr[k])
        CONTINUE
     ENDIF
     
; interpolate to be centered on the source. Does not work with Healpix
; maps.
; Resampling key word not set here (July2013)     
     IF KEYWORD_SET(resampling) THEN BEGIN
        delta_x = x_arr[k] - FLOOR(x_arr[k])
        delta_y = y_arr[k] - FLOOR(y_arr[k])
        IF delta_x GE 0.5 THEN delta_x = delta_x - 1
        IF delta_y GE 0.5 THEN delta_y = delta_y - 1
        xinterp    = FINDGEN(2*size+1) + delta_x
        yinterp    = FINDGEN(2*size+1) + delta_y
        img_tmp    = INTERPOLATE(img_tmp,    xinterp, yinterp, /GRID)
        wei_tmp = INTERPOLATE(wei_tmp, xinterp, yinterp, /GRID)
     ENDIF ELSE BEGIN
;remove the extra row
        img_tmp    = img_tmp[0:stackSize-1,0:stackSize-1]
        wei_tmp = wei_tmp[0:stackSize-1,0:stackSize-1]
     ENDELSE
     
; Rotate of 90deg ccw if asked
     IF KEYWORD_SET(rotate) THEN BEGIN
        img_tmp    = ROTATE(img_tmp,    k MOD 4)
        wei_tmp = ROTATE(wei_tmp, k MOD 4)
     ENDIF
     
; do the proper weighted stack
     IF MIN(FINITE(wei_tmp)) EQ 0 THEN $
        stop

; img_stack is an empty array of size 61*61, same for wei_stack
     img_stack = img_stack + img_tmp*1.   ;wei_tmp     ; puts image*weigth insto correct sibscript in img_stack
     wei_stack = wei_stack + wei_tmp             ; just puts weight from wei_tmp into wei_stack
; For each loop the next image is added onto the 61*61 array of img_stack.      
     
     
     
     
; If using the keyword DCUBE, which you need to when using the bootstrap program,
; this should create the datacube.
; Also needed if median is set.
; Puts, in each loop, the subimage array. then in the next loop adds the next to the 
; position below.      
     IF KEYWORD_SET(dcube) OR ARG_PRESENT(dcube) OR KEYWORD_SET(median) THEN $
        dcube[*,*,Nstacked]   = img_tmp
     
     IF KEYWORD_SET(wcube) OR ARG_PRESENT(wcube) THEN $
        wcube[*,*,Nstacked] = wei_tmp

; replace nstack with nstack+1 and use that in the next loop.
; The      
     Nstacked   = Nstacked + 1
     id_stacked = [id_stacked, k]
     
  ENDFOR
print, 'this is', Nstacked

; If there are more than one image to stack,
  IF Nstacked NE 0 THEN BEGIN
     
     img_stack  = img_stack
     
     wei_stack  = wei_stack  ; works but I dont know why.
     id_stacked = id_stacked[1:Nstacked]
     
     IF KEYWORD_SET(wcube) THEN $
        wcube = wcube[*,*,0:Nstacked-1]
     IF KEYWORD_SET(dcube) OR KEYWORD_SET(median) THEN $
        dcube = dcube[*,*,0:Nstacked-1]
  ENDIF

datcube = FLTARR(stackSize, stackSize, Nstacked-1)
datcube = dcube

;wei_stack= wei_stack/Nstacked
;print, 'size of dcube second time is',size(dcube)

; Return median (instead of mean) if asked by user
  IF KEYWORD_SET(MEDIAN) THEN BEGIN
     IF KEYWORD_SET(wei_stack) THEN BEGIN
        MESSAGE, 'Warning ! Uncertainty stacking is meaningless with median stacking.', /CONT
     ENDIF
     img_stack = MEDIAN(dcube,DIMENSION=3)
     wei_stack = MEDIAN(wcube,DIMENSION=3)
  ENDIF



END
