;;pro stack_KAH

;; load up the names for the files, positions, and number of objects in each file.
pos_all = ['FWHMCIV_slices/nobals/Stripe82_Helms_FWHMCIV_nobal_group10.txt']
number = 44

stack_sz = 20

;names for printed out images
nam250 =     ['FWHMCIV_slices/nobals/Stripe82_Helms_FWHMCIV_250_nobal.fits']
nam250_err = ['FWHMCIV_slices/nobals/Stripe82_Helms_FWHMCIV_250_nobalerr.fits']


;;; load the images
pic250    = readfits('HELMS_image_250_SANEPIC_v0.2.fits',hd250,ext=1)
wei250 = readfits('HELMS_image_250_SANEPIC_v0.2.fits',hd250,ext=2)


;get coords for 250 for Helms in terms of pixels

EXTAST, hd250, astrhd
pos1 = fltarr(4,number)
openr, 1, (pos_all)
readf,1,pos1
a = pos1(1,*)
d = pos1(2,*)
AD2XY,a,d,astrhd,x,y


;;;;;; FOR 250 for Helms
; gives the size of the subimage, where subimage will be 2*stack_sz + 1
;;; Give the stacked image using the mean and median
dostack_xy_part, pic250, x,y, stpic,wei_1,datcube,size=stack_sz,/DCUBE,WEIGHT=wei250,Nstacked=Nstacked

size1 = Nstacked

stackSize = 2*stack_sz+1
  img_stack  = FLTARR(stackSize, stackSize)
  img_stack = (stpic) / (size1)
  wei_stack =  FLTARR(stackSize, stackSize)
  wei_stack = (wei_1) / (size1)

wei_stdev = STDDEV(wei_stack)

;creates a datacube to be used in the bootstrap analysis
Full_cube = FLTARR(stackSize, stackSize, size1)
    Full_cube[*,*,0:size1-1]=datcube[*,*,*]
;    Full_cube[*,*,size1:(size1 )-1] = datcube2[*,*,*]

;;; write out the stacked postage stamp
fits_write,nam250,img_stack

;edit the header info
h1 = headfits(nam250)
sxaddpar, h1, 'XTENSION','IMAGE   ','IMAGE extension'
sxaddpar, h1, 'PCOUNT',  0, 'required keyword; must = 0'
sxaddpar, h1, 'GCOUNT',  1, 'required keyword; must = 1'
sxaddpar, h1, 'WCSAXES', 2, 'Number of coordinate axes'
sxaddpar, h1, 'NAXIS',   2, 'Numer of data axes'
sxaddpar, h1, 'NAXIS1', stack_sz*2+1, 'length of data axis 1'
sxaddpar, h1, 'NAXIS2', stack_sz*2+1, 'length of data axis 2'
sxaddpar, h1, 'CRPIX1',  1, 'Pixel coordinate of reference point'
sxaddpar, h1, 'CRPIX2',  1,  'Pixel coordinate of reference point'
sxaddpar, h1, 'CDELT1',  -0.00166666,  '[deg] Coordinate increment at reference point'
sxaddpar, h1, 'CDELT2',  0.001666667,  '[deg] Coordinate increment at reference point'
sxaddpar, h1, 'CUNIT1',  'deg', 'Units of coordinate increment and value'
sxaddpar, h1, 'CUNIT2',  'deg', 'Units of coordinate increment and value'
sxaddpar, h1, 'CTYPE1',  'RA---TAN', 'Right ascension, gnomonic projection'
sxaddpar, h1, 'CTYPE2',  'DEC--TAN', 'Declination, gnomonic projection'
sxaddpar, h1, 'CRVAL1',  2.85039189502, '[deg] Coordinate value at reference point'
sxaddpar, h1, 'CRVAL2', -0.235228650344, '[deg] Coordinate value at reference point'
sxaddpar, h1, 'LONPOLE', 180, '[deg] Native longitude of celestial pole'
sxaddpar, h1, 'LATPOLE', -0.235228650344, '[deg] Native latitude of celestial pole'
sxaddpar, h1, 'RESTFRQ', -1, '[Hz] Line rest frequency'
sxaddpar, h1, 'RESTWAV', -1, '[Hz] Line rest wavelength'
sxaddpar, h1, 'CNAME1', 'Right Ascension', 'Axis name for labelling purposes'
sxaddpar, h1, 'CNAME2', 'Declination', 'Axis name for labelling purposes'
sxaddpar, h1, 'RADESYS','ICRS', 'Equatorial coordinate system'
sxaddpar, h1, 'EQUINOX', 2000, '[yr] Equinox of equatorial coordinates'
sxaddpar, h1, 'DATE-OBS', '2012-01-15T22:23:25.000000', 'ISO-8601 observation date'
sxaddpar, h1, 'TIMESYS', 'UTC     ', 'All dates are in UTC time'
sxaddpar, h1, 'CREATOR', 'SPG v8.0.0', 'Generator of this product'
sxaddpar, h1, 'INSTRUME', 'SPIRE   ', 'Instrument attached to this product'
sxaddpar, h1, 'OBJECT', 'L1A-1-2 ', 'Target name'
sxaddpar, h1, 'TELESCOP', 'Herschel Space Observatory', 'Name of telescope';

modfits,nam250,0,h1

;reads image back in and reads in psf image
stacked_img2 = readfits(nam250)

;;;;;Creates psf
psfimg = readfits("psf_250_spire.fits",ext=0)

; original images are in 1arcsec sky bin size
; so need to convert to same as Herschel images
; i.e 250 = 1401arcsec /6pixel = 233.5
psf_scale = 233/2 
psfimg2 = frebin(psfimg,233.5,233.5,/total)

;cuts out correct size
;Must match the size of the stacked image
new_img_size = stack_sz*2 + 1

;Gives central pixel and the distance either side of that pixel
;to cut out
cent = stack_sz
minVal = psf_scale - stack_sz
maxVal = psf_scale + stack_sz

;cuts out the image
newarr= psfimg2[minVal:maxVal,minVal:maxVal]
cent_pix = newarr[cent,cent]
newarr = (newarr/cent_pix) * 1

;converts to float to match the stacked image type
newpsf = FLOAT(newarr)

fits_write,"test_psf_250.fits", newpsf

;fits and gives flux and clsutering
fit_flux_plus_clustering_stacking, stacked_img2, newpsf, source_flux, clustering_signal, bg, gamma = 1.8


;uses the datacube and psf to do a bootstrap,
stacking_bootstrap,Full_cube,eflux,Nboot=200,psf=newpsf

;calculates the error on the population using the equation in Bethermin 2012
;rearranged to give err_pop
sig_IC = wei_stdev^2
sig_boot = eflux^2
N_stack = size1

err_pop = sqrt((sig_boot)*(N_stack) - (sig_IC))


;prints out the important values
print, "the flux (mJy) is ", source_flux*10^3
print, "error (mJy) from bootstrap is ",eflux*10^3
print, "the clustering (mJy) is ", clustering_signal*10^3
print, "the number of objects stacked is ", N_stack

END

