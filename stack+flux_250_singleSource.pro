;;pro stack_KAH

;; load up the names for the files, positions, and number of objects in each file.



stack_sz = 20
;names for printed out images
nam250 = ['Stipe82_HELMS_HERS_250_zall_det.fits']
nam250_err = ['Stipe82_HELMS_HERS_250_zallerr_det.fits']

;;; load the images
pic250    = readfits('HELMS_image_250_SANEPIC_v0.2.fits',hd250,ext=1)
weight250 = readfits('HELMS_image_250_SANEPIC_v0.2.fits',hd250,ext=2)

;get coords for 250 for Helms in terms of pixels
EXTAST, hd250, astrhd
a = 354.70319
d = -6.51855
AD2XY,a,d,astrhd,x,y

;;;;;; FOR 250 for Helms
; gives the size of the subimage, where subimage will be 2*stack_sz + 1
;;stack_sz = 20

new_img_size = stack_sz*2 + 1

;Gives central pixel and the distance either side of that pixel
;to cut out
cent = stack_sz
minValx = x - stack_sz
minValy = y - stack_sz
maxValx = x + stack_sz
maxValy = y + stack_sz

;cuts out the image
stpic1 = pic250[minValx:maxValx,minValy:maxValy]

;;; write out the stacked postage stamp
fits_write,nam250,stpic1

;edit the header info
h1 = headfits(nam250)
sxaddpar, h1, 'XTENSION','IMAGE   ','IMAGE extension'
sxaddpar, h1, 'PCOUNT',  0, 'required keyword; must = 0'
sxaddpar, h1, 'GCOUNT',  1, 'required keyword; must = 1'
sxaddpar, h1, 'WCSAXES', 2, 'Number of coordinate axes'
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
; i.e 250 = 1401/6
psfimg2 = frebin(psfimg,233.5,233.5,/total)

;cuts out correct size
;Must match the size of the stacked image
new_img_size = stack_sz*2 + 1

;Gives central pixel and the distance either side of that pixel
;to cut out
cent = stack_sz
minVal = (233/2)- stack_sz
maxVal = (233/2)+ stack_sz

;cuts out the image
newarr= psfimg2[minVal:maxVal,minVal:maxVal]
cent_pix = newarr[cent,cent]
newarr = (newarr/cent_pix) * 1

;converts to float to match the stacked image type
newpsf = FLOAT(newarr)

fits_write,"test_psf_250.fits", newpsf

;fits and gives flux and clsutering
fit_flux_plus_clustering_stacking, stacked_img2, newpsf, source_flux, clustering_signal, bg, gamma = 1.8

;prints out the important values
print, "the flux is ", source_flux

END
