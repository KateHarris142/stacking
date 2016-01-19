;;pro stack_KAH

;; load up the names for the files, positions, and number of objects in each file.
pos_all =  ['newHersHelmsFiles/zslices_HersHelms/Stripe82_Helms_zall.txt']
number = 549


nam250 =     ['newHersHelmsFiles/zslices_HersHelms/Stripe82_Helms_500_testall.fits']
nam250_err = ['newHelmsHersFiles/zslices_HersHelms/Stripe82_Helms_500_testall.txt']

nam250pt1 =     ['newHersHelmsFiles/zslices_HersHelms/Stripe82_Helms_500_testpt1.fits']
nam250_errpt1 = ['newHelmsHersFiles/zslices_HersHelms/Stripe82_Helms_500_testpt1.txt']

nam250pt2 =     ['newHersHelmsFiles/zslices_HersHelms/Stripe82_Helms_500_testpt2.fits']
nam250_errpt2 = ['newHelmsHersFiles/zslices_HersHelms/Stripe82_Helms_500_testpt2.txt']


;;; load the images
pic250   =readfits('HELMS_image_500_SANEPIC_v0.2.fits',hd250,ext=1)
weight250=readfits('HELMS_image_500_SANEPIC_v0.2.fits',hd250,ext=2)


;get coords for 250 for helms
EXTAST, hd250, astrhd
pos1 = fltarr(4,number)
openr, 1, (pos_all)
readf,1,pos1
a = pos1(1,*)
d = pos1(2,*)
AD2XY,a,d,astrhd,x,y

;;;split sample
loop = 5

sf_both = fltarr(loop+1)
sf_p1 = fltarr(loop+1)
sf_p2 = fltarr(loop+1)

for i =0, loop de begin
arr = randomu(seed,number)*number
arr2 = arr[number/2:*]
x2 = x[arr2]
y2 = y[arr2]

arr3 = arr[0:number/2]
x3 = x[arr3]
y3 = y[arr3]



;;;;;; FOR 500 for helms half 1
; gives the size of the subimage, where subimage will be 2*stack_sz + 1
stack_sz = 20
;;; Give the stacked image using the mean and median
dostack_xy_part, pic250, x2,y2, stpic1,wei_1,size=stack_sz,WEIGHT=weight250

stackSize = 2*stack_sz+1
img_stack_p1 = stpic1/wei_1
fits_write,nam250pt1,img_stack_p1 

;edit the header info
h2 = headfits(nam250pt1)
sxaddpar, h2, 'XTENSION','IMAGE   ','IMAGE extension'
sxaddpar, h2, 'PCOUNT',  0, 'required keyword; must = 0'
sxaddpar, h2, 'GCOUNT',  1, 'required keyword; must = 1'
sxaddpar, h2, 'WCSAXES', 2, 'Number of coordinate axes'
sxaddpar, h2, 'CRPIX1',  1, 'Pixel coordinate of reference point'
sxaddpar, h2, 'CRPIX2',  1,  'Pixel coordinate of reference point'
sxaddpar, h2, 'CDELT1',  -0.001333333,  '[deg] Coordinate increment at reference point'
sxaddpar, h2, 'CDELT2',  0.003333333,  '[deg] Coordinate increment at reference point'
sxaddpar, h2, 'CUNIT1',  'deg', 'Units of coordinate increment and value'
sxaddpar, h2, 'CUNIT2',  'deg', 'Units of coordinate increment and value'
sxaddpar, h2, 'CTYPE1',  'RA---TAN', 'Right ascension, gnomonic projection'
sxaddpar, h2, 'CTYPE2',  'DEC--TAN', 'Declination, gnomonic projection'
sxaddpar, h2, 'CRVAL1',  2.85039189502, '[deg] Coordinate value at reference point'
sxaddpar, h2, 'CRVAL2', -0.235228650344, '[deg] Coordinate value at reference point'
sxaddpar, h2, 'LONPOLE', 180, '[deg] Native longitude of celestial pole'
sxaddpar, h2, 'LATPOLE', -0.235228650344, '[deg] Native latitude of celestial pole'
sxaddpar, h2, 'RESTFRQ', -1, '[Hz] Line rest frequency'
sxaddpar, h2, 'RESTWAV', -1, '[Hz] Line rest wavelength'
sxaddpar, h2, 'CNAME1', 'Right Ascension', 'Axis name for labelling purposes'
sxaddpar, h2, 'CNAME2', 'Declination', 'Axis name for labelling purposes'
sxaddpar, h2, 'RADESYS','ICRS', 'Equatorial coordinate system'
sxaddpar, h2, 'EQUINOX', 2000, '[yr] Equinox of equatorial coordinates'
sxaddpar, h2, 'DATE-OBS', '2012-01-15T22:23:25.000000', 'ISO-8601 observation date'
sxaddpar, h2, 'TIMESYS', 'UTC     ', 'All dates are in UTC time'
sxaddpar, h2, 'CREATOR', 'SPG v8.0.0', 'Generator of this product'
sxaddpar, h2, 'INSTRUME', 'SPIRE   ', 'Instrument attached to this product'
sxaddpar, h2, 'OBJECT', 'L1A-1-2 ', 'Target name'
sxaddpar, h2, 'TELESCOP', 'Herschel Space Observatory', 'Name of telescope';

modfits,nam250pt1,0,h2



;;;;;; FOR 500 for Helms pt 2
;;; Give the stacked image using the mean and median
dostack_xy_part, pic250_2, x3,y3, stpic2,wei_2,size=stack_sz,WEIGHT=wei250_2

img_stack_p2 = stpic2/wei_2
fits_write,nam250pt2,img_stack_p2 

;edit the header info
h3 = headfits(nam250pt1)
sxaddpar, h3, 'XTENSION','IMAGE   ','IMAGE extension'
sxaddpar, h3, 'PCOUNT',  0, 'required keyword; must = 0'
sxaddpar, h3, 'GCOUNT',  1, 'required keyword; must = 1'
sxaddpar, h3, 'WCSAXES', 2, 'Number of coordinate axes'
sxaddpar, h3, 'CRPIX1',  1, 'Pixel coordinate of reference point'
sxaddpar, h3, 'CRPIX2',  1,  'Pixel coordinate of reference point'
sxaddpar, h3, 'CDELT1',  -0.001333333,  '[deg] Coordinate increment at reference point'
sxaddpar, h3, 'CDELT2',  0.003333333,  '[deg] Coordinate increment at reference point'
sxaddpar, h3, 'CUNIT1',  'deg', 'Units of coordinate increment and value'
sxaddpar, h3, 'CUNIT2',  'deg', 'Units of coordinate increment and value'
sxaddpar, h3, 'CTYPE1',  'RA---TAN', 'Right ascension, gnomonic projection'
sxaddpar, h3, 'CTYPE2',  'DEC--TAN', 'Declination, gnomonic projection'
sxaddpar, h3, 'CRVAL1',  2.85039189502, '[deg] Coordinate value at reference point'
sxaddpar, h3, 'CRVAL2', -0.235228650344, '[deg] Coordinate value at reference point'
sxaddpar, h3, 'LONPOLE', 180, '[deg] Native longitude of celestial pole'
sxaddpar, h3, 'LATPOLE', -0.235228650344, '[deg] Native latitude of celestial pole'
sxaddpar, h3, 'RESTFRQ', -1, '[Hz] Line rest frequency'
sxaddpar, h3, 'RESTWAV', -1, '[Hz] Line rest wavelength'
sxaddpar, h3, 'CNAME1', 'Right Ascension', 'Axis name for labelling purposes'
sxaddpar, h3, 'CNAME2', 'Declination', 'Axis name for labelling purposes'
sxaddpar, h3, 'RADESYS','ICRS', 'Equatorial coordinate system'
sxaddpar, h3, 'EQUINOX', 2000, '[yr] Equinox of equatorial coordinates'
sxaddpar, h3, 'DATE-OBS', '2012-01-15T22:23:25.000000', 'ISO-8601 observation date'
sxaddpar, h3, 'TIMESYS', 'UTC     ', 'All dates are in UTC time'
sxaddpar, h3, 'CREATOR', 'SPG v8.0.0', 'Generator of this product'
sxaddpar, h3, 'INSTRUME', 'SPIRE   ', 'Instrument attached to this product'
sxaddpar, h3, 'OBJECT', 'L1A-1-2 ', 'Target name'
sxaddpar, h3, 'TELESCOP', 'Herschel Space Observatory', 'Name of telescope';

modfits,nam250pt2,0,h3


;;; put the 2 together to get the final stack.
  stackSize = 2*stack_sz+1
  img_stack  = FLTARR(stackSize, stackSize)
  img_stack = (stpic1 + stpic2) / (wei_1 + wei_2)

;;; write out the stacked postage stamp
fits_write,nam250,img_stack

;edit the header info
h1 = headfits(nam250)
sxaddpar, h1, 'XTENSION','IMAGE   ','IMAGE extension'
sxaddpar, h1, 'PCOUNT',  0, 'required keyword; must = 0'
sxaddpar, h1, 'GCOUNT',  1, 'required keyword; must = 1'
sxaddpar, h1, 'WCSAXES', 2, 'Number of coordinate axes'
sxaddpar, h1, 'CRPIX1',  1, 'Pixel coordinate of reference point'
sxaddpar, h1, 'CRPIX2',  1,  'Pixel coordinate of reference point'
sxaddpar, h1, 'CDELT1',  -0.001333333,  '[deg] Coordinate increment at reference point'
sxaddpar, h1, 'CDELT2',  0.003333333,  '[deg] Coordinate increment at reference point'
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


;;; find the flux using the bethermin code and the Herschel psf
psfimg = readfits("psf_250_spire5.fits")

; for the whole sample
img2 = readfits(nam250)
fit_flux_plus_clustering_stacking, img2, psfimg, source_flux, clustering_signal, bg, gamma = 1.8
sf_both[i] = source_flux

;for part 1
img3 = readfits(nam250pt1)
fit_flux_plus_clustering_stacking, img3, psfimg, source_flux1, clustering_signal, bg, gamma = 1.8
sf_p1[i] = source flux1

;for part 2
img4 = readfits(nam250pt2)
fit_flux_plus_clustering_stacking, img4, psfimg, source_flux2, clustering_signal, bg, gamma = 1.8
sf_p3[i] = source flux2

endfor

print, sf_both
print, sf_p1
print, sf_p2


END
