#! /usr/bin/env python

import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord

from skimage import measure
from skimage.segmentation import find_boundaries

from scipy.special import erf

import pyregion

from fluxtools import fitsutils

import logging
logging.basicConfig(format="%(levelname)s (%(module)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# https://stackoverflow.com/a/61629292/3293881 @Ehsan
def norm_method(arr, point):
    point = np.asarray(point)
    return np.linalg.norm(np.indices(arr.shape, sparse=True)-point)


def kvis_polygon(region, r_index=0):
    """Get coordinate list for polygon (CLINES) in annotation file.

    """

    f = open(region, "r")
    lines = f.readlines()

    polygons = []
    for line in lines:
        if "clines" in line.lower():
            polygons.append(line)

    poly = polygons[r_index]  # If there are more than one cline entries.

    bits = poly.split(" ")
    coords = []
    for i in range(1, len(bits)-1, 2):
        coords.append((float(bits[i]), float(bits[i+1])))

    f.close()

    return coords


def is_in(x, y, region_x, region_y, max_x):
    """Check if (x, y) is in region.

    (x, y) is considered in region if the followed is met:
    A line drawn horizontally to the right from (x, y) intersects
    with an odd number of region edges.

    # This adapts a C++ implementation found here:
    # http://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/

    """



    def orientation(p1, p2, p3):
        """Get orientation of ordered triplet.

        0 == collinear
        1 == clockwise
        2 == counter-clockwise
        """

        slope = (p2[1] - p1[1]) * (p3[0] - p2[0]) - \
                (p3[1] - p2[1]) * (p2[0] - p1[0])

        if slope < 0:
            orient = 2
        elif slope > 0:
            orient = 1
        else:
            orient = 0


        return orient


    def on_segment(p1, p2, p3):
        """Checks if p3 lies on segment (p1, p2)."""

        if ((p3[0] <= max(p1[0], p2[0])) and (p3[0] >= min(p1[0], p2[0])) and\
            (p3[1] <= max(p1[1], p2[1])) and (p3[1] >= min(p1[1], p2[1]))):

            return True

        else:

            return False

    def intersects(p1, q1, p2, q2):
        """Determine if line segment (p1, q1) intersects with (p2, q2)

        Uses orientation conditions.
        """

        inter = False

        # General case:
        if (orientation(p1, q1, p2) != orientation(p1, q1, q2)) and \
            (orientation(p2, q2, p1) != orientation(p2, q2, q1)):

            inter = True

        # Special cases:
        if (orientation(p1, q1, p2) == 0) and (on_segment(p1, q1, p2)):
            inter = True
        if (orientation(p1, q1, q2) == 0) and (on_segment(p1, q1, q2)):
            inter = True
        if (orientation(p2, q2, p1) == 0) and (on_segment(p2, q2, p1)):
            inter = True
        if (orientation(p2, q2, q1) == 0) and (on_segment(p2, q2, q1)):
            inter = True


        return inter


    coords_in_region = []

    for i in range(len(x)):

        p1, q1 = (x[i], y[i]), (x[i]+max_x, y[i])

        intersections = 0

        for j in range(len(region_x)):
            p2 = (region_x[j], region_y[j])
            if j == len(region_x)-1:
                q2 = (region_x[0], region_y[0])
            else:
                q2 = (region_x[j+1], region_y[j+1])

            if intersects(p1, q1, p2, q2):
                intersections += 1

        if intersections%2 == 0:
            in_region = False
        else:
            in_region = True


        coords_in_region.append(in_region)

    return coords_in_region



def get_source_pixels_in_aperture(ra, dec, radius, fitsobj):
    scale = int((radius / abs(fitsobj.cdelt1)))
    xr, yr = fitsobj.world2pix(ra, dec)
    dist = norm_method(fitsobj.data, (xr, yr))

    cond1 = dist <= scale
    indices_inside = np.where(cond1)

    return indices_inside


def get_source_pixels_in_polygon(ra_vert, dec_vert, fitsobj, pix_buff=0):
    x_vert, y_vert = fitsobj.world2pix(ra_vert, dec_vert, no_int=True)
    xlim1, xlim2 = min(x_vert), max(x_vert)
    ylim1, ylim2 = min(y_vert), max(y_vert)
    xlim1 = int(xlim1-1) - pix_buff
    xlim2 = int(xlim2+1) + pix_buff
    ylim1 = int(ylim1-1) - pix_buff
    ylim2 = int(ylim2+1) + pix_buff
    max_x = max(x_vert)
    if xlim1 < 0:
        xlim1 = 0
    if ylim1 < 0:
        ylim1 = 0
    if xlim2 > len(fitsobj.data[:, 0]):
        xlim2 = len(fitsobj.data[:, 0])
    if ylim2 > len(fitsobj.data[:, 1]):
        ylim2 = len(fitsobj.data[:, 1])

    x_all, y_all = [], []
    for i in range(xlim1, xlim2):
        for j in range(ylim1, ylim2):
            x_all.append(i)
            y_all.append(j)
    
    cir = is_in(x_all, y_all, x_vert, y_vert, max_x)

    return np.asarray(x_all)[cir], np.asarray(y_all)[cir]


def get_source_pixels_from_region(region, fitsobj, r_index=0):

    

    r = pyregion.open(region)[r_index]

    if r.name == "circle":
        ra, dec, radius = r.coord_list
        indices_x, indices_y = get_source_pixels_in_aperture(ra, dec, radius, fitsobj)


    elif r.name == "polygon":

        ra, dec = [], []

        for i in range(0, len(r.coord_list), 2):
            ra.append(r.coord_list[i])
            dec.append(r.coord_list[i+1])

        indices_x, indices_y = get_source_pixels_in_polygon(
            ra_vert=ra,
            dec_vert=dec,
            fitsobj=fitsobj
        )

    else:
        raise RuntimeError("Item {} in {} is not a circle or polygon.".format(
            r_index, region
        ))

    return indices_x, indices_y


def extract_largest_angular_scale(fitsobj, indices_x, indices_y, sigma=None):

    labelled_array = fitsobj.data.copy() 
    labelled_array[:] = 0
    labelled_array[indices_x, indices_y] = 1
    if sigma is not None and fitsobj.rms is not None:
        labelled_array[np.where(fitsobj.data < sigma*fitsobj.rms)] = 0

    labelled_array[np.where(np.isnan(fitsobj.data))] = 0

    boundary = find_boundaries(labelled_array.astype(int), mode="outer").astype(int)
    bx, by = np.indices(boundary.shape)

    boundary_indices = np.where(boundary.flatten() == 1)[0]
    bx = bx.flatten()[boundary_indices]
    by = by.flatten()[boundary_indices]

    las = 0.

    r, d = fitsobj.pix2world(bx, by)
    max_coords_ra = [np.nan, np.nan]
    max_coords_dec = [np.nan, np.nan]
    coords = SkyCoord(ra=r*u.deg, dec=d*u.deg)
    for i in range(len(coords)):
        
        seps = coords[i].separation(coords)
        idx = np.argsort(seps)
        seps_sorted = seps[idx]
        if seps_sorted.value[-1] > las:
            las = seps_sorted.value[-1]
            max_coords_ra[0] = coords[i].ra.value
            max_coords_dec[0] = coords[i].dec.value
            max_coords_ra[1] = coords[idx][-1].ra.value
            max_coords_dec[1] = coords[idx][-1].dec.value


    return las, max_coords_ra, max_coords_dec



def extract_centroid_radec(fitsobj, source_flux, indices_x, indices_y):
    """Get emission-weighted centroid coordinates."""

    fz = np.nansum(source_flux)
    fx = indices_x*source_flux
    fy = indices_y*source_flux
    centroid_x = np.nansum(fx) / fz
    centroid_y = np.nansum(fy) / fz

    centroid_ra, centroid_dec = fitsobj.pix2world(centroid_x, centroid_y)

    return centroid_ra, centroid_dec






def extract_flux_density(fitsobj, indices_x, indices_y, sigma=3,
    full_region_error=False,
    blob_correct=False,
    outname=None,
    do_las=False):
    """
    """

    source_flux = fitsobj.data[indices_x, indices_y].flatten()
    source_bpp = fitsobj.bpp[indices_x, indices_y].flatten()
    source_solid_angle = fitsobj.solid_angle[indices_x, indices_y].flatten()
    source_sb = source_flux/(source_solid_angle*3600.*3600.)

    if fitsobj.rms is not None:

        source_rms = fitsobj.rms[indices_x, indices_y].flatten()


        source_rms_full = source_rms.copy()
        source_bpp_full = source_bpp.copy()

        n_good_pixels = len(np.where(source_flux >= source_rms*sigma)[0])
        idx = np.where(source_flux < source_rms*sigma)
        
        source_bpp[idx] = np.nan
        source_rms[idx] = np.nan
        source_sb[idx] = np.nan
        source_flux[idx] = np.nan
        indices_x_clipped = indices_x.copy().astype(float)
        indices_y_clipped = indices_y.copy().astype(float)
        indices_x_clipped[idx] = np.nan
        indices_y_clipped[idx] = np.nan
        source_rms_sb = np.nanmean(source_rms/source_solid_angle)/(3600.*3600.)
    else:
        n_good_pixels = len(source_flux)
        source_rms_sb = 0.
        indices_x_clipped = indices_x
        indices_y_clipped = indices_y

    centroid_ra, centroid_dec = extract_centroid_radec(
        fitsobj=fitsobj,
        source_flux=source_flux, 
        indices_x=indices_x_clipped,
        indices_y=indices_y_clipped
    )




    source_int_flux = np.nansum(source_flux / source_bpp)
    # source_sb = (source_flux/source_bpp)/abs(fitsobj.cdelt1*fitsobj.cdelt2*3600.*3600.)
    
    # source_avg_sb = np.nanmean(
        # source_flux/source_bpp
    # )/abs(fitsobj.cdelt1*fitsobj.cdelt2*3600.*3600.)
    source_avg_sb = np.nanmean(source_sb)
    source_peak_sb = np.nanmax(source_sb)
    
    source_peak_flux = np.nanmax(source_flux)

    source_avg_flux = np.nanmean(source_flux)
    
    if fitsobj.rms is not None:
        if full_region_error:
            source_unc_flux = (np.nansum(source_rms_full / source_bpp_full) * \
                np.sqrt(np.nanmean(source_bpp_full) / float(len(source_rms_full))))
            source_rms_avg = np.nanmean(source_rms_full)
        else:
            source_unc_flux = (np.nansum(source_rms / source_bpp) * \
                np.sqrt(np.nanmean(source_bpp) / float(len(source_rms))))
            source_rms_avg = np.nanmean(source_rms)

        if blob_correct:
            # https://doi.org/10.1111/j.1365-2966.2012.21373.x
            snr = source_peak_flux / source_rms_avg
            eta = (erf(np.sqrt(-np.log((sigma/snr)))))**2
            source_int_flux /= eta

    else:
        source_unc_flux = np.nan
        source_rms_avg = np.nan

    
    source_area = n_good_pixels*abs(fitsobj.cdelt1*fitsobj.cdelt2)

    if do_las:
        source_las, max_ra, max_dec = extract_largest_angular_scale(
            fitsobj=fitsobj, 
            indices_x=indices_x, 
            indices_y=indices_y,
            sigma=sigma
        )
    else:
        source_las, max_ra, max_dec = np.nan, np.nan, np.nan
    
    # print(source_int_flux / (source_area*3600.*3600.))
    # source_peak_sb = (source_peak_flux/np.nanmean(source_bpp)) /abs(fitsobj.cdelt1*fitsobj.cdelt2*3600.*3600.)

    params = [
        source_int_flux,  # 0
        source_peak_flux, # 1
        source_unc_flux,  # 2
        source_rms_avg,   # 3
        source_area,      # 4
        source_las,       # 5
        max_ra,           # 6
        max_dec,          # 7
        source_avg_sb,    # 8
        source_peak_sb,   # 9
        source_avg_flux,  # 10
        source_rms_sb,    # 11
        centroid_ra,      # 12
        centroid_dec      # 13
    ]
    
    return params
    

def measure_flux_density(fitsimage, rms, region=None,
    coords=None,
    radius=None,    
    r_index=0, 
    sigma : float = 3,
    verbose=False,
    absolute_flux=None,
    modelimage=None,
    psfimage=None,
    residualimage=None,
    dOmegaimage=None,
    full_region_error=False,
    blob_correct=False,
    do_las=False,
    outname=None):
    """
    
    Measure integrated flux within polygon region.


    Requires pyregion if using ds9 region file. Inspired by `radioflux.py` by
    Martin Hardcastle: https://github.com/mhardcastle/radioflux


    
    """

    if verbose:
        logger.info("Working on {}".format(fitsimage))
    fitsobj = fitsutils.Fits(fitsimage)
    fitsobj.add_beam(psfimage=psfimage) 
    fitsobj.add_rms(rms=rms) 
    fitsobj.add_beams_per_pixel()


    # get source pixels
    if region is not None:
        indices_x, indices_y = get_source_pixels_from_region(
            region=region,
            fitsobj=fitsobj, 
            r_index=r_index
        )

    elif coords is not None and radius is not None:
        indices_x, indices_y = get_source_pixels_in_aperture(
            ra=coords[0], 
            dec=coords[1],
            radius=radius, 
            fitsobj=fitsobj
        )
    else:
        # total flux in image
        indices_x, indices_y = np.indices(fitsobj.data.shape)
    

    params = extract_flux_density(
        fitsobj=fitsobj, 
        indices_x=indices_x, 
        indices_y=indices_y, 
        sigma=sigma,
        full_region_error=full_region_error,
        blob_correct=blob_correct,
        do_las=do_las
    )

    if verbose:
        logger.info("Region {}: flux density       = {} ({}) Jy".format(r_index, params[0], params[2]))
        logger.info("Region {}: peak flux (rms)    = {} ({}) Jy/beam".format(r_index, params[1], params[3]))
        logger.info("Region {}: average flux       = {} Jy/beam".format(r_index, params[10]))
        logger.info("Region {}: surface brightness = {} Jy/arcsec^2".format(
            r_index, params[8]
        ))
        logger.info("Region {}: peak SB            = {} ({}) Jy/arcsec^2".format(
            r_index, params[9], params[11]
        ))
        if do_las:
            logger.info("Region {}: max LAS            = {} arcmin".format(r_index, params[5]*60.))
        logger.info("Region {}: Centroid           = {}, {}".format(
            r_index, params[12], params[13]
        ))

    if outname is not None:
        fitsobj.writeout_source(indices_x, indices_y, outname)

    return params


def measure_all_flux_density(fitsimage, rms, region,
    sigma=2,
    verbose=True,
    absolute_flux=False,
    psfimage=None,
    modelimage=None,
    residualimage=None,
    full_region_error=False,
    blob_correct=False):
    """Measure all regions within a file.
    """

    fitsobj = fitsutils.Fits(fitsimage)
    fitsobj.add_beam(psfimage=psfimage) 
    fitsobj.add_rms(rms=rms) 
    fitsobj.add_beams_per_pixel()

    region_data = pyregion.open(region)
    
    total_flux = 0
    total_unc = 0
    peak_flux = 0
    total_area = 0
    average_rms = np.nan

    for i in range(len(region_data)):


        try:
            indices_x, indices_y = get_source_pixels_from_region(
                region=region,
                fitsobj=fitsobj, 
                r_index=i
            )
            params = extract_flux_density(
                fitsobj=fitsobj, 
                indices_x=indices_x, 
                indices_y=indices_y, 
                sigma=sigma,
                full_region_error=full_region_error,
                blob_correct=blob_correct,
                do_las=False
            )

            if verbose:
                logger.info("Region {}: flux density       = {} ({}) Jy".format(i, params[0], params[1]))
                logger.info("Region {}: peak flux          = {} ({}) Jy/beam".format(i, params[1], params[2]))
                logger.info("Region {}: surface brightness = {} Jy/arcsec^2".format(
                    i, params[0]/(params[4]*3600.*3600.)
                ))

            total_flux += params[0]
            total_unc = np.sqrt(total_unc**2 + params[2]**2)
            if params[1] > peak_flux:
                peak_flux = params[1]
            average_rms = np.nanmean([average_rms, params[3]])

            total_area += params[4]
        except ValueError:
            logger.warn("Region {} has no detected source.".format(i))

    if verbose:
        logger.info("Total flux density = {} ({}) Jy".format(total_flux, total_unc))
        logger.info("Peak flux          = {} ({}) Jy/beam".format(peak_flux, average_rms))
        logger.info("Surface brightness = {} Jy/arcsec^2".format(
                i, total_flux/(total_area*3600.*3600.)
            ))

    return [total_flux, peak_flux, total_unc, average_rms]

    
