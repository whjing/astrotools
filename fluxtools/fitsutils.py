#! /usr/bin/env python

import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u

try:
    FLOAT_TYPES = [float, np.float, np.float32, np.float64]
except AttributeError: # AttributeError 是 Python 中的一种异常类型，用于指示某个对象（通常是一个模块、类或实例）没有期望的属性。
    FLOAT_TYPES = [float, np.float_, np.float32, np.float64]

class Fits(object):
    """A class to hold and find relevant FITS image information.

    Here we are generally only interested in 2-dimensional images, as we 
    are usually interested only in the continuum flux density of sources.

    """
    #  构造函数（constructor）。它是在创建类的实例时自动调用的方法，用于初始化该实例的属性。



    def __init__(self, fitsimage, extension=0):
        """Open and/or assign a FITS image a class.

        Parameters
        ----------
        fitsimage : str or fits.HDUList object
            Either a filepath or an already opened/created fits.HDUList.
        extension : int, optional
            Extension in the HDU to get header/data information. [Default 0]
        squeeze : bool, optional
            Select True if wanting to remove extra axes. [Default True]

        This initialises the Fits object. This class in general is used to make
        openining FITS files and getting relevant information for flux density
        measurements easier for code/users. At least it seems that way to me!

        """
        
        if isinstance(fitsimage, str): # 检查fitsimage是否为字符串实例
            # Open FITS from file:
            hdu = fits.open(fitsimage)
            opened = True
        elif isinstance(fitsimage, fits.HDUList): # 检查fitsimage是否为 fits.HDUList 类型
            # Already opened:
            hdu = fitsimage
            opened = False
        else:
            raise IOError("`fitsimage` must be one of: str or fits.HDUList.")

        self.hdr = strip_wcsaxes(hdu[extension].header)
        self.data = np.squeeze(hdu[extension].data)
        self.wcs = WCS(hdu[extension].header).celestial #A copy of the current WCS with only the celestial axes included.

        if "CDELT1" in self.hdr.keys():
            self.cdelt1 = self.hdr["CDELT1"]
            self.cdelt2 = self.hdr["CDELT2"]
        elif "CD1_1" in self.hdr.keys():
            self.cdelt1 = self.hdr["CD1_1"]
            self.cdelt2 = self.hdr["CD2_2"]
        else:
            raise ValueError("Pixel sizes cannot be determined.")

        if opened: #表达式会自动判断变量的真假，不需要显式地与 True 或 False 进行比较
            hdu.close()

        self.bmaj = None
        self.bmin = None



    def pix2world(self, x, y):
        """Wrapper for a wrapper of all_pix2world.""" #"all_pix2world"的包装器。
        return pix_to_world(self.wcs, x, y)



    def world2pix(self, ra, dec, no_int=False):
        """Wrapper of a wrapper of all_world2pix."""
        x, y =  world_to_pix(self.wcs, ra, dec)
        # if not no_int:
        #     try:
        #         x = int(x)
        #         y = int(y)
        #     except TypeError:
        #         x = x.astype("i")
        #         y = y.astype("i")

        return x, y
    

    def add_beams_per_pixel(self):
        """Calculate the fraction of the beam a pixel occupies.

        Used for calculation of integrated flux density over a number of pixels.
        """

        if self.bmaj is None:
            self.add_beam()
        self.bpp = self.solid_angle / abs(self.cdelt1*self.cdelt2)


    def add_beam(self, bmaj=None, bmin=None, psfimage=None):
        """Add beam information manually."""

        if psfimage is None and bmaj is None:
            self.find_beam()
            self.bmaj = np.full_like(self.data, self.bmaj)
            self.bmin = np.full_like(self.data, self.bmin)
        elif psfimage is not None:
            psfdata = fits.getdata(psfimage)
            self.bmaj = psfdata[0, :, :]
            self.bmin = psfdata[1, :, :]
        else:
            # full_like 创建一个与给定数组具有相同形状的新数组，并用指定的常数值填充。
            self.bmaj = np.full_like(self.data, bmaj)
            self.bmin = np.full_like(self.data, bmin)

        # 立体角的计算方法： ?
        self.solid_angle = (np.pi * self.bmaj * self.bmin) / (4. * np.log(2.))
        



    def find_beam(self):
        """Find restoring/synthesized beam information."""

        if ("BMAJ" in self.hdr.keys()) and ("BMIN" in self.hdr.keys()):
        
            self.bmaj = self.hdr["BMAJ"]
            self.bmin = self.hdr["BMIN"]
        
        elif ("CLEANBMJ" in self.hdr.keys()) and ("CLEANBMN" in self.hdr.keys()):
        
            self.bmaj = self.hdr["CLEANBMJ"]
            self.bmin = self.hdr["CLEANBMN"]

        # Some specific surveys have specific beam functions:
        # Check for NVSS:
        # 将字符串中的所有字符转换为小写字母形式。
        elif  "nvss" in repr(self.hdr).lower():

            self.bmaj = 45./3600.
            self.bmin = 45./3600.

        # Check for SUMSS:
        elif "sumss" in repr(self.hdr).lower():

            self.bmaj = 45./3600.
            self.bmin = self.sumss_minor(self.hdr["CRVAL2"])/3600.

        # Check for TGSS:
        elif "tgss" in repr(self.hdr).lower():

            self.bmaj = 25./3600.
            self.bmin = self.tgss_minor(self.hdr["CRVAL2"])/3600.


        # ？ 这里为什么要跳过呢
        elif "bmaj" in repr(self.hdr).lower():
            pass

        if self.bmin is None:
            raise ValueError("No beam information found - try adding manually.")



    def add_rms(self, rms):
        """Add an rms array to the Fits object."""

        self.rms = None

        if rms is None:
            self.rms = np.full_like(self.data, 0.)
        else:
            # 先试图读取一个rms值，再试图从rms map中读取rms
            try:
                rms = float(rms)
                self.rms = np.full_like(self.data, rms)
            except ValueError:
                self.rms = np.squeeze(fits.getdata(rms))


    def writeout_source(self, indices_x, indices_y, outname):
        mask = self.data.copy()
        mask[indices_x, indices_y] = np.nan
        data = self.data.copy()
        data[~np.isnan(mask)] = np.nan #  mask 的地方为空值，也不用相乘
        fits.writeto(outname, data, self.hdr, overwrite=True)

    # def add_bpp(self):
    #     if self.bmaj is None:
    #         self.add_beam()

    #     self.bpp = get_bpp(self.bmaj, self.bmin, self.cdelt1, self.cdelt2)


        

        



    @staticmethod
    def sumss_minor(declination):
        """Calculate the SUMSS restoring beam for a given declination."""
        return 45./(np.sin(np.radians(declination)))

    @staticmethod
    def tgss_minor(declination):
        """Calculate the TGSS restoring beam for a given declination."""
        if declination >= 19.: 
            return 25.
        else:
            return 25./(np.cos(np.radians(declination - 19.)))


def ensure_array(a):
    """Ensure an object is an array, and if not make it so.

    Parameters
    ----------
    a : any
        This is turned into an array if not one already.

    Returns
    -------
    np.ndarray
        `a` as a np.ndarray.

    """


    if isinstance(a, np.ndarray) or isinstance(a, np.ma.MaskedArray):
        
        pass

    else:

        if hasattr(a, "__iter__") and not isinstance(a, str):
            a = np.asarray(a)
        elif isinstance(a, str):
            a = np.asarray([a])

    return a



def pix_to_world(wcs, x, y):
    """Convert from pixel to world coordinates.

    Parameters
    ----------
    x, y : int or arrays of ints
        x, y pixel coordinates.
    
    Returns
    -------
    array
        Two arrays: one of RA. one of dec. coordinates.

    """

    x, y = ensure_array(x), ensure_array(y)

    # We are only interested in the celestial coordinates:
    return wcs.celestial.all_pix2world(y, x, 0)

def world_to_pix(wcs, ra, dec):
    """Convert from world coordinates to pixel coordinates.

    Parameters
    ----------
    ra, dec : float or arrays of floats
        RA, dec. coordinates.

    Returns
    -------
    array
        Two int arrays for the x and y pixel coordinates.

    """
 
    ra, dec = ensure_array(ra), ensure_array(dec)

    # We are only interested in the celestial coordinates:
    y, x = wcs.celestial.all_world2pix(ra, dec, 0)

    # return x.astype("i"), y.astype("i")
    return x, y


def strip_wcsaxes(hdr):
    """Strip extra axes in a header object."""


    remove_keys = [key+i for key in 
                   ["CRVAL", "CDELT", "CRPIX", "CTYPE", "CUNIT", "NAXIS"]
                   for i in ["3", "4", "5"]]

    for key in remove_keys:
        if key in hdr.keys():
            del hdr[key]

    return hdr

def get_bpp(bmaj, bmin, cd1, cd2):
    return (np.pi*bmaj*bmin) / (4.*abs(cd1*cd2)*np.log(2.))
