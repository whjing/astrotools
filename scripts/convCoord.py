#! /usr/bin/env python

# modified

from argparse import ArgumentParser

from fluxtools.measures import measure_flux_density, measure_all_flux_density

from convCoord import 


def get_args():
    # 前置下划线通常用于指示一个变量是"内部的"或"私有的"。这意味着该变量在模块外部不应该被直接访问或修改。这是一种约定，而不是Python解释器强制执行的规则。

    _description = """
Read the coords from the terminal or file(list).
Four values are printed: J2000 in hmsdms, J2000 in deg, GAL in deg.
"""

    ps = ArgumentParser(description=_description)

    # ps.add_argument("fitsimage", type=str,
    #     help="FITS image to measure flux density from.")
    # ps.add_argument("rms", type=str,
    #     help="Either name of FITS image of spatially varying rms noise or a single value.")
    # ps.add_argument("-r", "--region", type=str, default=None,
    #     help="Name of a region file containing polygon and/or circle regions to measure.")
    # ps.add_argument("-c", "--coords", nargs=2, type=float, default=None,
    #     help="Centre RA, DEC of aperture to measure if not using a region file.")
    ps.add_argument("coords", nargs=2, type=float, default=None,
        help="Input coords.")
    # ps.add_argument("-R", "--radius", type=float, default=None,
    #     help="Radius of aperture to measure if not using a region file.")
    # ps.add_argument("-i", "--r_index", type=int, default=0,
    #     help="Row index in region file of region to measure. 0-indexed (i.e. row 0 is the first row). [Default 0]")
    # ps.add_argument("-s", "-t", "--sigma", "--threshold", dest="sigma", type=float, default=3,
    #     help="SNR to mask pixels. Pixels below `sigma`*`rms` are set to NaN. [Default 3]")
    # ps.add_argument("-p", "--psfimage", default=None, type=str,
    #     help="PSF cube describing a spatially varying PSF. Default is to use the BMAJ, BMIN from the FITS image.")
    # ps.add_argument("-F", "--full_region_error", action="store_true",
    #     help="Determine measurement error from the full region rather than just the measured pixels. This should generally be used but may overestimate the error if there are many blanked pixels within the measured region. ")
    # ps.add_argument("-B", "--blob_correct", action="store_true")
    # ps.add_argument("-v", "--verbose", action="store_true",
    #     help="Print out measurements in a more human-readable form.")
    # ps.add_argument("-o", "--outname", default=None,
    #     help="Set an output image name to write out an image with everything NaN except for measured pixels (e.g. to check region and sigma threshold).")
    # ps.add_argument("-A", "--all", action="store_true", 
    #     help="Measure all regions in a given region file. If --verbose is used, then all components are printed them summed at the end.")

    args = ps.parse_args()

    return args


def cli():

    args = get_args()
    params = convCoord(
        
    )
    # if args.all and args.region is not None:
    #     params = measure_all_flux_density(
    #         fitsimage=args.fitsimage,
    #         rms=args.rms,
    #         region=args.region,
    #         sigma=args.sigma,
    #         verbose=args.verbose,
    #         psfimage=args.psfimage,
    #         full_region_error=args.full_region_error,
    #         blob_correct=args.blob_correct
    #     )
    # else:
    #     params = measure_flux_density(
    #         fitsimage=args.fitsimage,
    #         rms=args.rms,
    #         region=args.region,
    #         coords=args.coords,
    #         radius=args.radius,
    #         r_index=args.r_index,
    #         sigma=args.sigma,
    #         verbose=args.verbose,
    #         psfimage=args.psfimage,
    #         full_region_error=args.full_region_error,
    #         blob_correct=args.blob_correct,
    #         outname=args.outname
    #     )

    print("{} {} {} {}".format(params[0], params[1], params[2], params[3]))


if __name__ == "__main__":
    cli()