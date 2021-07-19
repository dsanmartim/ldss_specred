#!/usr/bin/env python
# -*- coding: utf8 -*-
from ldss_spec.tools import spec_red
import os


def main():

    args = parse_arguments()
    ldss_pipeline = spec_red.SpecRed(args.raw_path, args.red_path, args.obs_log_file,
                                     flat_correction=args.flat_correction, fit_continuum=args.fit_continuum)

    ldss_pipeline.run()


def parse_arguments():
    """
    Parse the argument given by the user in the command line.

    Returns:

        args (Namespace): A namespace containing all the parameters that will be used.

    """
    import argparse

    # Parsing Arguments ---
    parser = argparse.ArgumentParser(description="Spectroscopic Reduction for LDSS3 engineering data")

    parser.add_argument('raw_path', metavar='raw_path', type=str, nargs=1,
                        help="Full path to ccd raw data (e.g. /Users/jamesbond/e191204_05/RAWSPEC)")

    parser.add_argument('red_path', metavar='red_path', type=str, nargs=1,
                        help="Full path to working/reduction directory (e.g. /Users/jamesbond/e191204_05/REDSPEC)")

    parser.add_argument('obs_log_file', type=str, nargs=1,
                        help="Full path with the name of the json file containing the observation log information"
                             " (e.g. /Users/jamesbond/e191204_05/obs_log_19_12_04.json)")

    # if the paremeter action is store_true but is missing in the command line, the default value will be False
    # if the paremeter action is store_false but is missing in the command line, the default value will be True

    parser.add_argument("--flat_correction", action="store_false", help="Apply flat correction to science data")

    parser.add_argument("--fit_continuum", action="store_true", help="Run iraf continuum task to the wavelength "
                                                                     "calibrated spectrum and save the fitted continuum")

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    main()
