#!/usr/bin/env python
# -*- coding: utf8 -*-

# Loading a few python packages
import os
import glob
import warnings
from astropy import log
from astropy.io import fits as pyfits
import json

# Loading iraf packages
from pyraf import iraf
from pyraf.iraf import onedspec
from pyraf.iraf import twodspec, apextract

class SpecRed():
    """
    This class runs a few methods over a data set defined in an python dictionary.
    """

    def __init__(self, raw_path, red_path, obs_log_file, flat_correction=True, fit_continuum=False):

        project_dir = os.path.abspath(os.path.dirname(os.path.dirname((os.path.dirname(__file__)))))
        os.system('export PYTHONPATH={}'.format(project_dir))

        self.data_dir = '{}/ldss_spec/data/'.format(project_dir)

        # Get code path

        # Getting raw and red paths
        self.raw_path = str(os.path.join(raw_path[0], ''))
        self.red_path = str(os.path.join(red_path[0], ''))

        self.obs_log_file = obs_log_file
        self.fit_continuum = fit_continuum

        # Getting observation log  dictionary + sky definition
        self.obs_dict = self.get_obs_log(self.obs_log_file[0])

        self.database_dir = str(os.path.join(self.red_path+'database', ''))

        if fit_continuum:
            self.fit_continuum = True

        if flat_correction:
            self.flat_correction = True


        ## Getting disp solution dictionary.
        # Keys are the names of the files inside tye dispsol folder
        # Values are the names that files will have inside the database folder
        self.disp_solution_files = {'vph_all_open':  'id_vph_all_open.0001',
                                    'vph_blue_open': 'id_vph_blue_open.0001',
                                    'vph_red_og590': 'id_vph_red_og590.0001'}

        # Sci extraction parameters
        self.spec_ext_flags = {'format': 'onedspec', 'interactive': 'yes', 'find': 'yes', 'nfind': '1', 'resize': 'no',
                               'recenter': 'yes', 'edit': 'yes', 'nsum': '10', 'width': '10', 'lower': '-10',
                               'upper': '10', 'b_naverage': '-100', 'extract': 'yes', 'extras': 'no',
                               'trace': 'yes', 'fittrace': 'yes', 't_function': 'chebyshev', 't_order': '3',
                               't_nsum': '15', 't_step': '15', 'background': 'fit', 'b_function': 'chebyshev',
                               'b_order': '1', 'clean': 'yes', 'readnoise': 'enoise', 'gain': 'egain',
                               'weights': 'variance'}

        # Arc extraction parameters
        self.arc_ext_flags = {'format': 'onedspec', 'interactive': 'no', 'find': 'no', 'resize': 'no',
                              'recenter': 'no', 'edit': 'no', 'lower': '-10.0', 'upper': '10.0',
                              'extract': 'yes', 'extras': 'no', 'trace': 'no', 'fittrace': 'no',
                              'background': 'none', 'review': 'yes'}

        # Identify parameters
        self.identify_flags = {'ftype': 'emission', 'fwidth': '40', 'cradius': '30', 'thresho': '10.', 'minsep': '10.'}

        self.wrange = { 'vph_all': [3900.0, 10300.0], 'vph_red': [6050.0, 10300.0], 'vph_blue': [3850.0, 6580.0] }

        # About warnings
        warnings.filterwarnings('ignore')
        log.propagate = False

    def run(self, *args, **kwargs):

        # cleaning up the reduction dir
        self.clean_path(path=self.red_path, ext_list=['fits', 'log', 'txt', 'dat'])

        # Loading dispsol content dictionary
        self.disp_sol_content_dictionary = self.get_dict_of_dispsol_content(self.disp_solution_files,
                                                                            dispsol_dirname='dispsol/')

        # Changing to red_path dir
        self.change_dir(path=self.red_path)

        # Creating directory database if does not exist
        if not os.path.exists(self.database_dir):
            os.makedirs(str(self.database_dir))

        # Extracting spectra and corresponding arcs
        self.extract_spec(path=self.raw_path, obs_dict=self.obs_dict, out_prefix='e_')

        # Obtaining and applying dispersion function
        self.wcalibrate_spec(obs_dict=self.obs_dict, disp_sol_content_dict=self.disp_sol_content_dictionary,
                             prefix='e_', out_prefix='w')

    @staticmethod
    def get_obs_log(filename):

        """
        It reads a json file and return it as a python dictionary.

        Args:

            filename (str): The filename of a json file with the obslog definitions.

        Returns:

            obs_log_dict (dict): A python dictionary containing the obs log definitions.

        """

        filename = str(filename)

        with open(filename, 'r') as f:
            obs_log_dict = json.load(f)

        #except FileNotFoundError as file_not_found_message:
        #    print(file_not_found_message)

        return obs_log_dict

    def get_spec_setup(self, file):

        hdu = pyfits.open(file)
        header = hdu[0].header
        hdu.close()

        # Getting spec setup
        grism = header['GRISM'].replace('-', '_').lower()
        filter = header['FILTER'].lower()
        slit = header['APERTURE'].lower()

        return slit, grism, filter

    def check_spec_setup(self, sci_file, arc_file):

        sci_slit, sci_grism, sci_filter = self.get_spec_setup(sci_file)
        arc_slit, arc_grism, arc_filter  = self.get_spec_setup(arc_file)

        # Check Slits
        message = 'SLIT {} from SCI file {} does not match SLIT {} from ARC file {}'.format(sci_slit.upper(), sci_file,
                                                                                            arc_slit.upper(), arc_file)
        assert sci_slit == arc_slit, message

        # Check Grisms
        message = 'GRISM {} from SCI file {} does not match SLIT {} from ARC file {}'.format(sci_grism.upper(), sci_file,
                                                                                             arc_grism.upper(), arc_file)
        assert sci_grism == arc_grism, message

        # Check Filter
        if sci_grism != 'vph_all':

            message = 'FILTER {} from SCI file {} does not match FILTER {} from ARC file {}'.format(sci_filter.upper(), sci_file,
                                                                                                    arc_filter.upper(), arc_file)
            assert sci_filter == arc_filter, message

    @staticmethod
    def get_dispsol_content(fname):

        # Getting disp sol content from input disp file
        with open(fname, 'r') as f:
            content = f.read()
        f.close()

        return content

    def get_dict_of_dispsol_content(self, disp_solutions_files_dictionary, dispsol_dirname='dispsol'):

        data_dir = '{}{}'.format(self.data_dir, dispsol_dirname)

        keys = disp_solutions_files_dictionary.keys()
        values = disp_solutions_files_dictionary.values()

        content_dict = {}
        for i, key in enumerate(keys):

            content = self.get_dispsol_content(fname=data_dir + values[i])

            content_dict[key] = content

        return content_dict

    def write_modified_disp_sol_to_database(self, arc_file, dispsol_content, database_dir):


        # Getting arc_file setup
        slit, grism, filter = self.get_spec_setup(arc_file)

        # This is the file name format as given in the dispsol directory, which is based on the
        # grism - filter combination
        inp_disp_file = 'id_{}_{}.0001'.format(grism, filter)

        # Output file name with the dispersion solution (as required by IRAF)
        output_disp_file_name = 'id{}'.format(arc_file.replace('.fits', ''))

        # Correspoing txt inside the dispersion solution file
        txt_to_be_replaced = str(inp_disp_file)
        txt_new = output_disp_file_name.replace('id', '')

        modified_disp_content = dispsol_content.replace(txt_to_be_replaced, txt_new)

        # Writting new content to out_disp_file in the self.database_dir
        with open(database_dir + output_disp_file_name, 'w') as f:
            print >> f, modified_disp_content
        f.close()

    @staticmethod
    def change_dir(path):
        """
        Change to directory defined in the path variable.

        Args:
            path (str): absolute or relative path.

        """
        # Changing system to path
        if os.getcwd() is not path:
            os.chdir(path)

    @staticmethod
    def clean_path(path, ext_list):
        """
        Clean up files from path. If the path does not exist, it will be created. It's not recursive, so it
        will not clean up a sub-folder within the path.

        Args:
            path (string): path to be cleaned.
            ext_list (list): list of strings containing the extension of files to be deleted.

        """
        path = str(os.path.join(path, ''))

        if not os.path.exists(path):

            log.info('Provided path does not exist, but it was created!')
            os.makedirs(path)

        elif os.path.exists(path):
            log.info('Cleaning up ' + path + '.')
            iraf.imdelete('*tmp*', verify='no')
            for ext in ext_list:
                for _file in glob.glob(os.path.join(path, '*.' + str(ext))):
                    os.remove(_file)

    def extract_spec(self, path, obs_dict, out_prefix='e_'):
        """
        This method runs *apall* over the data set inside the path folder as defined by the obs_dict dictionary.

        Args:

            path (string): path that contaings the data to be extracted with apall.

            obs_dict (dict): dictionary containing the filename of the spectra to be extracted, the corresponding
            arc file and the sky background definition. The dictionary should be something like this:

                      self.obs_dict = {'LTTT1788_1':  {'spec': 'ccd0093c1', 'arc': 'ccd0095c1', 'sky': '-35:-20,20:35'},

                                       'LTTT1788_2':  {'spec': 'ccd0094c1', 'arc': 'ccd0095c1', 'sky': '-35:-20,20:35'},

                                       'LTTT1788_3':  {'spec': 'ccd0096c1', 'arc': 'ccd0098c1', 'sky': '-35:-20,20:35'},

                                       'LTTT1788_4':  {'spec': 'ccd0097c1', 'arc': 'ccd0098c1', 'sky': '-35:-20,20:35'}}

            out_prefix (string): letter that will be added to the filename of the data extracted.

        Returns:

            It returns 1d extracted spectra.

        """
        path = str(os.path.join(path, ''))

        # Set task parameters.
        twodspec.apextract.apall.unlearn()

        # Make sure that in "apextract" the dispaxis parameter is set to "2". If the commando below does not work,
        # do it manually.
        twodspec.apextract(dispaxis=2)

        for target, p in obs_dict.iteritems():

            spec_in = p['spec'] + '.fits'
            arc_in = p['arc'] + '.fits'

            if len(obs_dict) > 0:
                print('\n')
                log.info('Extracting star spectrum ' + p['spec'] + '...')
                apextract.apall(input=path + spec_in, output=out_prefix + spec_in, b_sample=p['sky'],
                                **self.spec_ext_flags)
                print('\n')
                log.info('Extracting arc spectrum ' + p['arc'] + '...')
                apextract.apall(input=path + arc_in, output=out_prefix + arc_in, reference=path + spec_in,
                                **self.arc_ext_flags)


    def wcalibrate_spec(self, obs_dict, disp_sol_content_dict, prefix='e_', out_prefix='w'):
        """
        This method runs a couple of IRAF routines to obtain the dispersion function for a arc files and then apply it
        to the corresponding spectra.

        Args:
            obs_dict (dict):

            prefix (str):

            out_prefix (str):

        Returns:

        """

        onedspec.identify.unlearn()
        onedspec.refspec.unlearn()
        onedspec.dispcor.unlearn()

        for target, p in obs_dict.iteritems():

            spec_in = prefix + p['spec'] + '.0001.fits'
            arc_ref = prefix + p['arc'] + '.0001.fits'

            # Checking sepctral setup
            self.check_spec_setup(spec_in, arc_ref)

            ##### Copying disp solution to 'database' dir
            # 1. Getting content dictionary with disp solutions of the corresponding arc
            slit, grism, filter = self.get_spec_setup(arc_ref)

            w1, w2 = self.wrange[grism][0], self.wrange[grism][1]

            # Getting specific disp sol content of the corresponding arc file
            key = '{}_{}'.format(grism,filter)
            content_dict = disp_sol_content_dict[key]

            # 2. Writting solution to database dir
            self.write_modified_disp_sol_to_database(arc_ref, content_dict, database_dir=self.database_dir)

            ##### Running iraf to get updated disp sol
            print('\n')
            log.info('Finding wavelength solution to reference arc ' + arc_ref + '...')
            onedspec.identify(arc_ref, **self.identify_flags)

            print('\n')
            log.info('Associating the obtained wavelength solution with the spectrum of the star:')
            log.info(spec_in + ' -----> REFSPEC = ' + arc_ref + '.')
            onedspec.refspec(spec_in, reference=arc_ref, sort='', group='')

            print('\n')
            log.info('Applying wavelength calibration to ' + spec_in + '.')
            onedspec.dispcor(spec_in, out=out_prefix + spec_in, w1=w1, w2=w2)

            if self.fit_continuum:

                onedspec.continuum.unlearn()

                print('\n')
                log.info('Fitting continuum to ' + out_prefix + spec_in + '.')

                input = out_prefix + spec_in
                output = 'cont_' + out_prefix + spec_in
                onedspec.continuum(input=input, output=output, type='fit', function='legendre', order=15, niterate=10,
                                   low_reject=2.0, high_reject=0.0)

