#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 12:24:45 2022

@author: hagen
"""
import pandas as pd
import pathlib as pl
import subprocess
import io
import atmPy.general.measurement_site as atmms
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.surfrad.surfrad as atmsrf
import numpy as np

class SpaceAndTime(object):
    def __init__(self):
        self._site = atmsrf.network.stations.Table_Mountain
        self._datetime = pd.to_datetime('2022-05-26 17:52:00') 
        self._timezone = 'UTC'
    
    @property
    def solar_zenith_angle(self):
        sza = 90 - np.rad2deg(self.site.get_sun_position(self.datetime)['altitude'])
        return sza
    
    @property
    def timezone(self):
        return self._timezone
    
    @timezone.setter
    def timezone(self, value):
        if value != 'UTC':
            raise ValueError("NO! UTC!! Nothing else! Don't even try!")
        return
    
    @property
    def site(self):
        return self._site 
    
    @site.setter
    def site(self, value):
        required = ['lat', 'lon']
        if isinstance(value, atmms.Station):
            site = value
        elif isinstance(value, dict):
            for req in required:
                assert(req in value.keys()), f'The provided dictionary is missing the key: "{req}".'
            if 'name' not in value.keys(): # this was necessary because atmPy wanted it. should probably be dealt with in atmPy
                value['name'] = 'station'
            site = atmms.Station(lat = value['lat'], lon = value['lon'], name = value['name'])
        else:
            ValueError('Dictionary with lat and lon keys or atmPy.general.measurement_site.Station instance are required.')
        self._site = site
        
    @property 
    def datetime(self):
        return self._datetime
    
    @datetime.setter 
    def datetime(self, value):
        if isinstance(value,pd._libs.tslibs.timestamps.Timestamp):
            dt = value
        elif isinstance(value, str):
            dt = pd.to_datetime(value)
        else:
            assert(False), 'pandas Timestamp or string required'
        self._datetime = dt
        return
                
            

class SolarSpectrum(object):
    def __init__(self, lrt_instance, path2spectrum = 'data/solar_flux/atlas_plus_modtran'):
        self.lrt_instance = lrt_instance
        self.path2spectrum = lrt_instance.path2libradtran.joinpath(path2spectrum)
        self._data = None
        self._header = None
    
    @property
    def data(self):
        if isinstance(self._data, type(None)):
            # self._data = pd.read_csv(lrtrun.solar_spectrum.path2spectrum, delim_whitespace=True, skiprows=6, names=['wavelength', 'irradiance'], index_col=0)
            self.header
        return self._data
    
    @property
    def header(self):
        if isinstance(self._header, type(None)):
            header = []
            with open(self.solar_spectrum.path2spectrum, 'r') as rein:
                # for i in range(5):
                maxheadlines = 20
                for i in range(maxheadlines):
                    line = rein.readline()
                    if line[0] != '#':
                        break
                    header.append(line)
                assert(i != range(maxheadlines)[-1]), 'reached and of maxheaderlines and header end is still not reached'

            self._header = ''.join(header)
            self._data = pd.read_csv(self.solar_spectrum.path2spectrum, delim_whitespace=True, skiprows=i+2, names=['wavelength', 'irradiance'], index_col=0)
        return self._header

class AtmosphereAndGases(object):
    def __init__(self):
        self.ozone = 300.
        self.h2o =  'default'#precipitable water in mm; default is close to 14.29475 (no idea where it comes from)

class LibRadTrans(object):
    def __init__(self):
        self.path2libradtran = pl.Path('/home/hagen/programms/libRadtran-2.0.4')
        self._p2lrtbin = self.path2libradtran.joinpath('bin/')
        # self._p2lrttmp = self.path2libradtran.joinpath('tmp/')
        # f2input = 'UVSPEC_CLEAR.INP'
        # list(p2lrttmp.glob('*'))
        self.solver = 'disort'
        self.wavelengthrange = (400, 800)
        self.solar_spectrum = SolarSpectrum(self,path2spectrum = 'data/solar_flux/atlas_plus_modtran')
        self.output_quantity = 'radiances' #brightness, reflectivity, transmittance
        self.atmosphere = AtmosphereAndGases() 
        self.spacetime = SpaceAndTime()
    
    @property
    def input_string(self):
#         input_str = f"""                         
#         # Location of atmospheric profile file. 
#         atmosphere_file ../data/atmmod/afglus.dat

#         # Location of the extraterrestrial spectrum
#         source solar {self.solar_spectrum.path2spectrum.as_posix()}

#         mol_modify O3 300. DU    # Set ozone column
#         day_of_year 170          # Correct for Earth-Sun distance
#         albedo 0.2               # Surface albedo
#         sza 32.0                 # Solar zenith angle
#         rte_solver {self.solver}        # Radiative transfer equation solver
#         number_of_streams  6     # Number of streams
#         wavelength {self.wavelengthrange[0]} {self.wavelengthrange[1]}   # Wavelength range [nm]
        
#         output_quantity transmittance
        
#         quiet
#         """
        input_list = []                        
        # Location of atmospheric profile file. 
        input_list.append('atmosphere_file ../data/atmmod/afglus.dat')

        # Location of the extraterrestrial spectrum
        input_list.append(f'source solar {self.solar_spectrum.path2spectrum.as_posix()}')

        input_list.append(f'mol_modify O3 {self.atmosphere.ozone:0.1f} DU')    # Set ozone column
        if self.atmosphere.h2o != 'default':
            input_list.append(f'mol_modify H2O {self.atmosphere.h2o:0.4f} MM')    # Set ozone column
        input_list.append(f'day_of_year {self.spacetime.datetime.day_of_year}')          # Correct for Earth-Sun distance
        input_list.append('albedo 0.2')               # Surface albedo
        input_list.append(f'sza {self.spacetime.solar_zenith_angle:0.4f}')                 # Solar zenith angle
        input_list.append(f'rte_solver {self.solver}')        # Radiative transfer equation solver
        input_list.append('number_of_streams  6')     # Number of streams
        input_list.append(f'wavelength {self.wavelengthrange[0]} {self.wavelengthrange[1]}')   # Wavelength range [nm]
        
        if self.output_quantity != 'radiances':
            input_list.append(f'output_quantity {self.output_quantity}')
        
        input_list.append('quiet')
        # """
        return '\n'.join(input_list)
    
    @property
    def command4execution(self):
        command = f"cd {self._p2lrtbin.as_posix()}; ../bin/uvspec"
        return command
    
    def run_model(self):
        bla = subprocess.Popen(self.command4execution, shell = True,  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.run_stdout, self.run_error = bla.communicate(input = self.input_string.encode())
        if len(self.run_error)>0:
            raise RuntimeError(self.run_error)
        self.result = LibRadTransResult(self,self.run_stdout)
            
class LibRadTransResult(object):
    def __init__(self, lrt_instance, stdout):
        # self.stdout = stdout
        self.lrt_instance = lrt_instance
        self.abbriviation_dict =    {'lambda':  'wavelength',  
                                     'edir':    'direct_beam_irradiance',
                                     'tdir':    'direct_beam_transmittance',
                                     'edn':     'diffuse_down_irradiance', 
                                     'tdn':     'diffuse_down_transmittance', 
                                     'eup':     'diffuse_up_irradiance', 
                                     'tup':     'diffuse_up_transmittance', 
                                     'uavgdir': 'direct_beam_contribution_to_mean_intensity', 
                                     'uavgdn':  'diffuse_down_radiation_contribution_to_mean_intensity',
                                     'uavgup':  'diffuse_up_radiation_contribution_to_mean_intensity'
                                    }
        self.dataset = self.stdout2xarray(stdout)
        
    def stdout2xarray(self,stdout):
        if self.lrt_instance.output_quantity == 'radiances':
            params = ['lambda',  'edir', 'edn', 'eup', 'uavgdir', 'uavgdn', 'uavgup']
        elif self.lrt_instance.output_quantity == 'transmittance':
            params = ['lambda',  'tdir', 'tdn', 'tup', 'uavgdir', 'uavgdn', 'uavgup']
        df = pd.read_csv(io.StringIO(stdout.decode()), names = [self.abbriviation_dict[p] for p in params], index_col=0,
                    delim_whitespace=True,)
        ds = df.to_xarray()
        ds['datetime'] = self.lrt_instance.spacetime.datetime
        #### TODO below will cause the zenith angle calculation to be executed a second time! consider dealing with this in the property getter.
        ds['solar_zenith_angle'] = self.lrt_instance.spacetime.solar_zenith_angle 
        return ds