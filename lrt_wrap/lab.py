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


class FrozenClass(object):
    __isfrozen = False
    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError( "%r is a frozen class" % self )
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True

class _GeometryFromSpaceAndTimeAtmpy(FrozenClass):
    def __init__(self):
        self._site = atmsrf.network.stations.Table_Mountain
        self._datetime = pd.to_datetime('2022-05-26 17:52:00') 
        self._timezone = 'UTC'
        self._freeze()
    
    
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
    
    @property
    def day_of_year(self):
        return self.datetime.day_of_year
    
    # @property
    # def longitude(self):
    #     return self.site.lon
    
    # @property
    # def latitude(self):
    #     return self.site.lat
    
    @property
    def altitude(self):
        return self.site.alt/1e3
    
    def __str__(self):
        out =     ['geometry by atmpy']
        out.append('===================================')
        out.append(f'solar_zenith_angle: {self.solar_zenith_angle}')
        out.append(f'altitude: {self.altitude}')
        out.append(f'day_of_year: {self.day_of_year}')
        return '\n'.join(out)
    
        
class _GeometryFromSpaceAndTime(FrozenClass):
    def __init__(self):
        self._site = atmsrf.network.stations.Table_Mountain
        self.datetime = pd.to_datetime('2022-05-26 17:52:00') 
        # self._timezone = 'UTC'          
        self.longitude = -105.2368
        self.latitude = 40.12498
        self.altitude = 1.689
        self._freeze()
    
    def __str__(self):
        out =     ['geometry from location and datetime']
        out.append('===================================')
        out.append(f'longtude: {self.longitude}')
        out.append(f'altitude: {self.altitude}')
        out.append(f'datetime: {self.datetime}')
        return '\n'.join(out)
    

    
        
class _GeometryByHand(FrozenClass):
    def __init__(self):
        # default values are for Table mountain on '2022-05-26 17:52:00'
        self.day_of_year = 146
        self.solar_zenith_angle = 23.548980874156157
        self.altitude = 1.689
        self._freeze()
        
    def __str__(self):
        out =     ['geometry by hand']
        out.append('===================================')
        out.append(f'solar_zenith_angle: {self.solar_zenith_angle}')
        out.append(f'altitude: {self.altitude}')
        out.append(f'day_of_year: {self.day_of_year}')
        return '\n'.join(out)
    
class Geometry(FrozenClass):
    def __init__(self):
        self._by_space_time = False
        self._by_hand = False
        self._geometry = None
        self.settings = None
        self._freeze()
    
    def __str__(self):
        return self.settings.__str__()
    
    def __repr__(self):
        # print(self.__str__())
        return self.__str__()

    def from_location_and_time(self, atmpy = False):
        if atmpy:
            self._by_hand=True
            self._by_space_time = False
            out = _GeometryFromSpaceAndTimeAtmpy
        else:
            self._by_hand=False
            self._by_space_time = True
            out = _GeometryFromSpaceAndTime
        self.settings = out()    
        return 
    
    def by_hand(self):
        self._by_hand=True
        self._by_space_time = False
        self.settings =  _GeometryByHand()
        return 
    
    

class SolarSpectrum(FrozenClass):
    def __init__(self, lrt_instance, path2spectrum = 'data/solar_flux/atlas_plus_modtran'):
        self.lrt_instance = lrt_instance
        self.path2spectrum = lrt_instance.path2libradtran.joinpath(path2spectrum)
        self.wavelengthrange = (400, 800)
        self._dataset = None
        self._freeze()
    
    @property
    def dataset(self):
        if isinstance(self._dataset, type(None)):
            header = []
            with open(self.path2spectrum, 'r') as rein:
                # for i in range(5):
                maxheadlines = 20
                for i in range(maxheadlines):
                    line = rein.readline()
                    if line[0] != '#':
                        break
                    header.append(line)
                assert(i != range(maxheadlines)[-1]), 'reached and of maxheaderlines and header end is still not reached'

            header = ''.join(header)
            data = pd.read_csv(self.path2spectrum, delim_whitespace=True, skiprows=i+2, names=['wavelength', 'irradiance'], index_col=0)
            data = data.to_xarray()
            data.attrs['header'] = header
            self._dataset = data
        return self._dataset

class AtmosphereAndGases(FrozenClass):
    def __init__(self):
        self.ozone = 300.
        self.h2o =  'default'#precipitable water in mm; default is close to 14.29475 (no idea where it comes from)
        self.parameterization = 'reptran'
        self.resolution = 'coarse'
        self.scattering_off = False
        self._freeze()

class Aerosols(FrozenClass):
    def __init__(self):
        self.scene = 'aerosol_default' #### TODO add vulcano, etc Manual p. 48
        self.ssa = None
        self.gg = None
        self.aod = None
        self._freeze()

class LibRadTrans(object):
    def __init__(self, aerosols = True):
        self.path2libradtran = pl.Path('/home/hagen/programms/libRadtran-2.0.4')
        self._p2lrtbin = self.path2libradtran.joinpath('bin/')
        # self._p2lrttmp = self.path2libradtran.joinpath('tmp/')
        # f2input = 'UVSPEC_CLEAR.INP'
        # list(p2lrttmp.glob('*'))
        self.solver = 'disort'
        self.solar_spectrum = SolarSpectrum(self,path2spectrum = 'data/solar_flux/atlas_plus_modtran')
        self.output_quantity = 'radiances' #brightness, reflectivity, transmittance
        self.atmosphere = AtmosphereAndGases() 
        self.geometry = Geometry()
        self.geometry.from_location_and_time(atmpy=True)
        self.test_input = None # This is a str that will be added to the input for testing
        
        if aerosols:
            self.aerosols = Aerosols()
        else: 
            self.aerosols = None
    
    @property
    def input_string(self):

        input_list = []                        
        ##### Location of atmospheric profile file. 
        input_list.append('atmosphere_file ../data/atmmod/afglus.dat')

        #### Solar spectrum and limits
        # Location of the extraterrestrial spectrum
        if self.output_quantity != 'transmittance':
            input_list.append(f'source solar {self.solar_spectrum.path2spectrum.as_posix()}')
            
        input_list.append(f'wavelength {self.solar_spectrum.wavelengthrange[0]} {self.solar_spectrum.wavelengthrange[1]}')   # Wavelength range [nm]

        #### atmosphere and gasses
        input_list.append(f'mol_abs_param {self.atmosphere.parameterization} {self.atmosphere.resolution}')
        if self.atmosphere.scattering_off:
            input_list.append(f'no_scattering mol')
        
        input_list.append(f'mol_modify O3 {self.atmosphere.ozone:0.6f} DU')    # Set ozone column
        if self.atmosphere.h2o != 'default':
            input_list.append(f'mol_modify H2O {self.atmosphere.h2o:0.6f} MM')    # 
        
    
        
        #### geometry
        if self.geometry._by_hand:
            input_list.append(f'day_of_year {self.geometry.settings.day_of_year}')          # Correct for Earth-Sun distance
            input_list.append(f'sza {self.geometry.settings.solar_zenith_angle:0.6f}')                 # Solar zenith angle
        else:
            # assert(False), 'noep'
            dt = self.geometry.settings.datetime
            time = f'{dt.year:4d} {dt.month:02d} {dt.day:02d} {dt.hour:02d} {dt.minute:02d} {dt.second:02d}'
            input_list.append(f'time {time}')  
            
            lat = self.geometry.settings.latitude 
            if lat < 0:
                hem='S'
                lat *=-1
            else:
                hem='N'
            input_list.append(f'latitude {hem} {lat:0.6f}')
            
            lon = self.geometry.settings.longitude 
            if lon < 0:
                hem='W'
                lon *=-1
            else:
                hem='E'
            input_list.append(f'longitude {hem} {lon:0.6f}')
            
        if not isinstance(self.geometry.settings.altitude, type(None)):
            input_list.append(f'altitude {self.geometry.settings.altitude:0.6f}')
        
        #### surface
        input_list.append('albedo 0.2')               # Surface albedo
        
        #### aerosols
        if not isinstance(self.aerosols, type(None)):
            if self.aerosols.scene == 'aerosol_default':
                input_list.append('aerosol_default')
            else:
                assert(False), 'other scenes need programming'
                
            if not isinstance(self.aerosols.ssa, type(None)):
                input_list.append(f'aerosol_modify ssa set {self.aerosols.ssa:0.6f}')
            if not isinstance(self.aerosols.gg, type(None)):
                input_list.append(f'aerosol_modify gg set {self.aerosols.gg:0.6f}')
            if not isinstance(self.aerosols.aod, type(None)):
                input_list.append(f'aerosol_modify tau set {self.aerosols.aod:0.6f}')
        
        #### model settings
        input_list.append(f'rte_solver {self.solver}')        # Radiative transfer equation solver
        input_list.append('number_of_streams  6')     # Number of streams
        
        if self.output_quantity != 'radiances':
            input_list.append(f'output_quantity {self.output_quantity}')
        
        input_list.append('quiet')
        
        if not isinstance(self.test_input,type(None)):
            input_list.append(self.test_input)
        
        # the original
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
        return LibRadTransResult(self,self.run_stdout)
            
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
        if self.lrt_instance.geometry._by_space_time:
            ds['datetime'] = self.lrt_instance.geometry.settings.datetime
        #### TODO below will cause the zenith angle calculation to be executed a second time! consider dealing with this in the property getter.
        elif self.lrt_instance.geometry._by_hand:
            ds['solar_zenith_angle'] = self.lrt_instance.geometry.settings.solar_zenith_angle 
        return ds