#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 15:28:01 2022

@author: hagen
"""

import atmPy.aerosols.size_distribution.sizedistribution as atmsd
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.surfrad.surfrad as atmsrf
import atmPy.general.measurement_site as atmms
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.cal_facility.lab as atmcal

import pathlib as pl
import io
import pandas as pd
import matplotlib.pyplot as plt


from importlib import reload

import xarray as xr
import pathlib as pl
import numpy as np
from scipy import integrate
import lrt_wrap.lab as lrtw


class Mfrsr(object):
    #### TODO redo with MFR ... to simulate an albedo reck
    def __init__(self):
        self.lrt_init = xr.Dataset()
        self._mfrsr_calibration = None
        # self._aerosol_modes = None
        # self._aerosol_optical_properties = None
        self.mfrsr_calibration = '/home/hagen/Dropbox/projecte/libradtran_pg/V026217216C.SPR.txt'
        self.aerosol_modes = {'center':  [100, 2000],
                              'width': [0.2, 0.2],
                              'particle_no':  [1e7, 1e3],
                              'n':    [1.5 + 0.001j, 1.5 + 0.001j],
                              }

        return

    @property
    def mfrsr_calibration(self):
        return self._mfrsr_calibration

    @mfrsr_calibration.setter
    def mfrsr_calibration(self, value):
        # read filter responds
        # TODO get here from date and installation files
        path2mfrsr_cal = pl.Path(value)
        ds = atmcal.read_mfrsr_cal(path2mfrsr_cal)
        ds.attrs['path2calibration'] = value
        self.lrt_init['filter_responds'] =  ds.responds
        self.lrt_init['filter_statistics'] =  ds.statistics
        self._mfrsr_calibration = ds

    @property
    def aerosol_modes(self):
        return self.lrt_init['aerosol_modes'] 

    @aerosol_modes.setter
    def aerosol_modes(self, aerosol_modes):
        def row2sizedist(row):
            diameter_range = [50, 20e4]
            numberOfDiameters = 100
            sd = atmsd.simulate_sizedistribution(diameter=diameter_range,
                                                 numberOfDiameters=numberOfDiameters,
                                                 centerOfAerosolMode=row.center.real,
                                                 widthOfAerosolMode=row.width.real,
                                                 numberOfParticsInMode=row.particle_no.real,
                                                 )
            return sd

        aerosol_modes = pd.DataFrame(aerosol_modes)

        aerosol_modes.index.name = 'aerosol_mode_idx'
        aerosol_modes.columns.name = 'aerosol_mode_param'

        # create the sizedistribution instancec and add to the mode list
        aerosol_modes['sizedist'] = aerosol_modes.apply(row2sizedist, axis=1)

        # add aerosol models to the dataset
        self.lrt_init['aerosol_modes'] = aerosol_modes
        
    @property 
    def aerosol_optical_properties(self):
        # if isinstance(self._aerosol_optical_properties, type(None)):
        if 'aerosol_opt_prop' not in self.lrt_init.variables.keys():
            df_aeros_input = pd.DataFrame()
            for channel in self.lrt_init.channel:
                if self.lrt_init.filter_responds.sel(channel = channel).dropna('wavelength').shape[0] == 0:
                    print(f'channel {channel.item()} has no valid data')
                    continue
            
                filter_central_wl = self.lrt_init.filter_statistics.sel(channel = channel, stats = 'CENT').item()
            
                df_aeros_opt_prop = pd.DataFrame()
                df_phase_func = pd.DataFrame()
                for modeid in  self.lrt_init.aerosol_mode_idx:
                    # break
            
                    aeromode = self.lrt_init.aerosol_modes.sel(aerosol_mode_idx = modeid)
                    sd = aeromode.sel(aerosol_mode_param = 'sizedist').item().copy()
            
                    # sd = ds_lrt_init.aerosol_modes.sel(aerosol_mode_idx = 0, aerosol_mode_param = 'sizedist').values.item().copy()
            
                    sd.optical_properties.parameters.wavelength = filter_central_wl
                    sd.optical_properties.parameters.refractive_index = aeromode.sel(aerosol_mode_param = 'n').item()
            
                    df_phase_func[modeid.item()] = sd.optical_properties.angular_scatt_func.iloc[0]
            
                    # sd.optical_properties.scattering_coeff + 
                    sd.optical_properties.extinction_coeff.values.item()
                    sd.optical_properties.extinction_coeff.iloc[0,0]
            
                    df_aeros_opt_prop[modeid.item()] = pd.Series({'AOD_scatt': sd.optical_properties.scattering_coeff.iloc[0,0], 'AOD_ext':sd.optical_properties.extinction_coeff.iloc[0,0]})
            
                asf = df_phase_func.sum(axis = 1)
                gg = integrate.simps(asf.values * np.cos(asf.index), asf.index) / integrate.simps(asf.values, asf.index)
                ssa = df_aeros_opt_prop.sum(axis = 1).AOD_scatt / df_aeros_opt_prop.sum(axis = 1).AOD_ext
                aod = df_aeros_opt_prop.sum(axis = 1).AOD_ext
                df_aeros_input[channel.item()] = pd.Series({'gg':gg, 'ssa':ssa, 'aod':aod})
                
            df_aeros_input.index.name = 'aerosol_opt_param'
            df_aeros_input.columns.name = 'channel'
            self.lrt_init['aerosol_opt_prop'] = df_aeros_input
            # self._aerosol_optical_properties = True
            
        return self.lrt_init['aerosol_opt_prop']
    
    def simulate_measurement(self):
        # convolved_list = []
        # convolved_int_list = []
        weighted_avg_list = []
        # keep_spectral = False
        # the following can be done more simple if the entrire spectral range is considered at once, however, it will take longer to calculate... still something to consider if broad band is considered?!?
        for channel in self.lrt_init.channel:
            # break
            channel = int(channel)
            filter_responds = self.lrt_init.sel(channel = channel)['filter_responds']
            filter_responds = filter_responds.dropna('wavelength')
            if filter_responds.shape[0] == 0:
                print(f'channel {channel} has no valid data')
                continue
            wllims = (float(filter_responds.wavelength.min()) , float(filter_responds.wavelength.max()))
            if 0:
                if max(wllims) > 800:
                    #### TODO take care of the 870, 940 and 1625 channel
                    print(f'channel {channel} has wl above 800 nm, deal with it later')
                    continue
            # set up the lrtrun
            lrtrun = lrtw.LibRadTrans()
            lrtrun.output_quantity = 'transmittance'
            # lrtrun.spacetime.datetime = dt
            lrtrun.solar_spectrum.wavelengthrange = wllims
            
            #### TODO atmosphere
            #### TODO Ozone, how? or fit it?
            #### TODO water, how? or fit it?
            
            
            # aerosols
            lrtrun.aerosols.gg = self.lrt_init.aerosol_opt_prop.sel(channel = channel, aerosol_opt_param = 'gg').item()
            lrtrun.aerosols.ssa = self.lrt_init.aerosol_opt_prop.sel(channel = channel, aerosol_opt_param = 'ssa').item()
            lrtrun.aerosols.aod = self.lrt_init.aerosol_opt_prop.sel(channel = channel, aerosol_opt_param = 'aod').item()
        
        
            # run model
            lrtres = lrtrun.run_model()
            ds_rtm = lrtres.dataset
        
            # interpolate filter function on irradiance specturm
            filter_responds = filter_responds.interp({'wavelength': ds_rtm.wavelength})
            weighted_avg = (filter_responds * ds_rtm[['direct_beam_transmittance','diffuse_down_transmittance']]).sum()/filter_responds.sum()
            weighted_avg_list.append(weighted_avg)
            # # convolve filter function with irradiances
            # da_convolv_resp = filter_responds * ds_rtm[['direct_beam_transmittance','diffuse_down_transmittance']]
            
            
            # convolved_int_list.append(da_convolv_resp.integrate('wavelength'))
            # if keep_spectral:
            #     convolved_list.append(da_convolv_resp)
            # break
        
        mfrsr_trans = xr.concat(weighted_avg_list, 'channel')
        return mfrsr_trans