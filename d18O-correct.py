#! /usr/bin/env python

"""
d18O-correct.py
(c) 2016-2019 Andrew D. Wickert
Licensed under GNU GPL v3
"""

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d 

class meltwaterCore(object):

  def __init__(self):
    pass

  def import_sheet(self, wb='MississippiAll.xlsx'):
    self.sheet = pd.read_excel(wb, sheet_name=self.sheet_name)
    self.sheet = self.sheet.drop(0) # drop row with units data
    self.z = self.sheet['Depth'].values.astype(float) # cm
    self.MgCa = self.sheet['Mg/Ca'].values.astype(float)
    self.d18O_VPDB = self.sheet[self.sheet.keys()[2]].values.astype(float) # different names

  def VPDB_to_VSMOW(self, d18O_VPDB):
    """
    We are following Hut (1987):
    del-18Owater VPDB = del-18Owater VSMOW - 0.27 (Hut, 1987),
    after Flower (2004)
    d18O_VSMOW = d18O_VPDB + 0.27
    ; why is this so different from the Coplen et al., 1983:
    d18O_VSMOW = 1.03091 * d18O_VPDB + 30.91
    ?
    Looks like the one we use is a conversion for CO2 equivalent...
    """
    d18O_VSMOW = d18O_VPDB + 0.27
    return d18O_VSMOW
    
  def MgCa_T(self):
    """
    Species may be 'G. ruber', 'G. bulloides', or 'N. pachyderma'
    
    Returns low error, median, high error. 1-sigma errors
    
    Mg/Ca in mmol/mol (so is Mg/Ca * 1000)

    G. Generic planktonic 350-500 um:
    Following Anand et al. (2003)
    with error bars from Khider et al. (2015)
    but not using any of their data on correlations with
    salinity and dissolution
        Anand's original equation for 350-500 um, used by Flower's group:
        Mg/Ca = (0.38 +/- 0.02) * exp(0.090 * T)
        with T in degrees C, and the error the original from Anand et al. (2003)

    G. RUBER white 350-500 um:
    Following Anand et al. (2003)
    with error bars from Khider et al. (2015)
    but not using any of their data on correlations with
    salinity and dissolution
        with T in degrees C, and the error the original from Anand et al. (2003)
        NOT YET ADDED -- FOLLOWING FLOWER IN USING GENERIC EQUATION

    G. RUBER pink 350-500 um:
    Following Anand et al. (2003)
    with error bars from Khider et al. (2015)
        but not using any of their data on correlations with
        salinity and dissolution
        with T in degrees C, and the error the original from Anand et al. (2003)
        NOT YET ADDED -- FOLLOWING FLOWER IN USING GENERIC EQUATION
        
    G. BULLOIDES:
      Mg/Ca = (1.006 +/- 0.032) * exp((0.065 +/- 0.003)*T) (R^2 = 0.82).
      Vasquez Riveiros et al., 2016
      
    N. PACHYDERMA:
      Mg/Ca = (0.580 +/- 0.016) * exp((0.084 +/- 0.006)*T) (R^2 = 0.70).
      Vasquez Riveiros et al., 2016
    
    O. UNIVERSA
      ln(Mg/Ca) = (-0.162 +/- 0.295) + 0.096 (+/- 0.014) * T
      Russel et al., 2004

    """
    if (self.foram_species == 'Generic planktonic 350-500 um') \
      or (self.foram_species[:8] == 'G. ruber'):
      # NOW FOR RUBER, LOCAL CALIB, RICHEY CORE TOP
      #T = (1/0.090) * np.log(self.MgCa/0.38)
      T = (1/0.090) * np.log(self.MgCa/0.449)
      Terr = 0.45
      Trange = np.array([T-Terr, T, T+Terr]).transpose()
    # ! CORRECTLY PROPAGATE ERROR WITH NUMERICAL TEST
    #if self.foram_species == 'G. ruber white 350-500 um':
    #  pass
    # ! CORRECTLY PROPAGATE ERROR WITH NUMERICAL TEST
    #if self.foram_species == 'G. ruber pink 350-500 um':
    #  pass
    # ! CORRECTLY PROPAGATE ERROR WITH NUMERICAL TEST
    elif self.foram_species == 'G. bulloides':
      T = (1/0.065) * np.log(self.MgCa/1.006)
      Terr = 1.0 # JUST A QUICK ESTIMATE!!!!!!!!!!!!!!!!!!!!!!
      Trange = np.array([T-Terr, T, T+Terr]).transpose()
    elif self.foram_species == 'N. pachyderma':
      T = (1/0.084) * np.log(self.MgCa/0.580)
      Terr = 1.1 # JUST A QUICK ESTIMATE!!!!!!!!!!!!!!!!!!!!!!
      Trange = np.array([T-Terr, T, T+Terr]).transpose()
    elif self.foram_species == 'O. universa':
      T = (1/0.096) * np.log(self.MgCa/0.85)
      Terr = 1.0 # JUST A QUICK ESTIMATE!!!!!!!!!!!!!!!!!!!!!!
      Trange = np.array([T-Terr, T, T+Terr]).transpose()
    else:
      sys.exit("WRONG SPECIES!")
    
    self.T = Trange
    
  def MgCa_T_composite(self, _sheets, age_models, ma_half_window=500, species='Generic planktonic 350-500 um'):
    """
    '_sheets' must be a list of sheet names
    """
    MgCa_T_all = []
    MgCa_ages_all = []
    for i in range(len(_sheets)):#3
      self.sheet_name = _sheets[i]
      self.age_model_filename = age_models[i]
      self.foram_species = species # Allow us to have variables species all together in the future
      self.import_sheet()
      self.apply_age_model()
      self.MgCa_T()
      MgCa_T_all.append(self.T.copy())
      MgCa_ages_all.append(self.ages.copy())
    MgCa_T_all = np.vstack(MgCa_T_all)
    MgCa_ages_all = np.hstack(MgCa_ages_all)

    # Moving average
    self.MgCa_ages_ma = np.arange(ma_half_window/10., \
                            np.round(np.nanmax(MgCa_ages_all) + ma_half_window), \
                            ma_half_window)
    MgCa_T_ma = []
    MgCa_T_STD_ma = [] # NOTE: not including conversion equation error -- just considering variability in data to be most important
    for age in self.MgCa_ages_ma:
      age_window = (MgCa_ages_all > (age - ma_half_window)) * \
                   (MgCa_ages_all < (age + ma_half_window))
      MgCa_T_ma.append(np.nanmean(MgCa_T_all[:,1][age_window]))
      MgCa_T_STD_ma.append(np.nanstd(MgCa_T_all[:,1][age_window]))
    self.MgCa_T_ma = np.array(MgCa_T_ma)
    self.MgCa_T_STD_ma = np.array(MgCa_T_STD_ma)
    self.MgCa_ages_ma = np.array(self.MgCa_ages_ma[np.isnan(self.MgCa_T_ma) == False])
    self.MgCa_T_ma = np.array(self.MgCa_T_ma[np.isnan(self.MgCa_T_ma) == False])
    self.MgCa_T_STD_ma = np.array(self.MgCa_T_STD_ma[np.isnan(self.MgCa_T_STD_ma) == False])

  def d18O_T_correction(self, d18O_uncorrected, T):
    """
    Bemis et al. (1998) for O. universa (high-light)
        shown by Thunell et al. (1999) to also apply to G. ruber
        
    T = 14.9 - 4.8 (d18O_uncorrected - d18O_corrected)
    # for later -- see p. 1172 of Taylor et al. (2015)

    d18O_corrected is given with respect to the same standard as 
    d18O_uncorrected
    """
    # !!!!!!!!!!!! I DON'T INCLUDE ERROR ESTIMATES HERE!!!!!
    if self.foram_species[:8] == 'G. ruber':
      d18O_T_corr = ((T - 14.9) / 4.8) + d18O_uncorrected
    # SHOULD HAVE A SQUARED TERM HERE -- SEE TAYLOR ET AL!
    # DROPPING IT FOR THIS FIRST PASS TO MAKE THE EQUATION SOLVABLE WITHOUT
    # A NUMERICAL SOLVER.
    elif self.foram_species == 'G. bulloides':
      d18O_T_corr = ((T - 16.9) / 4.38) + d18O_uncorrected
    elif self.foram_species == 'N. pachyderma':
      d18O_T_corr = ((T - 16.5) / 4.8) + d18O_uncorrected
    elif self.foram_species == 'O. universa':
      d18O_T_corr = ((T - 16.5) / 4.8) + d18O_uncorrected
    else:
      sys.exit("WRONG SPECIES!")
    return d18O_T_corr
    # LA

  # I should include the SL curve as an object
  # Also, I need to create a way from the outside to determine the SL curve
  #def ice_volume_correction(self, ages_d18O, SL_curve_file='Lambeck_SL_curve.txt', header_length=9):
  def ice_volume_correction(self, ages_d18O, SL_curve_file='Spratt_SL_curve.txt', header_length=1):
    """
    Schrag et al. (2002)
    d18O sw + 1.0 +/- 0.1 at LGM (Schrag et al., 
    I take this to nominally be 26-20 ka, within the full 26.5-(20-19) ka
      bounds suggested by Clark et al. (2009).
    From this, I rescale the sea-level curve, here from Lambeck et al. (2014),
      to 0-1, and add on the +/- 0.1 error, also rescaled.
    """
    # !!!!!!!!!!!!!!!! CURRENTLY JUST WORKING ON MEAN AGES, NOT RANGE OF AGES.
    # IS THIS A SOURCE OF UNCERTAINTY? HOW TO PROPAGATE?
    
    print "Using sea-level curve", SL_curve_file
    
    # Difference and error
    LGM_d18O_diff = 1
    LGM_d18O_diff_error = 0.1
    # Sea-level -- add point for 2020 (approx. today) at start
    SL = np.genfromtxt(SL_curve_file, skip_header=header_length)
    ages_ka = np.hstack((-70, SL[:,0]))
    esl = np.hstack((0, SL[:,2]))
    esl_err = np.hstack((0, SL[:,3]))
    # Interpolate to 100-year chunks
    _ages = ages_ka * 1000
    f = interp1d(_ages, esl, bounds_error=False)
    ages_i = 100 * np.arange( np.ceil(_ages[0]/100.), \
                              np.floor(_ages[-1]/100. ))
    esl_i = f(ages_i)
    # LGM ESL
    LGM_esl = np.mean(esl_i[(ages_i >= 20000) * (ages_i <= 26000)])
    # Correct
    SL_at_ages_d18O = f(ages_d18O)

    _mean = LGM_d18O_diff * SL_at_ages_d18O / LGM_esl
    _min = _mean - (LGM_d18O_diff_error * SL_at_ages_d18O / LGM_esl)
    _max = _mean + (LGM_d18O_diff_error * SL_at_ages_d18O / LGM_esl)
    
    output = np.vstack((_min, _mean, _max)).transpose()
    
    output *= -1 # correction, so flip sign
    
    return output
    

  def fill_missing_T(self): # _interp
    
    f = interp1d(self.MgCa_ages_ma, self.MgCa_T_ma, bounds_error=False)
    fSD = interp1d(self.MgCa_ages_ma, self.MgCa_T_STD_ma, bounds_error=False)
    # Age model for values -- CRUDE
    T_filled = self.T.copy()
    T_filled[:,1] = f(self.ages)
    SD = fSD(self.ages)
    T_filled[:,0] = T_filled[:,1] - 3*SD
    T_filled[:,2] = T_filled[:,1] + 3*SD
    # Replace values where needed
    self.T_filled = self.T.copy()
    self.T_filled[np.isnan(self.T)] = T_filled[np.isnan(self.T)]
    
  def fill_missing_T_window(self, half_window):
    T_filled = self.T.copy()
    for i in range(len(self.T)):
      Ti = self.T[i,1]
      if np.isnan(Ti):
        age_noT = self.ages[i]
        if np.isnan(age_noT):
          pass
        else:
          age_window = (self.ages > (age_noT - half_window)) * \
                       (self.ages < (age_noT + half_window))
          if age_window.any():
            T_local = np.nanmean(self.T[:,1][age_window])
            T_min = np.nanmean(self.T[:,0][age_window])
            T_max = np.nanmean(self.T[:,2][age_window])
            T_std_local = np.nanstd(self.T[:,1][age_window])
            T_filled[i] = [T_min - T_std_local, T_local, T_max + T_std_local]
    self.T_filled = T_filled

  def apply_age_model(self):
    """
    Age model (run outside with Bacon)
    """
    age_model = np.genfromtxt(self.age_model_filename, delimiter='\t', skip_header=1)
    depths_i = age_model[:,0]
    ages_median_i = age_model[:,4] # use weighted median
    ages_min_i = age_model[:,1]
    ages_max_i = age_model[:,2]
    self.fage = interp1d(depths_i, ages_median_i, bounds_error=False)
    self.fage_min = interp1d(depths_i, ages_min_i, bounds_error=False)
    self.fage_max = interp1d(depths_i, ages_max_i, bounds_error=False)
    # Age model for values
    self.ages = self.fage(self.z)

  def correct_d18O(self):
    # CORRECT OXYGEN ISOTOPE VALUES
    # PDB to SMOW
    d18O_VSMOW = self.VPDB_to_VSMOW(self.d18O_VPDB)
    # T correction -- COULD BE DONE BETTER (NO ANALYTICAL ERROR INCLUDED)
    # ANALYTICAL ERROR IS TYPICALLY VERY SMALL, BUT.... !!!!!!!!!!!!!!!!!!
    d18O_VSMOW_Tcorr_min = self.d18O_T_correction(d18O_VSMOW, self.T_filled[:,0])
    d18O_VSMOW_Tcorr_mean = self.d18O_T_correction(d18O_VSMOW, self.T_filled[:,1])
    d18O_VSMOW_Tcorr_max = self.d18O_T_correction(d18O_VSMOW, self.T_filled[:,2])
    d18O_VSMOW_Tcorr = np.vstack((d18O_VSMOW_Tcorr_min, \
                                  d18O_VSMOW_Tcorr_mean,
                                  d18O_VSMOW_Tcorr_max,)).transpose()
    # Ice volume correction
    SLcorr = self.ice_volume_correction(self.ages)

    # Final
    self.d18O_sw_ivc = d18O_VSMOW_Tcorr + SLcorr

  def d18O_moving_average(self, ma_half_window):
    # Moving average
    self.ages_ma = np.arange(ma_half_window/10., \
                            np.round(np.nanmax(self.ages) + ma_half_window), \
                            ma_half_window)
    d18O_sw_ivc_ma = []
    for age in self.ages_ma:
      age_window = (self.ages > (age - ma_half_window)) * \
                   (self.ages < (age + ma_half_window))
      d18O_sw_ivc_ma.append(np.nanmean(self.d18O_sw_ivc[age_window]))
    self.d18O_sw_ivc_ma = np.array(d18O_sw_ivc_ma)

  def error_for_plotting(self, sigma=2):
    """
    2 sigma error by default
    """
    self.err_minus = self.d18O_sw_ivc[:,1] - self.d18O_sw_ivc[:,0]
    self.err_plus = self.d18O_sw_ivc[:,2] - self.d18O_sw_ivc[:,1]

    self.d18O_err = np.vstack((self.err_minus, self.err_plus))
    self.d18O_err *= sigma
    self.ages_err  = np.vstack((self.ages - self.fage_min(self.z), \
                               self.fage_max(self.z) - self.ages))

  def output(self):
    output = {}
    output['ages'] = self.ages.copy()
    output['d18O_sw_ivc'] = self.d18O_sw_ivc.copy()
    output['ages_ma'] = self.ages_ma.copy()
    output['d18O_sw_ivc_ma'] = self.d18O_sw_ivc_ma.copy()
    output['ages_err'] = self.ages_err.copy()
    output['d18O_err'] = self.d18O_err.copy()
    return output
