%% CoastalTools classes and functions
% CoastalTools uses the <matlab:doc('muitoolbox') muitoolbox> and 
% <matlab:doc('dstoolbox') dstoolbox>. In addition it includes the
% following classes and functions.

%% CoastalTools Classes
% * *CoastalTools* – main UI
% * *CT_BeachAnalysis* – tools to analyse beach profile data
% * *CT_ModelPlots* – UI for CoastalToools plots (modifies ModelPlot class)
% * *CT_PlotFig* – extension of PlotFig class to plot profile data
% * *CT_SimUI* – UI for shoreline and profile simulation models
% * *CT_WaveModels* – models using wind, wave and water level data
% * *Sim_BMV* - runs the BMV beach profile simulation model
% * *Sim_BMVinput* – data input for the BMV model
% * *Sim_BMVmodel* – setup and running of BMV simulation
% * *Sim_YGOR* - runs the YGOR shoreline simulation model
% * *Sim_YGORinput* – data input for the YGOR model

%% CoastalTools Functions
% * _data_cleanup_ – set of functions to help clean up timeseries and beach profile data
% * _ct_sitevulnerability_ – function to compute Beach Vulnerability Index for a single site 
% * _ct_sitevulnerability_ – function to compute Beach Vulnerability Index for a set of beach profiles
% * _simBMVfitting_ – functions to fit the BMV model to a profile data set
% * _simYGORmodel_ – functions to run the YGOR model 

%% Additional Coastal Classes
% In the folder ../muiAppCoastalClasses.
%%
% * *ctBeachProfileData* – handle beach profiles as Chainage + Elevation and/or Easting, Northings + Elevation.
% * *ctBlueKenueData* – Blue Kenue file formats for data a point
% * *ctHindcastParameters* – input parameters for the wind-wave model
% * *ctShorelineData* – shoreline position data (e.g. distance to mean tide level)
% * *ctStructureInput* – definition of a beach or structure used to estimate overtopping quantities
% * *ctTidalAnalysis* – analyse water level data to extract tidal constituents and use constituents to construct tidal elevation timeseries.
% * *ctWaterLevelData* – import water level data
% * *ctWaveData* – import wave data
% * *ctWaveModel* – model of wave propagation to a nearshore or deep-water locations
% * *ctWaveParameters* – input parameters for the wave models
% * *ctWindData* – import wind data
% * *ctWindWaveModel* – model wind-wave generation over varying fetch lengths

%% Import format functions
% The folder ../muiAppCoastalClasses/FormatFiles contains a range of 
% format files for loading a number of different data formats (eg waves, 
% water levels, etc) and some generic QC functions.

%% Additional Coastal Functions
% In the folder ../muiAppCoastalFcns.
%%
% * _addslrtotides_ – add a linear, or exponentially varying, rate of sea level rise to a tidal record.
% * _beachtransportratio_ - compute the ratio of onshore to alongshore wave
% action given by tan(alp).
% * _celerity_ - calculate the wave celerity using Hunt's equation.
% * _coeff_AB_ - function called by overtopping functions otop_Q and
% otop_C.
% * _deanbeachprofile_ - find the bed slope across the surf zone.
%   the profile is based on a user defined slope between Hw and SWL (0mOD)
%   and a Dean profile below this level. 
% * _eff_fetch_ - compute the effective fetch for each mean direction based on the 
%   directional distribution function.
% * _get_profile_times_ - get the composite time intervals for all
% profiles.
% * _getalp_ - find the angle between the wave crest and the bed contour.
% * _hallermeier_zones_ - computes the limits of the surf and shoal zones and returning the
% surf zone limit and the closure depth or shoal limit.
% * _hb_break_ - wave height after breaking for given water depth.
% * _hs_break_ - significant wave height after breaking for given water
% depth.
% * _hs_surf_ - calculate the inshore wave height at the edge of the surf
% zone.
% * _iribarren_ - function to calculate the Iribarren number which characterise the 
%   wave breaker type.
% * _littoraldrift_ - sediment transport drift rates for sand and shingle.
% * _littoraldriftstats_ - estimate net drift over selected period from the
% time series of drift rates.
% * _otop_Q_ - function to calculate the overtopping discharge for a simple sloping
%   structure based on the Method proposed by Owen.
% * _posneg_dv_stats_ - computes the rate of change of variable and plot 
% histograms for positive and negative components.
% * _profileslope_ - calculate the bed slope at some depth within the surf
% zone.
% * _readfetchfile_ - read the file that contains the fetch lengths as a
% function of direction.
% * _refraction_ - plane bed wave refraction and shoaling using linear wave
% theory.
% * _rmswaveorbitalvelocity_ - calculate the root mean square wave orbital
% velocity.
% * _runup_ - wave runup magnitude.
% * _runup_slope_ - calculate runup beach slope using Reis A H and Gama C,
% 2010.
% * _section_centroid_ - find centroid and area of a cross-section defined
% by points (xi,zi).
% * _settling_velocity_ - calcualte the settling velocity using Soulsby
% equation.
% * _shoaling_ - plane bed wave shoaling using linear wave theory.
% * _shore_orientation_ - find the orientation for a series of coordinates
% that define a line.
% * _shore_profile_order_ - sort the profile order based on the E,N of the base point - min(Chainage)
% * _simple_tide_ - function to compute a tidal water level time series using the main
% constituents scaled to the required tidal amplitude.
% * _slope_points_ - find point and slope on line (eg a shore profile)
% nearest to zlevel.
% * _sortENdata2line_ - sort the eastings and northings into an order that makes the best
% continuous line
% * _tau_crit_ - calculate the critical erosion shear stress and erosion rate for sand, 
% mud or mixed sediments.
% * _tma_spectrum_ - calculate the TMA spectrum for waves that are depth
% limited.
% * _ucrit_ - compute the critical flow velocity for a given critical shear stress and
% wave conditions in the wave-currrent case.
% * _ut_constants.mat_ - binary file of constituents used by
% *ctTidalAnalysis* class.
% * _ut_reconstr_ - reconstruct superposed harmonics (code from Matlab(TM)
% Forum).
% * _ut_solv_ - execute harmonic tidal analysis (code from Matlab(TM)
% Forum).
% * _waterlevelfreqplots_ - create various water level exceedance and
% duration plots.
% * _wave_energyflux_ - function to calculate the wave energy flux using
% linear wave theory.
% * _wave_friction_ - compute the wave friction factor for rough and smooth turbulent
%   conditions and smooth laminar conditions.
% * _xshore_bailard_ - computes the cross-shore transport for given wave
% and beach conditions.

%% See Also
% See <matlab:doc('coastaltools') CoastalTools>, <matlab:doc('muitoolbox') muitoolbox>
% <matlab:doc('dstoolbox') dstoolbox>
