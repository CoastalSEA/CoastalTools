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
% * _ct_sitevulnerability_ – compute Beach Vulnerability Index for a single site 
% * _ct_beachvulnerability_ – compute Beach Vulnerability Index for a set of beach profiles
% * _ct_beachvariance_ - compute the beach variability over a set of profiles
% * _simBMVfitting_ – fit the BMV model to a profile data set
% * _simYGORmodel_ – run the YGOR model 

%% Additional Coastal Classes
% In the folder ../muiAppCoastalClasses.
%%
% * *ctBeachProfileData* – handle beach profiles as Chainage + Elevation and/or Easting, Northings + Elevation.
% * *ctHindcastParameters* – input parameters for the wind-wave model
% * *ctModelData* – load numerical model output (e.g. Blue Kenue file
% formats)
% * *ctShorelineData* – shoreline position data (e.g. distance to mean tide level)
% * *ctStructureInput* – definition of a beach or structure used to estimate overtopping quantities
% * *ctTidalAnalysis* – analyse water level data to extract tidal constituents and use constituents to construct tidal elevation timeseries.
% * *ctWaterLevelData* – import water level data
% * *ctWaveData* – import wave data
% * *ctWaveParameters* – input parameters for the wave models
% * *ctWindData* – import wind data

%% Import format functions
% The folder ../muiAppCoastalClasses/FormatFiles contains a range of 
% format files for loading a number of different data formats (eg waves, 
% water levels, etc) and some generic QC functions.

%% Coastal Classes and functions for wave models
% In the folder ../muiAppCoastalWaves.
%% 
% *Classes*
%%
% * *ctWaveModel* – model of wave propagation to a nearshore or deep-water
% locations.
% * *ctWaveSpectraPlots* - analyse wave spectra data held as spectral density as a
%   function of direction and frequency, or loaded from a file.
% * *ctWaveSpectrum* - creates and holds a wave spectrum in terms of spectral density 
%   as a function of direction and frequency.
% * *ctWindWaveModel*– model wind-wave generation over varying fetch
% lengths.
% * *waveModels* - Abstract class wave models to define data access methods.
%%
% *Functions*
%%
% * _addwaterlevels2waves_ - class function to add water levels to a 
% selected wave dataset. Used by *waveModels* and *ctWaveData* when being 
% used for derivative models such as runup in *CT_WaveModels*.
% * _datawell_directional_spectra_ - estimates the directional distribution of a wave spectrum for directions, dirs, given the mean, spread, skewness and kurtosis parameters as output by datawell buoys SPT file format.
% * _directional_spreading_ – sample a directional spreading function at selected direction intervals.
% * _extract_wave_data_ - extract Hs, Tp and Dir from a dataset that does not use default naming
%   convention (e.g. Copernicus re-analysis data).
% * _extract_wind_data_ - extract AvSpeed,MaxSpeed,Dir from a dataset that does not use default naming
%   convention.
% * _get_inshore_spectrum_ - construct the offshore and inshore spectra for given wave conditions or wave buoy spectral data.
% * _isangletol_ – boolean check of whether an angle lies between upper and lower bounds defined as specific angles, or a tolerance.
% * _setspectrum_ -reduce a detailed model spectrum to the format defined by the Datawell 
%   buoy spt file format.
% * _subsample_spectra_ts_ - create a timeseries by interpolating a spectrum timeseries using times
%   from another timeseries.
% * _taylor_plot_ts_ - add a timeseries of test points to a Taylor diagram (assumes common
%   reference point).
% * _trapz_weights_periodic_ - integration weights for periodic trapezoidal
% rule.
% * _wave_spectrum_ - calculate the spectral energy at a number of frequencies using a selection of spectrum definitions (Bretschneider open ocean, Pierson-Moskowitz fully developed, JONSWAP fetch limited, and TMA shallow water).
% * _wave_spectrum_gamma_ - estiamte the JONSWAP gamma from a frequency spectrum, S(f), or a
%   frequency–direction spectrum S(dir,f).
% * _wave_spectrum_params_ - integrate a 2-D spectra to obtain wave
% parameters.
% * _wrm_single_animation_ - animation of model spectra timeseries.

%% Additional Coastal Functions
% In the folder ../muiAppCoastalFcns.
%%
% * _addslrtotides_ – add a linear, or exponentially varying, rate of sea level rise to a tidal record.
% * _alterangle_ – adjust the direction angle of a wind or wave dataset. Applies a linear scaling of a direction shift between min and max directions (Dir0 and Dir1) and then the maximum shift outside of this range (Dir2).
% * _beachtransportratio_ - compute the ratio of onshore to alongshore wave
% action given by tan(alp).
% * _binned_wave_climate_ – compute equal energy flux bins for the wave height-direction scatter.
% * _celerity_ - calculate the wave celerity using Hunt's equation.
% * _coeff_AB_ - function called by overtopping functions otop_Q and
% otop_C.
% * _crenulate_bay_ - generate the shoreline for an equilibrium crenulate bay using the method of Hsu and Evans, 1989 (Hsu and Evans, 1989).
% * _ct_coastal_plots_ - functions to do provide additional bespoke plot options using the wave and beach process data in coastal tools.
% * _ct_data_cleanup_ – set of functions to help clean up timeseries and beach profile data
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
% * _lambertw_ - computes the Lambert W-Function. Author: Pascal Getreuer,
% Matlab(TM) Exchange Forum.
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
% * _scale_waterlevels_ - function to scale the water levels based on factors for above and below zero.
% * _section_centroid_ - find centroid and area of a cross-section defined
% by points (xi,zi).
% * _sediment_properties_ – calculate a set of sediment properties based on bulk density.
% * _sedprops_ - function to return one of a range of sediment properties based on selection defined in 'prop'.
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
% continuous line.
% * _tau_bed_ - bed shear stress under combined wave-current action.
% * _tau_crit_ - calculate the critical erosion shear stress and erosion rate for sand, 
% mud or mixed sediments.
% * _tma_spectrum_ - calculate the TMA spectrum for waves that are depth
% limited.
% * _tidehighlow_ - function to compute tidal high and low waters from water level timeseries.
% * _tidalrange_ - compute tidal range from a water level time series.
% * _tidalrange_nltc_ - fit trend and cycles to tidal range or HW/LW time series.
% * _ucrit_ - compute the critical flow velocity for a given critical shear stress and
% wave conditions in the wave-currrent case.
% * _ut_constants.mat_ - binary file of constituents used by
% *ctTidalAnalysis* class.
% * _ut_reconstr_ - reconstruct superposed harmonics. Author: D. Codiga,
% Matlab(TM) Exchange Forum.
% Forum).
% * _ut_solv_ - execute harmonic tidal analysis. Author: D. Codiga,
% Matlab(TM) Exchange Forum.
% * _waterlevelfreqplots_ - create various water level exceedance and
% duration plots.
% * _wave_energyflux_ - function to calculate the wave energy flux using
% linear wave theory.
% * _wave_friction_ - compute the wave friction factor for rough and smooth turbulent
% conditions and smooth laminar conditions.
% * _wave_scatter_ - plot H-T scatter diagram that uses the wave celerity based on water depth 
% to add contours of wave steepness, plot a 3D contoured histogram or
% various combinations of height, period, depth and steepness.
% * _wave_scatter_3d_ - plot H-T-D scatter diagram of wave height, period
% and direction.
% * _wave_steepness_ - compute the wave steepness for given wave height, 
% wave period and water depth.
% * _xshore_bailard_ - computes the cross-shore transport for given wave
% and beach conditions.

%% See Also
% See <matlab:doc('coastaltools') CoastalTools>, <matlab:doc('muitoolbox') muitoolbox>
% <matlab:doc('dstoolbox') dstoolbox>
