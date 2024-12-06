%% Menu Options
% Summary of the options available for each drop down menu.

%% File
% * *New*: clears any existing model (prompting to save if not already saved) and a popup dialog box prompts for Project name and Date (default is current date). 
% * *Open*: existing Asmita models are saved as *.mat files. User selects a model from dialog box.
% * *Save*: save a file that has already been saved.
% * *Save as*: save a file with a new or different name.
% * *Exit*: exit the program. The close window button has the same effect.

%% Tools
% * *Refresh*: updates Scenarios tab.
% * *Clear all > Project*: deletes the current project, including all Setup data and all Cases.
% * *Clear all > Figures*: deletes all results plot figures (useful if a large number of plots have been produced).
% * *Clear all > Data*: deletes all cases listed on the _Data_ tab but does not affect Site setup.
% * *Clear all > Models*: deletes all cases listed on the _Models_ tab but does not affect Site setup.

%% Project
% * *Project Info*: edit the Project name and Date
% * *Cases > Edit Description*: user selects a Case to edit Case description.
% * *Cases > Edit DS properties*: initialises the  UI for editing Data Set properties (DSproperties).
% * *Cases > Edit Data Set*: initialises the Edit Data UI for editing data sets.
% * *Cases > Save*: user selects a data set to be saved from a list box of scenarios and the is then prompted to name the file. The data are written to an Excel spreadsheet. 
% * *Cases > Delete*: user selects the Case(s) to be deleted from a list box of Cases and these are then deleted (model setup is not changed).
% * *Cases > Reload*: user selects a Case to reload as the current parameter settings.
% * *Cases > View settings*: user selects a Case to display a table listing the parameters used for the selected Case. 
% * *Export/Import > Export*: user selects a Case class instance to export as a mat file.
% * *Export/Import > Import*: user selects an exported Case class instance (mat file) to be loaded.
%%
% *NB*: to export the data from a Case for use in another application 
% (eg text file, Excel, etc), use the *Project>Cases>Edit Data Set* option 
% to make a selection and then use the ‘Copy to Clipboard’ button to paste 
% the selection to the clipboard.

%% Setup > Import Data
% Load data for Waves, Water Levels, Winds, Beach profiles, Shorelines,
% BlueKenue and User defined data sets. To create a new instance (e.g. for a
% different location or data source) use Load. To add data to an
% existing data set, use Add.
%%
% * *Load data*: prompts user for file format to be loaded. The options available vary with Data type. The user then selects one or more files to be loaded. Once files have been read, user is prompted for a description (working title) for the data set. 
% * *Add data*: prompts user for file to be loaded (only one file at a time can be added). Only files with the format used to create the data set can be used to Add data to a data record.
% * *Quality control*: runs a series of checks on the data. This is only available for waves and water levels in this version. See <matlab:ct_open_manual manual> for details of the checks applied. Once run, the last column of the Data tab is updated to show that qc has been applied.

%% Setup > Site Parameters
% * *Site parameters>Wave propagation*: input dialogue allows user to define characteristic properties of the site. The current values can be viewed on the _Site/Waves_ tab.
% * *Site parameters>Wind-wave hindcast*: input dialogue to define model parameters for the wind-wave hindcast model. The current values can be viewed on the _Site/Waves_ tab.
% * *Structure parameters*: input dialogue to define the parameters used when calculating wave overtopping. The current values can be viewed on the _Site/Waves_ tab.

%% Setup > Model Parameters
% * *YGOR simulation parameters*: used to define the fit coefficients to be used in YGOR simulations. The current values can be viewed on the _Site/Simulation_ tab.
% * *BMV simulation parameters*: used to define the parameters needed to run the BMV simulation. The current values can be viewed on the _Site/Simulation_ tab.
% * *Model constants*: Various constants are defined for use in models, such as the acceleration due to gravity, viscosity and density of sea water, and density of sediment. 

%% Setup > Data Clean-up
% * *Concatenate two timeseries*: allows two timeseries data sets to be joined. Two records, or two variables, of a similar type can be joined to form a single timeseries or timeseries collection. 
% * *Resample timeseries*: allows a selected timeseries to be resampled at user specified interval, using a user specified method (e.g. mean/max/min over the interval). 
% * *Patch timeseries*: allows gaps in a selected timeseries to be patched using the data from another timeseries that overlaps the primary timeseries (at least for some or all of the gaps). 
% * *Trim timeseries*: allows the start and end dates of a timeseries to be modified.
% * *Delete interval*: set an interval within a selected variable to NaN. This maintains the length of the record for plotting, etc but removes the data. Can be useful when cleaning up tidal data and similar time series where this option can be used in conjunction with the Patch option to create a "clean" timeseries.
% * *Merge cases*: some cleanup functions only work on one variable at a time and the data are then saved as a new case. This function compiles several variables back into a single case (e.g., when concatenating Hs, Tp and Dir for use in models).
% * *Delete multiple profiles*: allows all profile data sets with less than ‘N’ surveys available to be deleted. 
% * *Edit or Delete profile in timeseries*: options to scroll through all profiles in the data set of select individual profiles and then to either edit or delete individual profiles. For editing profiles this is easier to use than the _Project>Scenarios>Edit Data Set_ option.

%% Run > Wave properties
% Models that use the time series data to compute
% various wave related properties. All of these models require the Site parameters to have been defined. 
%%
% * *Deepwater Waves*: use an input wave data set (and water levels if
% avaialble to compute the deepwater waves using linear wave theory for shoaling and 
% plane bed refraction.
% * *Nearshore Waves*: use an input wave data set (and water levels if
% avaialble to compute the nearshore waves using linear wave theory for shoaling and 
% plane bed refraction. This can be to a specified still water depth, or
% the edge of the surf zone.
% * *Wind-waves*: uses an input wind data set and defined fetch lengths. The effective fetch ic calculated using the method of Donelan (1985), or the method given in the Shore Protection Manual (1992). Wave height and period are computed using the effective fetch lengths and the TMA spectrum (Hughes, 1984; Bouws et al., 1985; 1987). 
% * *Spectral Transfer*: uses a Spectral Transfer table from the
% WaveRayModel to transfer an offshore timeseries to an inshore point.
% Requires a SpectralTransfer table to have been imported from the WaveRayModel
% App.
% * *Wave Energy*: uses selected wave data set (input or modelled) to
% compute the wave energy flux using linear wave theory.
% * *Runup*: uses selected wave data set (input or modelled) to
% compute wave run-up using the formulation of Stockdon et al (2006). 
% * *Littoral Drift*: uses selected wave data set (input or modelled) to
% compute the littoral drift using one of 4 options: US Army Corps CERC
% formula ((1992); SANDS formula; Kamphuis formula (1991) and the Damgaard
% and Soulsby formula (1997).
% * *X-shore Transport*: uses selected wave data set (input or modelled) to
% compute cross-shore transport using the method proposed by Bailard and Inman (1981). 
% * *Overtopping*: uses selected wave data set (input or modelled) and water
% levels to compute the volume of water overtopping a sea wall, or beach crest, 
% using the overtopping formula proposed by Owen (1980). 
% * *Iribarren Number*:  uses selected wave data set (input or modelled) to
% compute the Iribarren number and breaker type.

%% Run > Beach properties
% Models that use the beach profile data to compute beach volume and shoreline change
%%
% * *Profiles*: analysis options for indivivual profiles
% * _Profiles > Volumes_: compute the volume per metre run of beach, within
% a user defined control box (minimum lower level and shoreward position).
% Also saves the non-dimesional volume and centroid position.
% * _Profiles > Shoreline position_: compute the change in shoreline position at
% a user specified elevation. Also saves the beach slope at the defined
% position.
% * _Profiles > Location plot_: plot the location of the profiles
% available for analysis.
% * _Profiles > Centroid plot_: plot the time varying movement of the
% profile centroid.
% * _Profiles > Space-time plot_: plot the time varying change in a beach
% property such as volume, beach slope or shoreline position as a surface with 
% axes of time and alongshore position. Requires the parameter to have been
% computed for all profiles to be included in the plot.
%%
% * *Shore change*: analysis options for a set of profiles that make a
% shoreline
% * _Shore change > Shoreline_: compute the change in shoreline position at
% a user specified elevation for a set of profiles.
% * _Shore change > Change plot_: plot of shoreline change as E,N line plot, or
% the distance, slope or orentation as a function of alongshore position
% * _Shore change > Rates plot_: plot rates of change in shoreline position
% and beach slope.
%%
% * *Beach type*: compute the Beach Type from the dimensionless fall velocity (Ω) for the grain size defined in Site Parameters and the selected wave time series 
% * *Shore profile*: plots the profile based on the parameters defined for
% the BMV analysis
% * *Dean profile*: an idealised beach profile based on a linear slope (upper beach slope) between the beach crest level and mean tide level (SWL=0mOD) and an equilibrium profile (Dean, 1977) out to deep water. 

%% Run > Tidal analysis
% Models to derive tidal constituents and generate a simulated tidal record
%% 
% * *Analysis*: a list box prompts the user to select a water level data
% set and the tidal consitutents are computed and saved along with the
% tidal record for the period of the input data set.
% * *Reconstruction*: a list box prompts the user to select an existing tidal analysis data set and then to define the data range and time interval to be used for the predictions. The saved constituents are used to predict the tides for the defined date range and 

%% Run (other options)
% Other modelling tools to derive new data sets, simulate shoreline change
% and estimate beach vulnerability
%%
% * *Derive Output*: initialises the Derive Output UI to select and define manipulations of the data or call external functions and load the result as new data set.
% * *Simulation*: opens the UI for two simulations tools. One fits the
% incident wave energy to shorline changes following the approach of Miller & Dean,
% 2004 and Yates, et al, 2009. The other seeks to track the centroid of an
% idealised beach profile against the actual movement of the beach
% centroid.
% * *Vulnerability*: tools to assess beach vulnerability index for a
% shoreline
% * _Vulnerability > BVI site_: uses ouput from the models for runup, littoral 
% drift, cross-shore transport and shoreline position change rates to
% estimate the average BVI for a site. All input models need to have been
% run before using this option
% * _Vulnerability > BVI profile set_: the user selects a number of
% profiles to define a shoreline and the models needed to determine the
% variation at each profile are run (runup, littoral drift, x-shore drift and
% shoreline change at defined elevation). The shoreline change outputs can
% be viewed using the _Change plot_ and _Rates plot_ options.
% * _Vulnerability > BVI set plot_: generates a plot of beach vulenrability
% index and the contributing indices, as a function of beach profile location along the shore.
% * *User model*: option provided for the user to add their own models. The 
% demonstration code retrieves inshore wave data and returns the inshore 
% wave energy and the ratio of alongshore to cross-shore transport.

%% Analysis
% * *Plots*: initialises the Plot UI to select variables and produce various types of plot. The user selects the Case, Dataset and Variable to used, along with the Plot type and any Scaling to be applied from a series of drop down lists, 
% * *Statistics*: initialiss the Statistics UI to select data and run a range of standard statistical methods.
% * *Output Definitions*: dialogue to edit the meta-data used to desctibe imported data or model outputs.

%% Help
% * *Help*: access the online documentation for CoastalTools.

%% See Also
% The <matlab:ct_open_manual manual> provides further details of setup and 
% configuration of the model.