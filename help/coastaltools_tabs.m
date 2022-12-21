%% Tab options
% The Tabs display various information such as a list of data sets available (imported data or model output) , a summary of what is 
% currently defined and rapid access to simple plots. Note: these only updated when clicked on using a mouse 
% and values cannot be edited from the Tabs.

%% Summary of Tabs
% * *Data*: lists the data sets that have been loaded with summary meta-data.
% * *Models*: lists data sets that have been derived using hte analysis and 
% modelling tools within CoastalTools.
% * *Site*: tabulates the properties required for model calculations, and 
% is split into two tabs:
% * _Site>Waves_: tabulates the site properties used in the various wave models 
% and structure properties, as used in overtopping calculations (only one 
% definition at any one time).
% * _Site>Simulation_: tabulates the parameter settings required for the YGOR 
% and BMV simulations models.
% * *Q-Plot*: provides rapid access to default plotting of any variable in any data set shown on the Data of Models tabs. These plots cannot be edited and are designed to check data records when loading.
% * *Calcs*: tabulates the results from model analysis, and is split into two tabs:
% * _Calcs>Volumes_: tabulates the results of the most recent beach volume calculations. Includes a button to copy results to the clipboard.
% * _Calcs>Shoreline_: tabulates the results of the most recent shoreline calculations. Includes a button to copy results to the clipboard.
% * _Calcs>Profile_: plot an idealised profile and closure depths, or the Dean profile
% * *Stats*: tabulates the results from statistics UI, and is split into two tabs:
% * _Stats>Descriptive_: tabulates the results of the most recent descriptive statistics. Includes a button to copy results to the clipboard.
% * _Stats>Extremes_: tabulates the results of the most recent extremes analysis. Includes a button to copy results to the clipboard.

%% Accessing Case meta-data
% On the *Data* and *Model* tabs, using the mouse to click on  the fist column of a case record generates a
% table figure with a summary of the meta-data for the selected case. The
% figure includes buttons to Copy the data to the clipboard, view the
% DSproperties of the selected dataset and examine the Source of the
% dataset, which may be alist of files for imported data or details of the
% model used.

%% See Also
% The <matlab:ct_open_manual manual> provides further details of setup and 
% configuration of the model.