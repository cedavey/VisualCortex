%  [layerconfig,doplot, fprefix] = parseVisualCortexInputs(config_file, userinputs)
% 
% Inputs:
%  doplot      - true to plot weights over time & final weights
%  fprefix     - prefix to append to output files
%  configFunctionHandle - function name if config provided in config fn
%                         rather than .mat file
% Outputs:
%  layerconfig - structure containing user network configuration
%  doplot      - plot request (y/n)
%  fprefix     - prefix to append to output files
% 
% Process config info & return layerconfig struct
function [ layerconfig, doplot, fprefix, genTuningCurves ] = parseVisualCortexInputs( config_file, userinputs )
   optargs = { true, 'visual_cortex', [], false }; % doplot, fprefix, configFileHandle
   inargs  = length(userinputs);
   optargs(1:inargs) = userinputs(:);
   [ doplot, fprefix, configFunctionHandle, genTuningCurves ] = optargs{:};
   %% Parse input parameters
   if isempty( genTuningCurves ), genTuningCurves = false; end
   outfile = []; outstruct = []; runSimInputs = true; % inputs allow runSim
   % Try loading most current param config from function 
   if ~isempty(configFunctionHandle)
      if isa(configFunctionHandle, 'function_handle')
         [layerconfig, config_file] = configFunctionHandle(config_file);
      elseif ~exist(config_file,'file') && ~exist([config_file '.m'],'file')
         cprintf('*comment','\nConfig function file - not found, exiting...\n');
         runSimInputs = false;
         return;
      else
         cprintf('*comment','\nConfig function handle type not recognised, exiting...\n');
         runSimInputs = false;
         return;
      end
   else
     % No function handle provided but have filename, so try loading .mat 
     % Config file may contain layerconfig info to run simulation, or it
     % may contain outstruct (& optional layerconfig info) to plot results 
     % only, in which case runSim should be set (or will default) to false.
     if exist('config_file','var') && ~isempty(config_file)
       try 
         % if simulation hasn't been run yet the file will only contain the
         % layerconfig struct. If it's already been run, & we're continuing
         % a simulation, or just plotting it, the file will contain the
         % outstruct structure instead.
         if ischar(config_file)
            load(config_file);
             if ~exist('layerconfig','var') && ~exist('outstruct','var')
               cprintf('*comment','\nInvalid mat file - no config or output structures found, exiting...\n');
               return;
             % layer config params not found, but output struct found -> 
             % look for config params in outstruct, but if not found can
             % plot results, but not run simulation
             elseif ~exist('layerconfig','var')
                if isfield(outstruct,'layerconfig')  && isfield(outstruct,'ntwkconfig')
                   outstruct2variables;
                end
                % if outstruct exists, but not layerconfig -> can only plot
                if layerconfig.runSim
                     cprintf('*comment','\nCannot run simulation because mat file contains output info only, try plotting...\n');
                     runSimInputs = false;
                end
             % don't have output struct, but have layer config -> run sim   
             elseif ~exist('outstruct','var')
                if ~layerconfig.runSim
                   cprintf('*comment','\nCannot plot sim results because mat file contains config info only, will try running simulation...\n');
                   runSimInputs = true;
                end
             end
         elseif isstruct(config_file)
            % config_file is a config structure
            if isfield(config_file,'T')
               layerconfig = config_file;
               config_file = 'visual_cortex';
               
            % config_file is an output structure - either plot results if
            % continueSim is false, or extend the simulation if it's true
            elseif isfield(config_file,'outweights')
               outstruct    = config_file;
               config_file  = outstruct.config_file;
               outstruct2variables;
               % extend the simulation
               if ~isfield(layerconfig,'continueSim')
                  layerconfig.continueSim = false; 
               end
               if layerconfig.continueSim
                  % new simulation time has to be greater than old time
                  if outstruct.outtime(end) >= outstruct.layerconfig.T
                     % PROBLEM: old sim time may not be correctly kept in
                     % outtime(end) if final output sample was not at the
                     % last time step of the simulation
                     str = sprintf('\nTo extend simulation new time (%d) must be greater than old (%d)...\n', ...
                                    outstruct.layerconfig.T, outstruct.outtime(end));
                     cprintf('*comment', str);
                     runSimInputs = false;
%                      return;
                  else
                     runSimInputs = true;
                  end
               % plot the simulation results
               else
                  runSimInputs = false;
               end
            end
         end
       catch ME
         cprintf('Errors','Error loading configuration parameters (%s), exiting...\n',ME.message);
         return;
       end
     % No function handle or .mat file provided, try using default values
     else 
        configFunctionHandle = @visualCortexNoInhibConfig;
        config_file = 'visual_cortex_paramConfig'; % filename to save config to
        try
          layerconfig = configFunctionHandle('visual_cortex_paramConfig');
        catch ME
          cprintf('*comment','Error loading default configuration parameters (%s), exiting...\n',ME.message);
          return;
        end
     end
   end
   % Parse configuration parameters if not done already
   if ~isfield( layerconfig, 'parse' )
     layerconfig = parseVisualCortexConfigFile(layerconfig);
   end
   layerconfig.runSimInputs = runSimInputs;
%    s = fprintf('Using visual cortex simulation configuration file %s\n', config_file); 
%    cprintf('strings','%s',s);
   

   end