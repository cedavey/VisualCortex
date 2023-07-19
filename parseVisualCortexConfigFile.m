% config = parseVisualCortexConfigFile(config_file,verbose)
% Inputs:
%   config  - config structure obtained from config file
%   verbose - if true tell user when missing default values
% Outputs:
%   config  - config structure containing all required values for simulation
function config = parseVisualCortexConfigFile(config,verbose)
  if nargin<2, verbose = true; end
  if verbose, fprintf('\nParsing visual cortex layer config options...\n'); end
  if ~exist('config','var')
      cprintf('*Errors',fprintf('Must input a config structure, exiting...\n'));
      return
  end
  layers  = config.layernames;
  
  % Condition checks for inputs
  nonNegInteger   = @(x) all(isnumeric(x)) && all(isint(x)) && all((x>=0)); % non-neg int, scalar
  nonNegFloat     = @(x) all(isnumeric(x)) && all((x>=0)); % non-neg float, scalar
  positiveInteger = @(x) all(isnumeric(x)) && all(isint(x)) && all(x>0); % positive int, scalar
  positiveFloat   = @(x) all(isnumeric(x)) && all((x>0));  % positive float
  [fromDistList,distList] = getDistListAndCondition;
  distCondition = @(x) ischar(x) && any(strcmpi(x,distList));

  [fromPSPList, pspList]  = gePSPListAndCondition;
  fromLayerList = @(x) ischar(x) && any(strcmpi(x,layers));
  fromLayerListOrEmpty = @(x) (ischar(x) && any(strcmpi(x,layers))) || isempty(x);
  
  % Optional fields - name of field, default value, printing, other conditions
  config = checkOptField(config,'runSim',       true,  verbose, @islogical);
  config = checkOptField(config,'continueSim',  false, verbose, @islogical);
  config = checkOptField(config,'evolveSep',    false, verbose, @islogical);
  config = checkOptField(config,'thresholdWgts',true,  verbose, @islogical);
  config = checkOptField(config,'estCov',       true,  verbose, @islogical);
  config = checkOptField(config,'initCovWithAnal',true,verbose,@islogical);
  config = checkOptField(config,'calcCorr',     false, verbose, @islogical);
  config = checkOptField(config,'initTime',     10,    verbose, nonNegInteger);
  config = checkOptField(config,'interval',     10,    verbose, nonNegInteger);
  config = checkOptField(config,'intervalQ',    10,    verbose, nonNegInteger);
  config = checkOptField(config,'normQ',        true,  verbose, @islogical);
  config = checkOptField(config,'endpoints',    true,  verbose, @islogical);
  config = checkOptField(config,'b',            1,     verbose, positiveFloat);
  config = checkOptField(config,'beta',         10,    verbose, positiveFloat);
  config = checkOptField(config,'delta',        1,     verbose, positiveFloat);
  config = checkOptField(config,'connparams',   [],    verbose, positiveFloat);
  config = checkOptField(config,'psptype',     'delta',verbose, fromPSPList);
  config = checkOptField(config,'pspparams',    0.15,  verbose, positiveFloat);
  config = checkOptField(config,'pspdelay',     1,     verbose, positiveInteger);
  config = checkOptField(config,'retina',      [],     verbose, fromLayerListOrEmpty);
  config = checkOptField(config,'retina_on',   [],     verbose, fromLayerListOrEmpty);
  config = checkOptField(config,'retina_off',  [],     verbose, fromLayerListOrEmpty);
  config = checkOptField(config,'velocity',     0,     verbose, nonNegFloat);
  config = checkOptField(config,'noisedist','gaussian',verbose, distCondition);
  config = checkOptField(config,'noiseparams',  [2 1], verbose, nonNegFloat);
  checkDistParams(config.noisedist, config.noiseparams,'noise', nonNegFloat);

  checkPSPParams(config.psptype,config.pspparams,config.pspdelay,positiveFloat);
 
  % Required fields
  checkReqField(config,'dt',     'config', verbose, positiveFloat);
  checkReqField(config,'T',      'config', verbose, positiveFloat);
  checkReqField(config,'inlayer','config', verbose, []);
  
  if isfield(config,'radialconf')
     radialconf = config.radialconf;
  else
     radialconf.periodic     = true;
     radialconf.radialConn   = []; 
%      radialconf.unequalVar   = false;           % radius configuration options
%      radialconf.rayleighDist = false;           % Rayleigh distributed distance (else Gaussian)
%      radialconf.radialSquRoot= true ;              
  end
  radialconf = checkOptField(radialconf,'periodic',     true, verbose, @islogical);
  radialconf = checkOptField(radialconf,'radialDist',   false,verbose, @islogical);
  radialconf = checkOptField(radialconf,'rayleighDist', true, verbose, @islogical);
  radialconf = checkOptField(radialconf,'radialSquRoot',false,verbose, @islogical);
  radialconf = checkOptField(radialconf,'unequalVar',   false,verbose, @islogical);
  radialconf = checkOptField(radialconf,'uniqueLoc',    false,verbose, @islogical);
  config.radialconf = radialconf;

  %% Check layers
  for li=1:length(layers)
    config = checkLayerInfo(config,layers{li},verbose,false);
  end
  if isfield(config,'inlayer')
    config = checkLayerInfo(config,config.inlayer,verbose,true);
  end
  
  %% Check connections
  % Get connections and check format of variable
  conns  = config.layerconns; % e.g. {'A','B'; 'B','C';};
  if ~isDim(conns,2) && size(conns,2)==2
    errStr = sprintf('Invalid format of layer connections variable (%d dimensions with second dimension size of %d',...
                      isDims(conns), size(conns,2));
    error('\n\t%s\n',errStr);
  end
  % Check that each layer appears somewhere in the connection list - else
  % the layer is not connected to anything & is useless!
  inlist = cellfun(@(lstr) any(strcmpi(lstr,conns(:))), layers);
  if ~all(inlist)
    errStr = sprintf('Layer %s does not appear in the connections',...
                     layers(find(inlist,1,'first')));
    error('\n\t%s\n',errStr);
  end
  % Check format of each connection
  for ci=1:size(conns,1)
    cxlabel = strjoin(conns(ci,:),'');
    config  = checkConnectionInfo(config,cxlabel,verbose);
  end
  
  config    = checkRFInfo(config,verbose);
  
  % If we've made it here then there haven't been any errors
  config.parse = true;
  if verbose, fprintf('...finished parsing layer configuration\n\n'); end

  % Check that there are sufficient neurons in a layer for the number of
  % required connections
end

% Check receptive field parameters of retina layer
function config = checkRFInfo(config,verbose)  % Condition checks for 
  if ~isfield(config,'RF') || isempty(config.RF)
     config.RF = [];
     return;
  end
  RF = config.RF; % parent receptive field structure
  positiveInteger = @(x) isnumeric(x) && isint(x) && all(x>0); % positive int, scalar
  nonNegInteger = @(x) isnumeric(x) && isint(x) && all(x>=0); % positive int, scalar
  positiveFloat = @(x) isnumeric(x) && all(x>0);  % positive float
  nonNegFloat   = @(x) isnumeric(x) && all(x>=0);  % positive float
  arrayFloat    = @(x) isnumeric(x) && isDim(x,2);
  biasDistStr   = @(x) ischar(x) && any(strcmpi(x,{'levick','none','uniform'}));
  normFnStr     = @(x) ischar(x) && any(strcmpi(x,{'power','sum','std','mean'}));
  angleVector   = @(x) iscell(x) && all(cellfun(@(in)  (ischar(in) && strcmpi(in,'radial')) ...
                                                    || (isnumeric(in) && isscalar(in) && (in<=180) && (in>=0)),...
                                        x)); 
  probVector    = @(x) isnumeric(x) && isvector(x) && all(x<=1) && all(x>=0);
  [fromDistList,distList] = getDistListAndCondition;

  % Required fields
  checkReqField(RF,'kernel_size',  'RF', verbose, positiveInteger);
  checkReqField(RF,'sd_surr_hor',  'RF', verbose, positiveFloat);
  checkReqField(RF,'sd_surr_ver',  'RF', verbose, positiveFloat);
  checkReqField(RF,'sd_cent_hor',  'RF', verbose, positiveFloat);  
  checkReqField(RF,'sd_cent_ver',  'RF', verbose, positiveFloat);
  checkReqField(RF,'coeff_cent',   'RF', verbose, positiveFloat);
  checkReqField(RF,'coeff_surr',   'RF', verbose, positiveFloat);
  checkReqField(RF,'ratio_surr_cent', 'RF', verbose, positiveFloat);
  checkReqField(RF,'retina_kernel','RF', verbose, arrayFloat);
  checkReqField(RF,'input_var',    'RF', verbose, positiveFloat);
  checkReqField(RF,'angleSet',     'RF', verbose, angleVector);
  checkReqField(RF,'angleProb',    'RF', verbose, probVector);
  if ~(length(RF.angleSet) == length(RF.angleProb))
    errStr = sprintf('All RF angles must have 1 probability assigned to them');
    error('\n\t%s\n',errStr);
  end
  if sum(RF.angleProb)<0.999 || sum(RF.angleProb)>1.001
    errStr = sprintf('RF angle probabilities must sum to 1');
    error('\n\t%s\n',errStr);
  end
   
  config.RF = checkOptField(RF,'normfn',  'sum',  verbose, normFnStr);
  config.RF = checkOptField(RF,'biasdist','none', verbose, biasDistStr);
end

function config = checkConnectionInfo(config,cxlabel,verbose)
  [fromDistList,distList]     = getDistListAndCondition;
  [fromConnList,connList]     = getConnectionListAndCondition;
  [fromRadiusList,radiusList] = getRadiusListAndCondition;
  % connList  = {'fixed_norm'};
  % Condition checks for inputs
  nonNegInteger   = @(x) isnumeric(x) && isint(x) && (x>=0); % non-negative scalar integer
  positiveInteger = @(x) isnumeric(x) && isint(x) && (x>0);  % positive scalar integer
  stndFloat       = @(x) isnumeric(x) && isscalar(x) && (x>=0) && (x<=1);
  nonNegFloat     = @(x) isnumeric(x) && all(x>=0); % non-negative float
  float           = @(x) isnumeric(x);
  positiveFloat   = @(x) isnumeric(x) && all(x>0);  % positive float
  scalarFloat     = @(x) isnumeric(x) && isscalar(x);
  fromDistList    = @(x) ischar(x) && any(strcmpi(x,distList));
  fromConnList    = @(x) ischar(x) && any(strcmpi(x,connList));
  fromRadiusList    = @(x) ischar(x) && any(strcmpi(x,radiusList));
    
  % allow radius to be randomised rather than constant
  if ~isfield(config, 'plastic')
     config.plastic.(cxlabel) = [];
  end
  config.plastic = checkOptField(config.plastic, cxlabel, true, verbose, @islogical);
  % Check delay variables exist, & check specification of distribution params
  if ~isfield(config, 'deldist')
     config.deldist.(cxlabel) = [];
  end
  config.deldist = checkOptField(config.deldist, cxlabel, 'constant', verbose, fromDistList);
  if ~isfield(config, 'delparams')
     config.delparams.(cxlabel) = [];
  end
  config.delparams = checkOptField(config.delparams, cxlabel, 0, verbose, nonNegFloat);
  
  % allow user to input a single layerdist for all layers
  if ~isfield(config, 'layerdist')
     config.layerdist.(cxlabel) = 0;
     defLayerDist     = 0;
  elseif isnumeric( config.layerdist )
     defLayerDist     = config.layerdist;
     config.layerdist = [];
     config.layerdist.default   = defLayerDist;
     config.layerdist.(cxlabel) = defLayerDist;
  elseif ~isfield(config.layerdist, cxlabel)
     if isfield(config.layerdist, 'default')
         defLayerDist = config.layerdist.default;
     else
         defLayerDist = 0;
     end
  end    
  config.layerdist = checkOptField( config.layerdist, cxlabel, defLayerDist, verbose, nonNegFloat );
  config           = checkOptField( config, cxlabel, 0, verbose, nonNegFloat );
  checkDistParams( config.deldist.(cxlabel), config.delparams.(cxlabel), 'delays', nonNegFloat );
  
  % Check synaptic connection deets
  checkReqField(config.N,cxlabel,'N',verbose,positiveInteger);
  checkReqField(config.r,cxlabel,'r',verbose,nonNegFloat);
  if ~isfield(config, 'radiusdist')
     config.radiusdist.(cxlabel) = [];
  end
  config.radiusdist   = checkOptField(config.radiusdist,  cxlabel,'constant',verbose, fromRadiusList);
  if ~isfield(config, 'radiusparams')
     config.radiusparams.(cxlabel) = [];
  end
  config.radiusparams = checkOptField(config.radiusparams,cxlabel, 0, verbose, nonNegFloat);
  if ~isfield(config, 'conndist')
     config.conndist.(cxlabel) = [];
  end
  config.conndist     = checkOptField(config.conndist,    cxlabel,'fixed_norm',verbose,fromConnList);
  if ~isfield(config, 'connparams')
     config.connparams.(cxlabel) = [];
  end
  config.connparams   = checkOptField(config.connparams,  cxlabel, 0, verbose, nonNegFloat);

  % Check weight initialisation distribution
  if ~isfield(config, 'wgtdist')
     config.wgtdist.(cxlabel) = [];
  end
  config.wgtdist   = checkOptField( config.wgtdist, cxlabel, 'constant', verbose, fromDistList );
  if ~isfield(config, 'wgtparams')
     config.wgtparams.(cxlabel) = [];
  end
  config.wgtparams = checkOptField( config.wgtparams, cxlabel, 0.5, verbose, float );
  checkDistParams( config.wgtdist.(cxlabel), config.wgtparams.(cxlabel), 'weights', float );

   % Check learning parameters
  if config.plastic.(cxlabel)
     checkReqField(config.ka, cxlabel, 'ka', verbose,scalarFloat);
     checkReqField(config.kb, cxlabel, 'kb', verbose,positiveFloat);
     checkReqField(config.k1, cxlabel, 'k1', verbose,scalarFloat);
     checkReqField(config.k2, cxlabel, 'k2', verbose,scalarFloat);
     checkReqField(config.n, cxlabel,  'n',  verbose,stndFloat);
     checkReqField(config.Rb, cxlabel, 'Rb', verbose,positiveFloat);
  end
  
  % Check number of synapses does not exceed number of presynaptic neurons
  if length(cxlabel)>2
    cprintf('*Errors','\n\tparseVisualCortexConfigFile assumes layer names are single characters...cannot validate number of synapses\n');
  end
  pre  = cxlabel(1);
  post = cxlabel(2);
  if prod(config.sz.(pre)) <= config.N.(cxlabel)
    warnStr = sprintf('Connection %s has more (or equal) synapses than there are neurons in layer %s (%d neurons)',...
                      cxlabel, pre, prod(config.sz.(pre)));
    disp(sprintf('\n\t%s\n',warnStr));
  end
end

% Check that the number of parameters provided for the specified
% distribution is appropriate. Optionally provide additional conditions in
% a function handle. Note that all parameters must satisfy the condition
% function.
function checkDistParams(dist,params,paramStr,conditions)
  % check that the parameter format is appropriate for the distribution type 
  valid = true;
  errStr = [];
  switch lower(dist)
    case 'constant'
      if ~isscalar(params)
        expect = 'scalar'; 
        valid = false;
      end
    case 'gaussian'
      if ~length(params)==2
        expect = '2 elmt vector'; 
        valid = false;
      end
    case 'natural'
      if ~length(params)==1
        expect = '1 elmt vector'; 
        valid = false;
      end
    case 'uniform'
      if ~length(params)==2
        expect = '2 elmt vector'; 
        valid = false;
      end
  end
  if ~valid
    errStr = sprintf('Invalid format for %s %s distribution parameters (expected %s, received %d)',...
                      upper(paramStr),dist,expect,isDim(params));
  end
  if exist('conditions','var') && ~isempty(conditions) && ~all(arrayfun(conditions,params))
    valid = false;
    errStr = sprintf('%s\n\tInvalid param type for %s %s distribution parameters',...
                      errStr,upper(paramStr),dist);
  end
  if ~valid
    error('\n\t%s\n',errStr);
  end
end

% Check that the number of parameters provided for the EPSP function
% generator distribution is appropriate. Optionally provide additional 
% conditions in a function handle. Note that all parameters must satisfy 
% the condition function.
function checkPSPParams(dist,params,delay,conditions)
  % check that the parameter format is appropriate for the distribution type 
  valid = true;
  errStr = [];
  switch lower(dist)
    case 'delta'
      if ~isscalar(params)
        expect = 'scalar'; 
        valid = false;
      end
      if delay>1
        errStr = cprintf('*Errors','\n\tEPSP is a delta fn, but pspdelay is longer than 1 time sample (%d samples)\n',...
                         delay,true);
        error('\n%s\n',errStr);
      end
    case {'exp','single'}
      if ~isscalar(params)
        expect = 'scalar'; 
        valid = false;
      end
      if delay==1
        errStr = crintf('*Errors','\n\tEPSP is an exponential, but pspdelay is only 1 time sample\n',true);
        error('\n%s\n',errStr);
      end
    case 'double'
      if numel(params)~=2
        expect = '2 elmt vector'; 
        valid = false;
      end
      if delay==1
        errStr = cprintf('*Errors','\n\tEPSP is a double exponential, but pspdelay is only 1 time sample\n',true);
        error('\n%s\n',errStr);
      end
  end
  if ~valid
    errStr = sprintf('Invalid format for %s %s distribution parameters (expected %s, received %d dimensional input with %d elmts)',...
                      upper('psp'),dist,expect,isDim(params),numel(params));
  end
  if exist('conditions','var') && ~isempty(conditions) && ~all(arrayfun(conditions,params))
    valid = false;
    errStr = sprintf('%s\nInvalid param type for %s %s distribution parameters',...
                      errStr,upper('psp'),dist);
  end
  if ~valid
    error('\n\t%s\n',errStr);
  end
end


% Check that all layer info is available for each layer. If it's input
% rather than a fully fledged layer slightly less info is required
function config = checkLayerInfo(config, lname, verbose, isinput)
  % Different condition checks on inputs
  posIntegerVector = @(x) isnumeric(x) && isint(x) && length(x)==2; % non-negative scalar
  nonNegFloat   = @(x) isnumeric(x)  && (x>=0); % non-negative float

  % Required fields
  checkReqField(config.sz,    lname,'sz',    verbose,posIntegerVector);
  checkReqField(config.lambda,lname,'lambda',verbose,nonNegFloat);
  
  % Optional fields
  if ~isinput
     if ~isfield(config, 'spatialJitter')
        config.spatialJitter.(cxlabel) = [];
     end
    config.spatialJitter = checkOptField(config.spatialJitter, lname, 0, verbose, nonNegFloat);
  end
  if ~isfield(config.radialconf, 'radialConn')
     config.radialconf.radialConn = [];
  end
  config.radialconf.radialConn = checkOptField( config.radialconf.radialConn, lname, true, verbose, @islogical);
end

function config = checkOptField(config, name, defVal, verbose, conditions)

	if ~isfield(config, name)
      config.(name) = defVal;
      if verbose
         if     isnumeric(defVal), defStr = num2str(defVal); 
         elseif islogical(defVal), defStr = num2str(defVal); 
         elseif ischar(defVal),    defStr = defVal; 
         end
         errStr = sprintf('%s not found, set to %s',name, defStr);
         cprintf('Comments',sprintf('\t%s\n',errStr));
      end
   elseif isempty( config.(name) )
      config.(name) = defVal;
      if verbose
         if     isnumeric(defVal), defStr = num2str(defVal); 
         elseif islogical(defVal), defStr = num2str(defVal); 
         elseif ischar(defVal),    defStr = defVal; 
         end
         errStr = sprintf('%s empty, set to %s',name,defStr);
         cprintf('Comments',sprintf('\t%s\n',errStr));
      end
   elseif exist('conditions','var') && ~isempty(conditions)
      if ~conditions( config.(name) )
         if isnumeric( config.(name) ) || islogical( config.(name) )
            defStr = num2str( config.(name) ); 
         elseif ischar( config.(name) )
           defStr = config.(name); 
         end
         errStr = cprintf('*Strings','%s value of %s is invalid by condition %s)\n',...
                          name,defStr,func2str(conditions),true);
         error('\n\t%s\n',errStr);
      end
	end
end

% config     - configuration structure with field values
% name       - name of field in structure
% label      - name of config structure (for embedded fields it's useful for
% verbose    - if true print warnings
%              identifying which struct level the problem is at)
% conditions - fn handle checking that conditions are met by field value
function checkReqField(config,name,label,verbose,conditions)
  if ~isfield(config,name)
    errStr = sprintf('%s.%s not found',label,name);
    error('\n\t%s\n',errStr);
  elseif exist('conditions','var') && ~isempty(conditions)
    if ~conditions(config.(name))
      errStr = sprintf('%s.%s value is invalid by condition %s)\n',label,name,func2str(conditions));
      error('\n\t%s\n',errStr);
    end
  end
end

function [connCondition,connList] = getConnectionListAndCondition()
  connList      = {'uniform', 'constant', 'gaussian'};
  connList      = {'fixed_norm', 'unfixed_norm'};
  connCondition = @(x) ischar(x) && any(strcmpi(x,connList));
end

function [fromPSPList,psptypeList] = gePSPListAndCondition()
  psptypeList   = {'delta','exp','single','double'};
  fromPSPList   = @(x) ischar(x) && any(strcmpi(x,psptypeList));
end

function [distCondition,distList] = getDistListAndCondition()
  % natural = natural image input rather than gauss or poisson, for example
  distList      = {'constant','deterministic','gaussian','natural','poisson','uniform'};
  distCondition = @(x) ischar(x) && any(strcmpi(x,distList));
end

function [distCondition,distList] = getRadiusListAndCondition()
  % natural = natural image input rather than gauss or poisson, for example
  distList      = {'constant','gaussian','uniform'};
  distCondition = @(x) ischar(x) && any(strcmpi(x,distList));
end















