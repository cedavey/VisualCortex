curr_dir = pwd;

doFF_v2 = true ;
doL_v2  = true ;
doFF_v1 = false; 
doL_v1  = false;

% calculate number of seconds rather than percentage, since should take
% same amount of time for layer B to evolve regardless of FF or L

amplitude    = 100;
dc_offset    = 200;
feedforward  = {'A', 'B', 'C'};
genMovies    = false;
iterationsL  = { 'start', 'perc23', 'end' };
iterationsFF = { 'start', 'perc6', 'end' };
reducedItL   = { 'perc10','perc15' };
reducedItFF  = { 'perc9', 'perc12' };
lateral      = {'A', 'B', 'C', 'L'};
max_freq     = 0.3;
min_phases   = 7;
num_angles   = 18;
num_freqs    = 10;
plottypes    = {'orient','freq'};
responsefn   = { 'f1', 'max' };

optargsL     = { 'genMovies',  genMovies,   ...
                 'responsefn', responsefn, 'iterations', iterationsL, ...
                 'num_freqs',  num_freqs,  'max_freq',   max_freq,    ...
                 'min_phases', min_phases, 'num_angles', num_angles,  ...
                 'dc_offset',  dc_offset,  'amplitude',  amplitude };
optargsFF    = { 'genMovies',  genMovies,   ...
                 'responsefn', responsefn, 'iterations', iterationsFF,...
                 'num_freqs',  num_freqs,  'max_freq',   max_freq,    ...
                 'min_phases', min_phases, 'num_angles', num_angles,  ...
                 'dc_offset',  dc_offset,  'amplitude',  amplitude };
redargsL     = { 'genMovies',  genMovies,   ...
                 'responsefn', responsefn, 'iterations', reducedItL, ...
                 'num_freqs',  num_freqs,  'max_freq',   max_freq,    ...
                 'min_phases', min_phases, 'num_angles', num_angles,  ...
                 'dc_offset',  dc_offset,  'amplitude',  amplitude };
redargsFF    = { 'genMovies',  genMovies,   ...
                 'responsefn', responsefn, 'iterations', reducedItFF,...
                 'num_freqs',  num_freqs,  'max_freq',   max_freq,    ...
                 'min_phases', min_phases, 'num_angles', num_angles,  ...
                 'dc_offset',  dc_offset,  'amplitude',  amplitude };

if doFF_v2
   try
      fprintf('\n\nStart generating feedforward tuning curves v2\n');
      cd /Users/cedavey/Dropbox/code/visual_cortex/results/feedforward/catDistribution/random_radius_N
   %    [ outstructFF, outfileFF, tuningFF, curvesFF ] = visual_cortex( [], true, 'feedforward', @visualCortexFeedForwardConfig, false );
%       tuningFF  = processVisualCortexInput( outstructFF, false, [], optargsFF );
      redtunFF  = processVisualCortexInput( outstructFF, false, [], redargsFF );
      curves    = generateVisualCortexTuningCurves( redtunFF, feedforward, responsefn, reducedItFF, true, plottypes, true );
      fprintf(' Finished simulating feed forward network and generating tuning curves v2\n');
   catch ME
      fprintf(' Error simulating feed forward network and generating tuning curves v2\n');
      getCatchMEstring(ME);
   end
end

if doL_v2
   try
      fprintf('\n\nStart generating lateral tuning curves v2\n');
      cd /Users/cedavey/Dropbox/code/visual_cortex/results/lateral/catDistribution/random_radius_N
   %    [ outstructL, outfileL, tuningL, curvesL ] = visual_cortex( [], true, 'lateral', @visualCortexWithInhibConfig, false );
%       tuningL   = processVisualCortexInput( outstructL, false, [], optargsL );
%       curves    = generateVisualCortexTuningCurves( tuningL, lateral, responsefn, iterationsL, true, plottypes, true );
      redtunL   = processVisualCortexInput( outstructL, false, [], redargsL );
      curves    = generateVisualCortexTuningCurves( redtunL, lateral, responsefn, reducedItL, true, plottypes, true );
      fprintf(' Finished simulating lateral network and generating tuning curves v2\n');
   catch ME
      fprintf(' Error simulating lateral network and generating tuning curves v2\n');
      getCatchMEstring(ME);
   end
end

if doL_v1
   try
      fprintf('\n\nStart generating lateral tuning curves v1\n');
      cd /Users/cedavey/Dropbox/code/visual_cortex/results/lateral/catDistribution/randomradius
      % fnameL    = '/Users/cedavey/Dropbox/code/visual_cortex/results/lateral/catDistribution/randomradius/outstruct_AB_BC_BL_LC_sz100.mat';
      tuningL_  = processVisualCortexInput( outstructL_, false, [], optargsL );
      curves    = generateVisualCortexTuningCurves( tuningL_, lateral,  responsefn, iterations, true, plottypes, true );
      fprintf(' Finished generating lateral tuning curves v1\n');
   catch ME
      fprintf(' Error generating lateral tuning curves v1\n');
      getCatchMEstring(ME);
   end
end

if doFF_v1
   try
      fprintf('\n\nStart generating feedforward tuning curves v1\n');
      cd /Users/cedavey/Dropbox/code/visual_cortex/results/feedforward/catDistribution/randomradius
      % fnameFF   = '/Users/cedavey/Dropbox/code/visual_cortex/results/feedforward/catDistribution/randomradius/outstruct_AB_BC_sz100.mat';
      tuningFF_ = processVisualCortexInput( outstructF_, false, [], optargsFF );
      curves    = generateVisualCortexTuningCurves( tuningFF_, feedforward, responsefn, iterations, true, plottypes, true );
      fprintf(' Finished generating feed forward tuning curves v1\n');
   catch ME
      fprintf(' Error generating feed forward tuning curves v1\n');
      getCatchMEstring(ME);
   end
end

cd( curr_dir );
