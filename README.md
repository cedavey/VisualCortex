All files required for the simulation have been uploaded. 

The main simulation is run with visual_cortex.m. However, before running the simulation, network and neuron parameters must be set. Example parameter configuration files are  visualCortexFeedForwardConfig.m, and visualCortexWithInhibConfig.m, which parameterises a network with only feed forward connections, and a network with an additional inhibitory neuron population, respectively. Running a config file generates a network parameter structure. This can either be saved, or passed directly into the visual_cortex.m file, which runs the simulation. 

Once the simulation is completed, processVisualCortexInput can be run with a variety of input sinusoids, to determine the systems response to input with various spatial and temporal frequency characteristics. Finally, tuning curves can be plotted using generateVisualCortexTuningCurves.m
