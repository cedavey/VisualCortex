% plasticLayers = getPlasticLayers( tuning )
%  OR
% plasticlayers = getPlasticLayers( outstruct )
%  OR
% plasticlayers = getPlasticLayers( layerconfig )
%
%    Returns a list of layers that have postsynaptic cells in a plastic 
%    conn, e.g. for outstruct.layerconfig.layerconns of {AB, BC, BL, LC},
%    with plastic connections between AB and BC, getPlasticLayers returns a
%    2x1 character array of ['B'; 'C'] to show that connections to layers B
%    and C are plastic.
function plasticLayers = getPlasticLayers( input )
   isplastic     = [];
   plasticLayers = [];

   % user passed in a tuning struct (output of processVisualCortexInput)
   if isfield( input, 'iterations' )
      if ~isfield( input, 'isplastic' )
         cprintf('Error',  sprintf('\nstruct error: got to add isplastic field to tuning struct\n\n'));
         return;
      end
      isplastic = input.isplastic;
      
   % user passed in an outstruct struct (output of visual_cortex)
   elseif isfield( input, 'layerconfig' )
      if ~isfield( input.layerconfig, 'plastic' )
         cprintf('Error',  sprintf('\nstruct error: got to add plastic field to outstruct layerconfig \n\n'));
         return;
      end
      isplastic = input.layerconfig.plastic;
      
   % user passed in a layerconfig struct (output of parseVisualCortexConfigFile)
   elseif isfield( input, 'layerconns' )
      if ~isfield( input, 'plastic' )
         cprintf('Error',  sprintf('\nstruct error: got to add plastic field to layerconfig struct\n\n'));
         return;
      end
      isplastic = input.plastic;
      
   else
      % perhaps user's input the isplastic boolean array from outstruct.layerconfig.plastic
      try
         % e.g  outstruct.layerconfig.plastic.AB = 0 -> input.AB = 0
         %      outstruct.layerconfig.plastic.BC = 1 -> input.BC = 1
         tmp = cell2mat( struct2cell( input ) );
         if islogical( tmp )
            isplastic = input;
         end
      catch
         isplastic = [];
      end
   end
      
   if isempty( isplastic ), return; end
   cxnames = fieldnames( isplastic ); % get all plastic connection names
   plastic = struct2array( isplastic );
   cxnames(~plastic) = [];
   plasticLayers = cell2mat( cxnames ); % convert from cell to string array 
   plasticLayers = unique( plasticLayers(:,2) ); % get unique second letters, since they're postsynaptic
end