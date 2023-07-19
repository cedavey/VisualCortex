% co_layers = getJointInputLayers(pre,post,connlabels,r)
function co_layers = getJointInputLayers(layerconfig,pre,post,connlabels,r)
% function co_layers = getJointInputLayers(pre,post,connlabels,r,conns,delayinds)
   % TO DO: generalise calculation of effective radius for covariance as a
   % function of distance from input layer, & distance btwn joint inputs to
   % postsynaptic layer (see pg. 24 of Plasticity #3 notebook)
   
   % If multiple inputs to a layer, we need cov of all layers inputting
   % to the layer, even if the layer is not itself plastic. 
   
   % Note that if there are common inputs to a layer, we need covariance
   % between the 2 presynaptic layers. Then we need number of connections,
   % & k2 from other presynaptic layer to postsynaptic layer connection. So
   % there are now 2 additional connections of interest. 
   % E.g. if B->C & L->C, & one of these is plastic (say B->C), then we are 
   % interested in any connections btwn B & L (B->L or L->B) to determine 
   % covariance between inputs to layer C. We are also interested in the 
   % non-plastic connection to C (L->C in this case) as this will determine
   % the new value of k2 for the B->C weight udpate equation.
   % I've called these 'co_input' (B->L) & 'other_input' (L->C).
   cxlabel   = [pre post];
   co_layers = {};
   [sz,N,lambda,Rb,kb,k2] = struct2v(layerconfig,'sz','N','lambda','Rb','kb','k2'); 
   for cc=1:size(layerconfig.layerconns,1)
      % Check for joint layered inputs to postsynaptic layer
      % i.e. pre -> post & co_pre -> post
      if (connlabels(cc,1)~=pre) && (connlabels(cc,2)==post)
         co_pre = connlabels(cc,1); % another input layer to postsyn layer
         % 2 layers both inputting to the same postsynaptic layer will only
         % effect each other's weights if they're correlated themselves
         % Find out if they're connected, or if they have a common input
         % Connection possibilities I'm considering are:
         % [pre -> co_pre] & [co_pre -> pre]
         % TO DO: [input -> pre; input -> co_pre]
         
         % Note: I've assumed that whichever direction the connection is
         % that's causing correlation btwn pre & co_pre (i.e. co_pre -> pre
         % or pre -> co_pre), the presynaptic one is also the one getting
         % the common input from an earlier layer

         % case: co_pre -> pre with input -> co_pre
         if strmatch([co_pre pre], connlabels) 
            % E.g. B->C (plastic), B->L (not plastic), L->C (not plastic)
            jt_input    = co_pre;
            co_post     = pre;
            co_pre      = co_pre; 
            colabel     = [co_pre co_post];
%             pre2postCxs = true;
            r_post      = r.(colabel);
            % to calc radius need to know radius of connections
            % coming into pre, since that determines cov(pre,pre) component
            preinput    = connlabels(connlabels(:,2)==co_pre,:);
            % only calculate effective radius for a single pre-input, o/w
            % too complicated!
            if size(preinput,1) > 1
               preinput = preinput(1,:);
            end
            if isempty(preinput)
               r_pre    = 0;
            else
               r_pre    = r.(preinput);
            end
            r_post      = r.(colabel);
            r_eff       = r_pre + r_post/2; % effective radius for cov calc. 

         % case: pre -> co_pre with input -> pre
         elseif strmatch([pre co_pre], connlabels)
            % E.g. B->C (plastic), L->B (not plastic), L->C (not plastic)
            jt_input    = co_pre;
            co_post     = co_pre;
            co_pre      = pre;
            colabel     = [co_pre co_post];
%             pre2postCxs = false;
            % to calc radius need to know radius of connections
            % coming into pre, since that determines cov(pre,pre) component
            preinput = connlabels(connlabels(:,2)==pre,:);
            if isempty(preinput)
               r_pre    = 0;
            else
               r_pre    = r.(preinput);
            end
            r_post      = r.(colabel);
            r_eff       = r_pre + r_post/2; % effective radius for cov calc. 

         % case: input -> pre & input -> co_pre 
         elseif any(connlabels(:,2)==pre) && any(connlabels(:,2)==co_pre)
            % make sure input layer connecting to pre & co_pre is the same
            if any(connlabels(connlabels(:,2)==pre,1) == connlabels(connlabels(:,2)==co_pre,1))
               % Note: assumes pre layer is only postsyn layer once 
               % TO DO: make more general
               inlayer  = connlabels(connlabels(:,2)==pre,1);
               jt_input = co_pre;
               r_eff    = r.([inlayer pre]/2) + r.([inlayer co_pre]/2);
               co_post  = co_pre;
               co_pre   = inlayer;
%                pre2postCxs = false;
               % TO DO: need to calc cov. btwn 2 layers that are only
               % indirectly connected via common input. This doesn't fit
               % with the way I've written the code, which calculates
               % covariance btwn directly connected common input layers
               cprintf('comment','\n\tgetJointInputLayers: if using analytical covariance, note that this function is not yet coded for indirectly correlated joint inputs\n');
            else
               cprintf('comment','\n\tgetJointInputLayers: if using analytical covariance, note that this function is not yet coded for indirectly correlated joint inputs, setting joint cov to 0\n');
               co_post     = [];
            end
         else
            co_post     = [];
         end
         if ~isempty(co_post)
            colabel     = [co_pre co_post];
            if layerconfig.estCov % estimate from rates (keep track of current estimate)
               mu_in    = ones(prod(sz.(co_pre)),prod(sz.(co_post)))*lambda.(co_pre);
               mu_out   = ones(1,prod(sz.(co_post)))*lambda.(co_post);
               sig_in   = ones(prod(sz.(co_pre)),prod(sz.(co_post)))*lambda.(co_pre);
               sig_out  = ones(1,prod(sz.(co_post)))*lambda.(co_post);
               Q_co     = zeros(prod(sz.(co_pre)),prod(sz.(co_post)));
            else % determine analytically 
               mu_in    = [];
               mu_out   = [];
               sig_in   = [];
               sig_out  = [];
               Q_co     = setBetweenLayerCov(sz.(co_pre), sz.(co_post),r_eff);
%                [Q_co, conns_co, del_co] = ...
%                       setBetweenLayerCov(sz.(co_pre), sz.(co_post),...
%                                          reff, conns_, delay_, pre2postCxs);
            end
            % calculate k2 contribution to covariance
            % First find presynlayer background rate from eqn 2 in 
            % Linsker (F_0^M), & the definition of k2 in eqn 4
            % k_2 = lambda(lambda - bckgrnd)*N*R_b*k_b
            other_input = [jt_input post]; % other input to postsyn layer
            bckgrnd = lambda.(pre) - k2.(cxlabel)/(kb.(cxlabel)*lambda.(pre));
            k2_co   = lambda.(jt_input)*(lambda.(pre) - bckgrnd)*...
                      N.(other_input)*kb.(cxlabel)*Rb.(other_input);
%             k2_co = k2.(cxlabel)/lambda.(pre)*lambda.(jt_input);
k2_co = k2.([co_post post]);
            co_layer = v2struct(k2_co, Q_co, sig_in, sig_out, mu_in, mu_out, r_eff, colabel, other_input,...
                               {'fieldnames', 'k2', 'Q', 'sig_in', 'sig_out', 'mu_in', 'mu_out', 'r', 'colabel', 'other_input'});
%             co_layer = v2struct(Q_co,mu_co,sig_co,r_eff,conns_co,colabel,connlabels(cc,1),del_co,...
%                                 {'fieldnames','Q','mu','sig','r','conns','colabel','co_pre','delayinds');
            co_layers{end+1} = co_layer;
         end
      end
   end
end


