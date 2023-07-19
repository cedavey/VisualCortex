% [mu_anal, mu_emp] = getVisualCortexMeanWeight(outstruct)
% Compares analytic & empirical mean weights
%

% TO DO: make it so returns analytic only if no empirical data found
%   - input config file name or layerconfig
function [mu_anal, mu_emp] = getVisualCortexMeanWeight(out)
   % if input an output struct
   if isfield(out,'layerconfig') && isfield(out,'ntwkconfig')
      conns = fieldnames(out.layerconfig.plastic); 
      for ci=1:length(conns)
         if out.layerconfig.plastic.(conns{ci})
            cxlabel = conns{ci};
            mu_emp.(cxlabel)  = mean(cellfun(@(w) mean(w(end,:)), out.outweights.(cxlabel))); 

            % calculate analytic mean weight
            k1 = out.layerconfig.k1.(cxlabel);
            k2 = out.layerconfig.k2.(cxlabel);
            N  = out.layerconfig.N.(cxlabel);
            Rb = out.layerconfig.Rb.(cxlabel);
            m  = -N*Rb*k1/k2;
            % no co-inputs into postsynaptic layer
            mu_anal.(cxlabel) = m;
            if ~isempty(out.ntwkconfig.co_inputs.(cxlabel))
               % have co-inputs into postsynaptic layer
               for li=1:length(out.ntwkconfig.co_inputs.(cxlabel))
                  colabel    = out.ntwkconfig.co_inputs.(cxlabel){li}.colabel;
                  otherinput = out.ntwkconfig.co_inputs.(cxlabel){li}.other_input;
                  co_k2      = out.ntwkconfig.co_inputs.(cxlabel){li}.k2;
                  % TO DO: adjust this calc of inhib wgt wrt distribution
                  co_wgt     = mean(out.layerconfig.wgtparams.(otherinput));
                  co_Rb      = out.layerconfig.Rb.(otherinput);
                  co_N       = out.layerconfig.N.(otherinput);
                  co_m(li)   = co_k2*co_N*co_Rb*co_wgt/k2;
                  
                  co_Qav     = mean(toVec(out.co_outputs.(cxlabel){li}.Q(end,:,:)));
                  co_m2(li)  = co_m(li) + co_N*co_Rb*co_wgt*co_Qav/k2;
%                   acovM(ni,pi) = (-k1_BC(ni,pi) - w_LC(ni,pi)*k2_LC(ni,pi) - w_LC(ni,pi)*Q_BL(ni,pi)) / k2_BC(ni,pi);

               end
               mu_anal.(cxlabel) = m - sum(co_m);
               mu_anal.inclQ.(cxlabel) = m - sum(co_m) - sum(co_m2);
            end
         end
      end
   end
end

