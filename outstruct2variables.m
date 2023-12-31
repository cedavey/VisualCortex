% Convert visual cortex output structure to individual variables - enables
% generating tuning curves

[outQ, outmu, outrates, outtime, outweights, ...
 layerconfig, ntwkconfig] ...
     = struct2v(outstruct,'outQ', 'outmu', 'outrates', 'outtime', 'outweights', ...
                          'layerconfig', 'ntwkconfig');

[T,dt,continueSim,...
 layernames, layerconns, inlayer, retina,...
 sz,spatialJitter,N,r,conndist,lambda,noisedist,noiseparams,RF,...
 wgtdist, wgtparams,b,beta,delta,...
 psptype,pspparams,pspdelay,psp,...
 plastic,deldist,delparams,layerdist,...
 ka,kb,n,k1,k2,Ra,Rb,endpoints,...
 runSim,evolveSep,estCov,calcCorr] = struct2v(layerconfig,...
        'T','dt','continueSim',...
        'layernames','layerconns','inlayer','retina',...
        'sz','spatialJitter','N','r','conndist','lambda',...
        'noisedist','noiseparams','RF',...
        'wgtdist','wgtparams','b','beta','delta',...
        'psptype','pspparams','pspdelay','psp',...
        'plastic','deldist','delparams','layerdist',...
        'ka','kb','n','k1','k2','Ra','Rb','endpoints',...
        'runSim','evolveSep','estCov','calcCorr');
 % continueSim is not in old code so will be returned empty if not there
 if isempty(continueSim), continueSim = false; end

if ~isempty(ntwkconfig)
   [cell_loc,conns,weights,delayinds,...
    time,interval,intervalQ,initTime,...
    totaldelay,maxdelay,synloc] ...
              = struct2v(ntwkconfig,'cell_loc','conns','weights','delayinds',...
                                    'time','interval','intervalQ','initTime',...
                                    'totaldelay','maxdelay','synloc');
end
