% Evaluate analytic & empirical mean as a parameter changes
% Model assumed to have A, B, C, L, where L is the inhibitory layer

% TO DO: impact of n.BC on mean wgt - do n=0.5 & n=1 have same mean wgt??

% With inhib, mean weight equation is:
% w_BC ~= (-k1_BC - w_LC*k2_LC) / k2_BC
% w_BC  = (-k1_BC - w_LC*k2_LC - w_LC*E(Q_BC)) / k2_BC
% k2_LC/k2_BC = nu_L / nu_B (where nu indicates mean rate)
runSim    = false;
calcAnal  = true;

if runSim
   nN     = 10; % number of repetitions for each parameter value
   param  = 'weight';
   config = visualCortexWithInhibConfig_Errol('inhib',false);

   switch param
      case 'weight'
         minp = -2; maxp = 0; deltap = 0.5;
         paramfield = 'config.wgtparams.LC';
      case 'rateL'
         minp = 10; maxp = config.lambda.B*2; deltap = 10;
         paramfield = 'config.L';
      case 'rateB'
         minp = 10; maxp = config.lambda.L*2; deltap = 10;
         paramfield = 'config.B';
      case 'rateRatio'
         rateB = 20;
         minp = 0.5*rateB; maxp = 1.75*rateB; deltap = 0.25*rateB;
         config.B = rateB;
         paramfield = 'config.L';
      otherwise
         cprintf('errors','\nHaven''t coded for this parameter...exiting\n\n');
         return;
   end
   pvals = minp:deltap:maxp;
   nP    = length(pvals);
   alloc = @() zeros(nN,nP);
   empM  = alloc();
   analM = alloc();
   acovM = alloc();
   k1_BC = alloc();
   k2_BC = alloc();
   k2_LC = alloc();
   w_LC  = alloc();
   Q_BL  = alloc();

   for pi=1:nP
   %    eval([paramfield '=' num2str(pvals(pi)) ';']);
      eval(sprintf('%s = %g;',paramfield,pvals(pi)));
      for ni=1:nN
         cprintf('comments',sprintf('Processing %s = %g (rep %d)\n',param,pvals(pi),ni));
         inhib        = visual_cortex(config, 0, 'inhib');
         k1_BC(ni,pi) = config.k1.BC;
         k2_BC(ni,pi) = config.k2.BC;
         k2_LC(ni,pi) = inhib.ntwkconfig.co_inputs.BC{1}.k2;
         w_LC(ni,pi)  = config.wgtparams.LC;
         Q_BL(ni,pi)  = mean(toVec(inhib.co_outputs.BC{1}.Q(end,:,:)));
         analM(ni,pi) = (-k1_BC(ni,pi) - w_LC(ni,pi)*k2_LC(ni,pi)) / k2_BC(ni,pi);
         acovM(ni,pi) = (-k1_BC(ni,pi) - w_LC(ni,pi)*k2_LC(ni,pi) - w_LC(ni,pi)*Q_BL(ni,pi)) / k2_BC(ni,pi);
         empM(ni,pi)  = mean(cellfun(@(c) mean(c(end,:)), inhib.outweights.BC));
      end
   end
end

% Allows re-calculation of analytical mean weights, just in case eqns change
% Note: assumes all values are set already - k1, k2, Q, etc
if calcAnal
   for pi=1:nP
   %    eval([paramfield '=' num2str(pvals(pi)) ';']);
      eval(sprintf('%s = %g;',paramfield,pvals(pi)));
      for ni=1:nN
         analM(ni,pi) = (-k1_BC(ni,pi) - w_LC(ni,pi)*k2_LC(ni,pi)) / k2_BC(ni,pi);
         acovM(ni,pi) = (-k1_BC(ni,pi) - w_LC(ni,pi)*k2_LC(ni,pi) - w_LC(ni,pi)*Q_BL(ni,pi)) / k2_BC(ni,pi);
      end
   end
end

figure; nr=2; nc=2; fi=1;
subplot(nr,nc,fi),fi=fi+1;
plot(pvals,empM');
xlabel('Repetition');
% legend(str2legend([param '='],pvals));
title('Empirical');

subplot(nr,nc,fi),fi=fi+1;
plot(pvals,analM(1,:));
xlabel(param);
title('Simple analytical');

subplot(nr,nc,fi),fi=fi+1;
plot(pvals,acovM');
title('Full analytical');
xlabel('Repetition');

subplot(nr,nc,fi),fi=fi+1;
plot(pvals,[mean(empM)' mean(analM)' mean(acovM)']);
legend('Empirical','Simple anal','Full anal');
xlabel(param);




