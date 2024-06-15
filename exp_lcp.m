function exp_lcp(number,n,n_trials)

save_results = 1; 
rng('default');
% number = 2;
% n = 300; 
% 
% n_trials = 2;%no. of random problems
    if (number ~= 3) || (nargin<3)
        n_trials = 1;
    elseif (number == 3) && (nargin<3)
        n_trials = 10;
    end
%Choose methods
pdmc  = 1; 
fb    = 1;
ps    = 1;
bpa   = 1*(number~=2);
ega   = 1;

%Parameters of PDMC/FB/PS algorithms
Qmat   = [0]; %1 does not converge
acclrt = [0 1]; 
idntfy = [0 1];
n_mthds = length(Qmat)*length(acclrt)*length(idntfy)*(pdmc+fb+ps)+bpa+ega;

%Other (default) parameters
LQ_coeff   = [1,0.999,1 1 1];     
sigma      = 1e-2;
eps        = 1e-6;
% stable_mins= [100,20,50,Inf];
iter_max   = 1e4; 


%% Solve using selected methods
index = method_number(pdmc,fb,ps,bpa,ega);

graphs = 1;
err = []; iters = []; cpu_time = []; methods = {};
times = {}; res = {}; objctv = {}; legends = {};  
identify_counters = [];
cvx_identifieds = [];
cpu_time_ls = [];
    for trial = 1:n_trials
        trial
        w0 = zeros(2*n,1);
        [M,d,~] = lcp(number,n);
              
        for mthdnmbr=index
            
            LQ = LQ_coeff(mthdnmbr);
%             stable_min = stable_mins(mthdnmbr);
            
            if mthdnmbr == 4 || mthdnmbr == 5
                Qmat_val = 0; idntfy_val = 0; acclrt_val = 0;
            else
                Qmat_val = Qmat; idntfy_val = idntfy; acclrt_val = acclrt;
            end
            
            for Q=Qmat_val
%                 stable_min = 10*(1-Q) + 10*Q; %good for lcp1
                stable_min = 50*(1-Q) + 100*Q; 
                for identify=idntfy_val
                    for acce=acclrt_val
                        stbl_min = acce*stable_min/2 + (1-acce)*stable_min;
                        [xsol,residuals,iter,obj,time,Times,method_label,identify_counter,cvx_identified,T_ls] = lcp_solver(mthdnmbr,M,d,w0,Q,LQ,iter_max,eps,acce,sigma,identify,stbl_min); 
  
                        %Store results
                        err(end+1)       = residuals(end);
                        iters(end+1)     = iter;
                        cpu_time(end+1)  = time;
                        methods{end+1}   = method_label;
                        identify_counters(end+1) = identify_counter;
                        cvx_identifieds(end+1) = cvx_identified;
                        cpu_time_ls(end+1)  = T_ls;
                        
                        times{end+1}    = Times;
                        res{end+1}      = residuals;
                        objctv{end+1}   = obj;
                        
                        %Graph of first instance
                        if  graphs && trial==1 %&& residuals(end)<eps

                            legends{end+1}  = method_label;
                            
                            residuals(end) = max(eps,residuals(end));
                            
                            figure(1);
                            semilogy(Times,residuals,'LineWidth',1.1); hold on;
                            
                            figure(2);
                            semilogy(residuals,'LineWidth',1.1); hold on;
                            

                        end
                    end
                end
            end
        end
    end

if graphs
    figure(1); 
    legend(legends);
    xlabel('CPU Time');
    ylabel('Residuals');
    
    figure(2); 
    legend(legends);
    xlabel('Iteration');
    ylabel('Residuals');
end

%Summary of Results
Method           = methods(1:n_mthds)';
MeanError       = mean(reshape(err',n_mthds,n_trials),2);       %Residuals
MeanIterations  = mean(reshape(iters',n_mthds,n_trials),2);
MeanRunningTime = mean(reshape(cpu_time',n_mthds,n_trials),2);
StdError        = std(reshape(err',n_mthds,n_trials),0,2);
StdIter     = std(reshape(iters',n_mthds,n_trials),0,2);
StdTime     = std(reshape(cpu_time',n_mthds,n_trials),0,2);
MeanRunningTimeLS = mean(reshape(cpu_time_ls',n_mthds,n_trials),2);
StdTimeLS     = std(reshape(cpu_time_ls',n_mthds,n_trials),0,2);

MeanIdentify  = mean(reshape(identify_counters',n_mthds,n_trials),2);
StdIdentify   = std(reshape(identify_counters',n_mthds,n_trials),0,2);
TotalConvexSetIdentified = sum(reshape(cvx_identifieds',n_mthds,n_trials),2);

Results = table(Method,MeanError,StdError,MeanIterations,StdIter,MeanIdentify,StdIdentify,MeanRunningTime,StdTime,MeanRunningTimeLS,StdTimeLS)

%For saving results
%Table
summary.Method          = Method;
summary.MeanError       = MeanError;
summary.MeanIterations  = MeanIterations;
summary.MeanRunningTime = MeanRunningTime;
summary.StdError        = StdError;
summary.StdIter         = StdIter;
summary.StdTime         = StdTime;
summary.MeanRunningTimeLS = MeanRunningTimeLS;
summary.StdTimeLS         = StdTimeLS;

summary.MeanIdentify    = MeanIdentify;
summary.StdIdentify     = StdIdentify;
summary.TotalConvexSetIdentified = TotalConvexSetIdentified;

summary.Results         = Results;
summary.Error       = reshape(err',n_mthds,n_trials);
summary.Iterations  = reshape(iters',n_mthds,n_trials);
summary.RunningTime = reshape(cpu_time',n_mthds,n_trials);

summary.residuals   = res; 
summary.objctv      = objctv;
summary.legends     = legends; 
summary.times       = times;

%Save
    if save_results
        filename = sprintf('lcp%d_n%d',number,n);
        save(filename,'summary');
    end
    

end    

%Assign index to solvers
function index= method_number(pdmc,fb,ps,bpa,ega)
index = [];
    if pdmc == 1
        index = 1;
    end
    
    if fb == 1
        index(end+1) = 2;
    end
    
    if ps == 1
        index(end+1) = 3;
    end
    
    if bpa == 1
        index(end+1) = 4;
    end
    
    if ega == 1
        index(end+1) = 5;
    end
    
end


%% Test problems

function [M,d,xsol] = lcp(number,n)

xsol = [];
          
    switch number
        case 1 
            M = diag(4*ones(1,n)) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1);
            d = ones(n,1);
            
        case 2 %P-matrix but not PD; MCPLIB
            M = 2*triu(ones(n),1) + eye(n);
            d = ones(n,1);
            xsol = zeros(n,1);
            xsol(n) = 1;
            
        case 3
            A = -5+10*rand(n);
            B = triu(-5+10*rand(n),1);
            B = B-B';
            eta = 0.3*rand(n,1);
            M = A*A' + B+ diag(eta);
            d = -500+1000*rand(n,1);
    end

end

%References
%1: Example 2 of 
%Qi, Liqun and Sun, Defeng and Zhou, Guanglu, A new look at smoothing
%Newton method for nonlinear complementarity problems and box constrained
%variational inequalities, Mathematical Programming, 87, 1--35, 2000.

%2 and 3: Example 7.1, 7.3 of 
%Kanzow, Christian, Some noninterior continuation methods for linear
%complementarity problems, SIAM Journal on Matrix Analysis and
%Applications, 17 (4), 851--868, 1996.