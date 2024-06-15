%% Experiments on randomly generated sparse affine feasibility problems

save_results = 1; 
rng('default');
n = 10000; m = round(n/4); s = round(m/4); 
alpha = 5; %Parameter for test problem (see "saf" at the bottom part)

n_trials = 10;%no. of random problems

%Choose methods
pdmc  = 1; 
fb    = 1;
ps    = 1;
pg_bt = 1;

%Parameters of PDMC/FB/PS algorithms
Qmat   = [0 1];
acclrt = [0 1]; 
idntfy = [0 1];
n_mthds = length(Qmat)*length(acclrt)*length(idntfy)*(pdmc+fb+ps)+pg_bt;

%Other (default) parameters
LQ_coeff   = [1,0.999,0.999 1.1];     
sigma      = 1e-2;
eps        = 1e-6;
% stable_mins= [100,20,50,Inf];
iter_max   = 1e4; 


%% Solve using selected methods
index = method_number(pdmc,fb,ps,pg_bt);

graphs = 1;
err = []; iters = []; cpu_time = []; methods = {};
times = {}; objctv = {}; legends = {};  
identify_counters = [];
cvx_identifieds = [];
cpu_time_ls = [];
    for trial = 1:n_trials
        trial
        [A,b,~] = saf(m,n,s,alpha);
        w0 = A'*b; 
              
        for mthdnmbr=index
            LQ = LQ_coeff(mthdnmbr);
%             stable_min = stable_mins(mthdnmbr);
            
            if mthdnmbr == 4
                Qmat_val = 0; idntfy_val = 0; acclrt_val = 0;
            else
                Qmat_val = Qmat; idntfy_val = idntfy; acclrt_val = acclrt;
            end
            
            for Q=Qmat_val
                stable_min = 50*(1-Q) + 100*Q;
                for identify=idntfy_val
                    for acce=acclrt_val
                        stbl_min = acce*stable_min/2 + (1-acce)*stable_min;
                        [~,iter,obj,time,Times,method,identify_counter,cvx_identified,T_ls] = saf_solver(mthdnmbr,A,b,s,w0,Q,LQ,iter_max,eps,acce,sigma,identify,stbl_min);
                        
                        err(end+1)      = obj(end);
                        iters(end+1)    = iter;
                        cpu_time(end+1) = time;
                        methods{end+1}  = method;
                        identify_counters(end+1) = identify_counter;
                        cvx_identifieds(end+1) = cvx_identified;
                        cpu_time_ls(end+1) = T_ls;
                        
                        times{end+1}    = Times;
                        objctv{end+1}   = obj;
                        
                        %Graph of first instance 
                        if graphs && trial==1 %&& obj(end)<eps                    
                            obj(end) = max(eps,obj(end));
                            legends{end+1}  = method;
                            
                            figure(1); 
                            semilogy(Times,obj,'LineWidth',1.1); hold on;
                            
                            figure(2);
                            semilogy(obj,'LineWidth',1.1); hold on;
                            

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
    ylabel('0.5||Aw^k-b||^2+0.5dist(w^k,S_2)^2');
    
    figure(2); 
    legend(legends);
    xlabel('Iteration');
    ylabel('0.5||Aw^k-b||^2+0.5dist(w^k,S_2)^2');
end

%Summary of Results
Method           = methods(1:n_mthds)';
MeanError       = mean(reshape(err',n_mthds,n_trials),2);
MeanIterations  = mean(reshape(iters',n_mthds,n_trials),2);
MeanRunningTime = mean(reshape(cpu_time',n_mthds,n_trials),2);
StdError        = std(reshape(err',n_mthds,n_trials),0,2);
StdIter     = std(reshape(iters',n_mthds,n_trials),0,2);
StdTime     = std(reshape(cpu_time',n_mthds,n_trials),0,2);

MeanIdentify  = mean(reshape(identify_counters',n_mthds,n_trials),2);
StdIdentify   = std(reshape(identify_counters',n_mthds,n_trials),0,2);
TotalConvexSetIdentified = sum(reshape(cvx_identifieds',n_mthds,n_trials),2);
MeanRunningTimeLS = mean(reshape(cpu_time_ls',n_mthds,n_trials),2);
StdTimeLS     = std(reshape(cpu_time_ls',n_mthds,n_trials),0,2);

 
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

summary.MeanIdentify    = MeanIdentify;
summary.StdIdentify     = StdIdentify;
summary.TotalConvexSetIdentified = TotalConvexSetIdentified;
summary.MeanRunningTimeLS = MeanRunningTimeLS;
summary.StdTimeLS         = StdTimeLS;

summary.Results         = Results;
summary.Error       = reshape(err',n_mthds,n_trials);
summary.Iterations  = reshape(iters',n_mthds,n_trials);
summary.RunningTime = reshape(cpu_time',n_mthds,n_trials);

summary.objctv      = objctv;
summary.legends     = legends; 
summary.times       = times;

%Save
    if save_results
        filename = sprintf('saf_random_n%d_m%d_s%d_alpha%d',n,m,s,alpha);
        save(filename,'summary');
    end
    

%Test problem
function [A,b,x_sol] = saf(m,n,s,alpha)
A = randn(m,n);

ind = randsample(n,s);
x_sol = zeros(n,1); 
    switch alpha
        case 0
            x_sol(ind) = randn(s,1);
        otherwise
            x_sol(ind) = sign(2*rand(s,1)-1).*10.^(alpha*rand(s,1));
            %Considered in NESTA paper; "high dynamic range"
    end
    
b = A*x_sol; 
    
    end

    
%Assign index to solvers
function index= method_number(pdmc,fb,ps,pg_bt)
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
    
    if pg_bt == 1
        index(end+1) = 4;
    end
end