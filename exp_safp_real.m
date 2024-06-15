%% Experiments on sparse affine feasibility problems with real-life datasets

function exp2(filename,save_results)
%filename: enter the filename of the .mat file containing the data

load(filename,'A','b');
w0 = A'*b;
[m,n]=size(A);
% w0 = zeros(n,1);

s_vals = round(0.05*n); %round(n*[0.005 0.01 0.05]);


%Choose methods
pdmc  = 1; 
fb    = 1;
ps    = 1;
pg_bt = 1;

%Parameters of PDMC/FB/PS algorithms
Qmat   = [0 1];
acclrt = [0 1]; 
idntfy = [1];
n_mthds = length(Qmat)*length(acclrt)*length(idntfy)*(pdmc+fb+ps)+pg_bt;

%Other (default) parameters
LQ_coeff   = [1 0.999 0.999 1.1];     
sigma      = 1e-1;
eps        = 1e-6;
% stable_min = 1e2;
iter_max   = 1e4; 


%% Solve using selected methods
index = method_number(pdmc,fb,ps,pg_bt);
graphs = 1;

    for j = 1:length(s_vals)
        s = s_vals(j)
       
        err = []; iters = []; cpu_time = []; methods = {};
        times = {}; objctv = {}; legends = {};  
        identify_counters = [];
        cvx_identifieds = [];
              
        for mthdnmbr=index
            LQ = LQ_coeff(mthdnmbr);
            
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
                        [~,iter,obj,time,Times,method,identify_counter,cvx_identified] = saf_solver(mthdnmbr,A,b,s,w0,Q,LQ,iter_max,eps,acce,sigma,identify,stbl_min);
                        
                        err(end+1)      = obj(end);
                        iters(end+1)    = iter;
                        cpu_time(end+1) = time;
                        methods{end+1}  = method;
                        identify_counters(end+1) = identify_counter;
                        cvx_identifieds(end+1) = cvx_identified;
                        
                        %Graph
                        if graphs %&& trial==1 %&& obj(end)<eps
                            times{end+1}    = Times;
                            objctv{end+1}   = obj;
                            legends{end+1}  = method;
                            
                            obj(end) = max(eps,obj(end));

                            
                            figure(2*j-1); 
                            semilogy(Times,obj,'LineWidth',1.1); hold on;
                            
                            figure(2*j);
                            semilogy(obj,'LineWidth',1.1); hold on;
                            

                        end
                    end
                end
            end
        end

        
        if graphs
            figure(2*j-1); 
            legend(legends);
            xlabel('CPU Time');
            ylabel('0.5||Aw^k-b||^2+0.5dist(w^k,S_2)^2');

            figure(2*j); 
            legend(legends);
            xlabel('Iteration');
            ylabel('0.5||Aw^k-b||^2+0.5dist(w^k,S_2)^2');
        end

        %Summary of Results
        Method      = methods';
        Error       = err';
        Iterations  = iters';
        LS_Iter    = identify_counters';
        Cvx_Identified = cvx_identifieds';
        RunningTime = cpu_time';

        Results = table(Method,Error,Iterations,LS_Iter,Cvx_Identified,RunningTime)

        %For saving results
        %Table
        summary.Method      = Method;
        summary.Error       = Error;
        summary.Iterations  = Iterations;
        summary.Identify    = LS_Iter;
        summary.RunningTime = RunningTime;
        summary.Cvx_Identified = Cvx_Identified;
        summary.Results     = Results;

        summary.objctv      = objctv;
        summary.legends     = legends; 
        summary.times       = times;

        %Save
            if save_results
                fname = sprintf('saf_real_%s_s%d',erase(filename,'.mat'),s);
                save(fname,'summary');
            end
    
    end
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