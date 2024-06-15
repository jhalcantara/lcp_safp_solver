function  [xsol,residuals,iter,objs,time,Times,method_label,identify_counter,cvx_identified,T_ls] = lcp_solver(method_index,M,d,w0,Q,LQ_coeff,iter_max,eps,acce,sigma,identify,stable_min)
    switch method_index
        case 1
            method = 'PDMC';
            [xsol,residuals,iter,objs,time,Times,method_label,identify_counter,cvx_identified,T_ls] = pm_lcp(method,M,d,w0,Q,LQ_coeff,iter_max,eps,acce,sigma,identify,stable_min);
        case 2
             method = 'FB';
            [xsol,residuals,iter,objs,time,Times,method_label,identify_counter,cvx_identified,T_ls] = pm_lcp(method,M,d,w0,Q,LQ_coeff,iter_max,eps,acce,sigma,identify,stable_min);
        case 3
             method = 'PS';
            [xsol,residuals,iter,objs,time,Times,method_label,identify_counter,cvx_identified,T_ls] = pm_lcp(method,M,d,w0,Q,LQ_coeff,iter_max,eps,acce,sigma,identify,stable_min);
        case 4
            x0 = w0(1:length(d));
            identify_counter = Inf;
            cvx_identified = 0;
            T_ls = 0;
            [xsol,residuals,iter,time,Times,method_label] = bpa_lcp(M,d,x0,iter_max,eps); objs = [];
        case 5
            x0 = w0(1:length(d));
            identify_counter = Inf;
            cvx_identified = 0;
            T_ls = 0;
            [xsol,residuals,iter,time,Times,method_label] = ega_lcp(M,d,x0,iter_max,eps); objs = [];

    end
            
end

%% PDMC, FB, and PS
function [xsol,residuals,iter,objs,time,T,method_label,identify_counter,cvx_identified,T_ls] = pm_lcp(method,M,d,w0,Q,LQ_coeff,iter_max,eps,acce,sigma,identify,stable_min)
start = tic; 
%Parameters
n = length(d); 
  
    switch Q
        case 1 %Q=I_n
            S = @(x) M*(M'*x)+x; 
            dS = 1;
            LQ = eigs(S,n,1,'largestreal','Tolerance',1e-3); %norm(S)
            
        otherwise
            f = norm(M,1)/sqrt(n);
            M = M./f;
            d = d./f;
            S = M*M'+eye(n); 
            dS = decomposition(S,'chol');
            LQ = 1;
    end

lambda = LQ_coeff/LQ;
    
%Initialization
iter = 1;
unchanged = 0; 
identify_counter = 0;
cvx_identified = 0;    
T_ls = 0;

%First iteration
[PS2w0,I0] = PS2(w0);
    
    if strcmp(method,'PS')
        w0 = PS2w0;
    end
    
[nabla0,affine0,err] = nabla_fQ(w0,M,d,Q,dS,n);

residuals = err;
objs = 0.5*norm(affine0)^2 + 0.5*norm(w0-PS2w0)^2;
T = toc(start);

[w1,PS2w1,I1] = newpoint(method,w0,lambda,nabla0,PS2w0);
w0 = w1; 

%Main Loop    
    while (err>eps) && (iter<iter_max)
        iter = iter+1;
        
        [nabla1,affine1,err] = nabla_fQ(w1,M,d,Q,dS,n);
        
        residuals = [residuals; err]; 
        obj = 0.5*norm(affine1)^2 + 0.5*norm(w1-PS2w1)^2;
        objs = [objs; obj];
        
        T = [T; toc(start)];

        %
        if identify || acce
            stable = isequal(I0,I1);
            unchanged = stable*(unchanged +1);
            I0 = I1;
        end
        
        if acce && stable
            p = w1 - w0;
                
                if Q == 1
                    old_obj = obj;
                else
                    old_obj = 0.5*(norm(nabla1)^2+norm(w1-PS2w1)^2);
                end
                
            [z,PS2z,nabla] = linesearch(method,Q,sigma,w1,p,old_obj,affine0,nabla0,affine1,nabla1);
  
            [w,PS2w1,I1] = newpoint(method,z,lambda,nabla,PS2z);
            
            %Update new values
            w0 = w1; 
            w1 = w ;  
        
        else%non-accelerated version
            w0 = w1; 
            [w1,PS2w1,I1] = newpoint(method,w1,lambda,nabla1,PS2w1);
        end    
        
        nabla0 = nabla1; 
        affine0 = affine1;
        

        %If convex set is identified
        if identify && (unchanged>=stable_min) %|| (err<1e-2))
            iter = iter+1;
            
            identify_counter = identify_counter+1;
            
            LS_start = tic; 
            w = lin_solve(M,d,n,I1);
            T_ls = T_ls + toc(LS_start);
            
            err = norm(min(w(1:n),M*w(1:n)-d));
            
            if err<eps 
                residuals = [residuals; err]; 
                objs = [objs; 0.5*norm(w-PS2(w))^2];    %w1 lies on S1
                T = [T; toc(start)];
                w1 = w;
                cvx_identified = 1;
                break; 
            else 
                residuals = [residuals; residuals(end)];
                objs = [objs; objs(end)]; 
                T = [T; toc(start)];
                unchanged = -1; %restart
                if mod(identify_counter,10)==0 || identify_counter == 1
                    fprintf('     Incorrectly identified convex set, identify_counter = %d \n',identify_counter);
                end
            end
            
            
            
        end    
    end%end while

xsol = w1(1:n);
time = toc(start);
  
%Results
method_label = methodname(method,Q,acce,identify);

    if err<eps
        if ~identify
            fprintf('Solution obtained via %s \n',method_label);
        else
            fprintf('Solution obtained via %s ; ',method_label);
            
            if cvx_identified
                fprintf('Convex set identified; ');
            end
            
        end
    else
        if ~identify
            fprintf('Rerun %s \n',method_label);
        else
            fprintf('Rerun %s ; ',method_label);
        end
    end     
    
    if identify
        fprintf('Total number of incorrect identification = %d \n',identify_counter-cvx_identified);  
    end
end


%% Formulas
%Compute new iterate
function [w1,PS2w1,I1] = newpoint(method,w0,lambda,nabla0,PS2w0)
    switch method
        case 'PDMC'
            w1 = (w0-lambda*nabla0+lambda*PS2w0)/(1+lambda);
            [PS2w1,I1] = PS2(w1);
            
        case 'FB'
            f_step = w0-lambda*nabla0; 
            [PS2f_step,~] = PS2(f_step);
            w1 = (lambda*PS2f_step + f_step)/(1+lambda); 
            [PS2w1,I1] = PS2(w1);
            
        case 'PS'
            [w1,I1] = PS2(w0-lambda*nabla0); 
            PS2w1 = w1; 
    end
    
end

%Line search
function [z,PS2z,nabla] = linesearch(method,Q,sigma,w1,p,old_obj,affine0,nabla0,affine1,nabla1)
norm_p = norm(p)^2;
    switch method
        case {'PDMC','FB'}
            t = 10; diff = 1;
                switch Q
                    case 1
                        dif      = affine1-affine0;
                        norm_obj = 0.5*norm(affine1)^2;
                        inn      = affine1'*dif;
                    otherwise
                        dif      = nabla1-nabla0;
                        norm_obj = 0.5*norm(nabla1)^2;
                        inn      = nabla1'*dif;
                end
                
                norm_dif = 0.5*norm(dif)^2;
                
                while (diff>0) && (t>1e-12)
                        t = 0.5*t;
                        %Denote f_Q(new) = 0.5*norm(new)^2;
                        %find "new"
                        z = w1 + t*p;
                        [PS2z,~] = PS2(z);
                        diff = (norm_obj + t*inn + t^2*norm_dif) + 0.5*norm(z-PS2z)^2  - old_obj + sigma*t^2*norm_p/2;
                end

        case 'PS'
            if any(p<0)
                t1 = min(-w1(p<0)./p(p<0)); 
        	else
                t1 = Inf;
            end
                    
            descent = nabla1'*p;

            if t1~=0
                %Exact
                switch Q
                    case 1
                        t2 = max(-2*descent,0)/(norm(affine1-affine0)^2+sigma*norm_p);
                    otherwise
                        t2 = max(-2*descent,0)/(norm(nabla1-nabla0)^2+sigma*norm_p);
                end
                t  = min(t1,t2);
            else
                t = 0;
            end   
            
            z = w1+t*p;
            PS2z = [];
    end
    
nabla = (1+t)*nabla1 - t*nabla0;
        
end

%Gradient of f_Q
function [nabla,affine,err] = nabla_fQ(w,M,d,Q,dS,n)
u = w(1:n); v = w(n+1:2*n);

%Compute Tw-d
affine = M*u-v-d;       

%Compute nabla f_Q = T'*Q*(Tw-d)
    if Q==1
        Qobj = affine;
    else
        Qobj = dS\affine;
    end
    
nabla = [M'*Qobj; -Qobj];

%Compute err
err = norm(min(u,affine+v));
end

%Projection onto S2
function [z,I] = PS2(w)
N = length(w);
z = reshape(w,N/2,2);
% a = z + rand(N/2,2).*(abs(z-z(:,[2 1]))<1e-12);
% I = reshape(a>a(:,[2 1]),N,1);
E = z==z(:,[2,1]);
E = E(:,1);
I = reshape(z>z(:,[2 1]),N,1);
I(E) = 1;
z = max(w,0).*I;
end

%Linear solver after convex set identification
function w = lin_solve(M,d,n,I1)
T = [M -eye(n)]; 
U = T(:,I1~=0);

nzero = U\d; 
w = zeros(2*n,1);
w(I1~=0) = nzero;
end


%% Two traditional projection methods for variational inequalities 

%% Basic Projection Algorithm (Algorithm 12.1.1, Facchinei and Pang)
%For positive definite A only (not necessarily symmetric)
function [xsol,residuals,iter,time,T,method_label] = bpa_lcp(M,d,x0,iter_max,eps)
start = tic;

%Parameters
S = 0.5*(M+M');
mu = min(eig(S)); 

%Initialization
err = 1;
iter = 0;
residuals = [];
T = [];

%Main Loop
    if mu>0
        
        norm_M = svds(M,1,'largest','Tolerance',1e-3);
        tau = 1.001*(0.5*norm_M^2/mu);
        
        while err>eps && iter<iter_max
            iter = iter+1; 
            y0 = M*x0-d;
            
            err = norm(min(x0,y0));
            residuals = [residuals; err];
            T = [T; toc(start)];
            
            x0 = max(x0-(1/tau)*y0,0);
        end
        xsol = x0;
    else
        residuals = Inf; 
        fprintf('Method not applicable \n');
        xsol =[];
    end
    
time = toc(start);
method_label = 'BPA';
    if err<eps
        fprintf('Solution obtained using %s \n',method_label);
    else
        fprintf('Rerun %s \n',method_label);
    end  
end

%% Extragradient algorithm (Algorithm 12.1.9, Facchinei and Pang)
function [xsol,residuals,iter,time,T,method_label] = ega_lcp(M,d,x0,iter_max,eps)
start = tic;

%Parameters
norm_M = svds(M,1,'largest','Tolerance',1e-3); %norm(M);
tau = 0.999/norm_M;

%Initialization
err = 1;
iter = 0;
T = [];
residuals = [];

%Main Loop
    while err>eps && iter<iter_max
        iter = iter+1; 
        y0 = M*x0-d;
        
        err = norm(min(x0,y0));
        residuals = [residuals; err];
        T = [T; toc(start)];
        
        x0_half = max(x0-tau*y0,0);
        x0 = max(x0-tau*(M*x0_half-d),0);
    end    
xsol = x0;
time = toc(start);

method_label = 'EGA';

    if err<eps
        fprintf('Solution obtained using %s \n',method_label);
    else
        fprintf('Rerun %s \n',method_label);
    end  
    
end


%% Miscellaneous
function method_label = methodname(method,Q,acce,identify)
method_label = method;
    if Q==0
        switch method
            case 'PDMC'
                method_label = 'MAveP';
            case 'FB'
                method_label = 'MARP';
            case 'PS'
                method_label = 'MAP';
        end
    end
    
    if acce
        method_label = strcat('A',method_label);
    end
    
    if identify
        method_label = strcat(method_label,'+');
    end
end





