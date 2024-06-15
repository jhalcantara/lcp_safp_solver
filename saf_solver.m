function [w,iter,obj,time,T,method_label,identify_counter,cvx_identified,T_ls] = saf_solver(mthdnmbr,A,b,s,w0,Q,LQ_coeff,iter_max,eps,acce,sigma,identify,stable_min)
    switch mthdnmbr
        case 1
            method = 'PDMC';
            [w,iter,obj,time,T,method_label,identify_counter,cvx_identified,T_ls] = pm_saf(method,A,b,s,w0,Q,LQ_coeff,iter_max,eps,acce,sigma,identify,stable_min);
        case 2
            method = 'FB';
            [w,iter,obj,time,T,method_label,identify_counter,cvx_identified,T_ls] = pm_saf(method,A,b,s,w0,Q,LQ_coeff,iter_max,eps,acce,sigma,identify,stable_min);
        case 3
            method = 'PS';
            [w,iter,obj,time,T,method_label,identify_counter,cvx_identified,T_ls] = pm_saf(method,A,b,s,w0,Q,LQ_coeff,iter_max,eps,acce,sigma,identify,stable_min);
        case 4
            eta = LQ_coeff; 
            identify_counter = Inf; 
            cvx_identified = 0;
            T_ls = 0;
            [w,iter,obj,time,T,method_label] = pg_bt_saf(A,b,s,w0,eta,iter_max,eps);
    end
            
end 
    
%% PDMC, FB and PS
function [w1,iter,objs,time,T,method_label,identify_counter,cvx_identified,T_ls] = pm_saf(method,A,b,s,w0,Q,LQ_coeff,iter_max,eps,acce,sigma,identify,stable_min)
start = tic; 
%Parameters
n = size(A,2);
S = @(x) A*(A'*x);
    
    switch Q
        case 1 %Q=I_n
            dS = 1;
			if (size(A,1) > size(A,2))
				S = @(x) A' * (A*x);
				len = n;
			else
				len = size(A,1);
			end
            LQ = eigs(S,len,1,'largestreal','Tolerance',1e-3); %norm(S)
        otherwise
			S = A*A';
            dS = decomposition(S,'chol');
            LQ = 1;
    end

lambda = LQ_coeff/LQ;
    
%Initialization
iter = 1;
unchanged = 0; 
identify_counter = 0; 
cvx_identified = 0;
T_ls = 0; %time spent on linear system iterations

%First iteration
[PS2w0,I0] = PS2(w0,s);
    
    if strcmp(method,'PS')
        w0 = PS2w0;
    end
    
[nabla0,affine0] = nabla_fQ(w0,A,b,Q,dS);

err = 0.5*norm(affine0)^2 + 0.5*norm(w0-PS2w0)^2;
objs = err;
T = toc(start);

[w1,PS2w1,I1] = newpoint(method,s,w0,lambda,nabla0,PS2w0);
w0 = w1;
  

%Main Loop    
    while (err>eps) && (iter<iter_max)
        iter = iter+1;
        
        [nabla1,affine1] = nabla_fQ(w1,A,b,Q,dS);
        
        err= 0.5*norm(affine1)^2 + 0.5*norm(w1-PS2w1)^2;
        objs = [objs; err];
        
        T = [T; toc(start)];

       
        if identify || acce
            stable = isequal(I0,I1);
            unchanged = stable*(unchanged +1);
            I0 = I1;
        end
        
        if acce && stable
            p = w1 - w0;
                
                if Q == 1
                    old_obj = err;
                else
                    old_obj = 0.5*(norm(nabla1)^2+norm(w1-PS2w1)^2);
                end
                
            [z,PS2z,nabla] = linesearch(method,s,Q,sigma,w1,p,old_obj,affine0,nabla0,affine1,nabla1);
  
            [w,PS2w1,I1] = newpoint(method,s,z,lambda,nabla,PS2z);
            
            %Update new values
            w0 = w1;
            w1 = w;
        
        else%non-accelerated version
			w0 = w1;
            [w1,PS2w1,I1] = newpoint(method,s,w1,lambda,nabla1,PS2w1);
        end    
		nabla0 = nabla1;
		affine0 = affine1;
        

        %If convex set is identified
        if identify && (unchanged>=stable_min) %|| (err<1e-2))
            iter = iter+1;
            
            LS_start = tic; 
            
            identify_counter = identify_counter+1; 
            if (n > 3000 || issparse(A)) %3000 is something that might be altered depending on the computational power
				w = lin_solve(A,b,n,I1,1,sqrt(eps));
            else
				w = lin_solve(A,b,n,I1);
            end
            
            T_ls = T_ls + toc(LS_start); 
            
            err = 0.5*norm(w-PS2(w,s))^2 + 0.5*norm(A*w-b)^2;

            if err<eps 
                T = [T; toc(start)];
                objs = [objs; err];    
                w1 = w;
                cvx_identified = 1;
                break; 
            else 
                T = [T; toc(start)];
                objs = [objs; objs(end)];
                unchanged = -1; %restart
                if mod(identify_counter,10)==0 || identify_counter == 1
                    fprintf('     Incorrectly identified convex set, identify_counter = %d \n',identify_counter);
                end
            end
        end    
    end%end while
% unchanged
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
function [w1,PS2w1,I1] = newpoint(method,s,w0,lambda,nabla0,PS2w0)
    switch method
        case 'PDMC'
            w1 = (w0-lambda*nabla0+lambda*PS2w0)/(1+lambda);
            [PS2w1,I1] = PS2(w1,s);
            
        case 'FB'
            f_step = w0-lambda*nabla0; 
            [PS2f_step,~] = PS2(f_step,s);
            w1 = (lambda*PS2f_step + f_step)/(1+lambda); 
            [PS2w1,I1] = PS2(w1,s);
            
        case 'PS'
            [w1,I1] = PS2(w0-lambda*nabla0,s); 
            PS2w1 = w1; 
    end
    
end

%Line search
function [z,PS2z,nabla] = linesearch(method,s,Q,sigma,w1,p,old_obj,affine0,nabla0,affine1,nabla1)
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
                        [PS2z,~] = PS2(z,s);
                        diff = (norm_obj + t*inn + t^2*norm_dif) + 0.5*norm(z-PS2z)^2  - old_obj + sigma*t^2*norm_p/2;
                end

        case 'PS'
            descent = nabla1'*p;
                %Exact
                switch Q
                    case 1
                        t = max(-2*descent,0)/(norm(affine1-affine0)^2+sigma*norm_p);
                    otherwise
                        t = max(-2*descent,0)/(norm(nabla1-nabla0)^2+sigma*norm_p);
                end 
             
            z = w1+t*p;
            PS2z = [];
    end
    
nabla = (1+t)*nabla1 - t*nabla0;
        
end

%Gradient of f_Q
function [nabla,affine] = nabla_fQ(w,A,b,Q,dS)
%Compute Aw-b
affine = A*w-b;       

%Compute nabla f_Q = A'*Q*(Aw-b)
    if Q==1
        Qobj = affine;
    else
        Qobj = dS\affine;
    end
    
nabla = A'*Qobj;
end

%Projection onto S2
function [z,I] = PS2(w,s)
[~,Ind]=sort(abs(w),'descend');
w(Ind(s+1:length(w)))=0;
z = w; 
I = sort(Ind(1:s));
end

%Linear solver after convex set identification
function w = lin_solve(A,b,n,I1,iterative,eps)
U = A(:,I1);
if (nargin >= 5 && iterative == 1)
	if (nargin < 6)
		eps = sqrt(1e-6/2);
	end
	if (length(I1) >= length(b))
		damping = 1e-10;
		UTUmap = @(x) U'*(U*x) + damping * x; %Ensuring that we can find a solution, and damping should be small so that the output is not too far away from a real solution
		%	D = (1./(full(sum(U.^2))+ damping));%Diagonal entries of U'U as the preconditioner
	else
		UTUmap = @(x) U'*(U*x);
		%	D = (pinv(full(sum(U.^2))+ damping));%Diagonal entries of U'U as the preconditioner
	end
	D = 1;
	nzero = PCG(UTUmap,U'*b,D,eps);
else
	if (length(I1) < length(b))
		nzero = (U'*U)\(U'*b);
	else
		nzero = U\b;
	end
end
w = zeros(n,1);
w(I1) = nzero;
end

%% Projected gradient method with backtracking (Beck and Teboulle)
function [w1,iter,objs,time,T,method_label] = pg_bt_saf(A,b,s,w0,eta,iter_max,eps)
start = tic; 
%Parameters
n = size(A,2); 
    if isempty(eta)
        eta = 1.1;
    end

%Get estimate of initial stepsize T0 (Remark 3.10 [BT11])
ind = randsample(n,s);
x0  = zeros(n,1);
x0(ind) = rand(s,1);
T0 = (norm(A*x0)/norm(x0))^2;

         
%Initialization
[w0,~] = PS2(w0,s);
iter = 0;
err = 1;
T = [];
objs = [];

%Main Loop
    while (err>eps) && (iter<iter_max) 
        iter = iter+1;
        
        c0 = A*w0-b;
        nabla0 = A'*c0;
        
        %Backtracking
        diff = 1;
        
            while diff>0 
                w1 = PS2(w0-(1/T0)*nabla0,s);
                c1 = A*w1-b;
                diff = norm(c1-c0)-sqrt(T0)*norm(w1-w0);
                T0 = eta*T0;
            end
            
        T0 = T0/eta;
        err = 0.5*norm(c1)^2;
        T = [T; toc(start)];
        objs = [objs; err];
        w0 = w1;
    end

%Results
    if err<eps
        fprintf('Solution obtained (PG_BT) \n');
    else
        fprintf('Rerun pg_bt_saf \n');
    end      
    
time = toc(start);

%For graph legend
method_label = 'PG-BT';
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


