function [norm,sol,v,Qq,res] = sf_rayleigh3(sigma,A,B,E,C,D,Pi,mu,order)

    n = size(A{1},1);
    m = size(E{1},2);
    p = size(C{1},1);
    q = size(B{1},2);
    N = numel(A);

    eps = 1e-6  ;
    h = sdpvar; % independent variable
    X = cell(1,N);
    M1 = cell(1,N);
    W = cell(1,N);
    Y = cell(1,N);
    Z = cell(N,N);
    c = [];
    
    d = cell(1,N);
    ddot = cell(1,N);
    for i = 1:N
        d{i} = 1+h^order;
        ddot{i} = order*h^(order-1);
    end
    dim_d = order;

    for i = 1:N
        W{i} = sdpvar(m,m,'symmetric');
        X{i} = zeros(n, n, 'like', sdpvar);
        for j = 1:n
            for tk = 1:j
                [X{i}(j,tk),ct] = polynomial(h,order);
                c = [c;ct];
                if(j ~= tk)
                    X{i}(tk,j) = X{i}(j,tk);
                end
            end
        end
        Y{i} = zeros(q, n, 'like', sdpvar);
        for j = 1:q
            for tk = 1:n
                [Y{i}(j,tk),ct] = polynomial(h,order-1);
                c = [c;ct];
            end
        end
        M1{i} = zeros(2*n+p, 2*n+p, 'like', sdpvar);
        for j = 1:2*n+p
            for tk = 1:j
                [M1{i}(j,tk),ct] = polynomial(h,order+dim_d+2);
                c = [c;ct];
                if(j ~= tk)
                    M1{i}(tk,j) = M1{i}(j,tk);
                end
            end
        end  
    end
    
    for i = 1:N
        for j = 1:N
            if(j == i)
                continue
            end
            Z{i,j} = zeros(n, n, 'like', sdpvar);
            for ti = 1:n
                for tj = 1:ti
                    [Z{i,j}(ti,tj),ct] = polynomial(h,2*order);
                    c = [c;ct];
                    if(tj ~= ti)
                        Z{i,j}(tj,ti) = Z{i,j}(ti,tj);
                    end
                end
            end
        end
    end

    F = [];

    obj = 0;
    
    rho = zeros(n, 1, 'like', sdpvar);
    eta = zeros(n, 1, 'like', sdpvar);
    
    for i = 1:N 

        rho(i) = 1;
        eta(i) = h/(sigma(i)^2);
        
        T = zeros(2*n+p,2*n+p,'like', sdpvar);
        T(1:n,1:n) = rho(i)*d{i}*(A{i}*X{i}+X{i}*A{i}'+B{i}*Y{i}+Y{i}'*B{i}'-jacobian(X{i},h))+(rho(i)*ddot{i}-eta(i)*d{i})*X{i};
        for j = 1:N
            if (i == j)
                continue;
            end
            T(1:n,1:n) = T(1:n,1:n) + eta(i)*Pi(i,j)*Z{i,j};
        end
        T(n+1:n+p,1:n) = rho(i)*(C{i}*X{i} + D{i}*Y{i});
        T(1:n,n+1:n+p) = transpose(T(n+1:n+p,1:n));
        T(n+1:n+p,n+1:n+p) = -rho(i)*eye(p);
        T(n+p+1:2*n+p,1:n) = rho(i)*X{i};
        T(1:n,n+p+1:2*n+p) = transpose(T(n+p+1:2*n+p,1:n));
        T(n+p+1:2*n+p,n+p+1:2*n+p) = -eps^(-1)*rho(i)*eye(n);
        
        T2 = [W{i}, E{i}'; E{i}, replace(X{i},h,0)];
        F = [F, sos(X{i}-eps*d{i}*eye(n)), sos(eps^(-1)*d{i}*eye(n)-X{i}), sos(M1{i}-eps*eye(2*n+p)), sos(-T - h*M1{i} - eps*eye(2*n+p)), T2 >= eps*eye(m+n)];
        obj = obj + mu(i)*trace(W{i});
        for j = 1:N
            if (j == i)
                continue
            end
            F = [F, sos(Z{i,j}-eps*eye(n)), sos([Z{i,j}, X{i}; X{i}, replace(X{j},h,0)] - eps*eye(2*n))];
        end
    end
    options = sdpsettings('solver','mosek');
    options.sos.model = 1;
    options.sos.scale = 1;
    options.sos.congruence = 1;
    options.sos.numblk = 1e-6;
    
    [sol,v,Qq,res] = solvesos(F,obj,options,c);
    norm = value(obj);
end
