function  norm = rayleigh_h2_norm(sigma, A, E, C, P, mu, order)

    n = size(A{1},1);
    m = size(E{1},2);
    p = size(C{1},1);
    N = numel(A);
    
    eps = 1e-6;
    h = sdpvar; % independent variable
    Q = cell(1,N);
    M1 = cell(1,N);
    c = [];

    for i = 1:N
        Q{i} = zeros(n, n, 'like', sdpvar);
        for j = 1:n
            for tk = 1:j
                [Q{i}(j,tk),ct] = polynomial(h,order);
                c = [c;ct];
                if(j ~= tk)
                    Q{i}(tk,j) = Q{i}(j,tk);
                end
            end
        end
        
        M1{i} = zeros(n, n, 'like', sdpvar);
        for j = 1:n
            for tk = 1:j
                [M1{i}(j,tk),ct] = polynomial(h,2*order+1);
                c = [c;ct];
                if(j ~= tk)
                    M1{i}(tk,j) = M1{i}(j,tk);
                end
            end
        end
        
    end

    F = [];
    d = (1+h^order);
    ddot  = order*h^(order-1);

    rho = zeros(n, 1, 'like', sdpvar);
    eta = zeros(n, 1, 'like', sdpvar);
    
    obj = 0;
    for i = 1:N 
        rho(i) = 1;
        eta(i) = h/(sigma(i)^2);
        
        T = rho(i)*d*(A{i}'*Q{i}+Q{i}*A{i}+jacobian(Q{i},h)+d*C{i}'*C{i}) - rho(i)*ddot*Q{i} - d*eta(i)*Q{i} + eps*rho(i)*d^2*eye(n);
        for j = 1:N
            if (i == j)
                continue;
            end
            T = T + d^2*eta(i)*P(i,j)*replace(Q{j},h,0);
        end
        
        F = [F, sos(Q{i}), sos(M1{i}), sos(-T - h*M1{i})];
        obj = obj + mu(i)*trace(E{i}'*replace(Q{i},h,0)*E{i});
    end
    
    options = sdpsettings('solver','mosek');
    options.sos.model = 2;
    options.sos.scale = 0;
    options.sos.congruence = 1;
    options.sos.numblk = 1e-6;
    
    [sol,v,Qq,res] = solvesos(F,obj,options,c);
    if(sol.problem == 0)
        norm = value(obj);
    else
        norm = inf;
    end
end
