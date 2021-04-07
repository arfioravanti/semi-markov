function [cost,status] = simulation_rayleigh_sf(A,B,C,D,K,P,sigma,x0,theta0)
    
    %Configuration
    tmax = 1e5;
    norm_min = 1e-5;
    norm_max = 1e6;
    opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
     
    %Initial Conditions
    t = 0;
    x = x0;
    cost = 0;
    theta = theta0;
    
    %Main Loop
	while(norm(x) >= norm_min && norm(x) <= norm_max && t <= tmax)
        tspan = [t t+sigma(theta)*sqrt(-2*log(1-rand(1)))];
        t0 = tspan(1);
        [tout,xout] = ode45(@(t,x) [(A{theta}+B{theta}*K{theta}(t-t0))*x(1:end-1); x(1:end-1)'*(C{theta}+D{theta}*K{theta}(t-t0))'*(C{theta}+D{theta}*K{theta}(t-t0))*x(1:end-1)], tspan, [x;cost], opts);
        t = tout(end);
        x = xout(end,1:end-1)';
        cost = xout(end,end);
        theta = find(rand(1) < cumsum(P(theta,:)), 1, 'first'); % New Mode
	end
    
    if(norm(x) < norm_min)
        status = 0; % OK
    elseif(norm(x) > norm_max)
        status = 1; % Error: Max Norm
    else
        status = 2; % Error: Time Limit
    end
    
end