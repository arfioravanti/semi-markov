% Example article Rational
clear, close, clc

A = cell(1,2);
A{1} = [0.7 -4; 0 -7];
A{2} = [-7 4; 0 0.7];

C = cell(1,2);
C{1} = eye(2);
C{2} = eye(2);

E = cell(1,2);
E{1} = ones(2,1);
E{2} = ones(2,1);

P = [0 1; 1 0];

mean_ST_vec = 0.01:0.01:2;
mu = [1 0];

order = 4;

%%
exponential = zeros(size(mean_ST_vec));

for i = 1:length(mean_ST_vec)
    mean_ST = mean_ST_vec(i);
    Lambda = [-1/mean_ST 1/mean_ST; 1/mean_ST -1/mean_ST];
    exponential(i) = markov_h2_norm(Lambda, A, E, C, mu);
end

figure
semilogy(mean_ST_vec,exponential)
grid on
%%
rayleigh = zeros(size(mean_ST_vec));

for i = 1:length(mean_ST_vec)
    mean_ST = mean_ST_vec(i);
    sigma = mean_ST*sqrt(2/pi)*ones(2,1);
    rayleigh(i) = rayleigh_h2_norm(sigma, A, E, C, P, mu, order);
end

figure
semilogy(mean_ST_vec,rayleigh)
grid on
%%
erlang = inf(size(mean_ST_vec));
for j = 1:2
    k = [3 3];
    for i = 1:length(mean_ST_vec)
        if(isinf(erlang(i)))
            mean_ST = mean_ST_vec(i) + randn(1)*0.00001;
            lambda = k./mean_ST;
            erlang(i) = erlang_h2_norm(lambda, k, A, E, C, P, mu, 4);
        end
    end
end
figure
semilogy(mean_ST_vec,erlang)
grid on
%%
figure
semilogy(mean_ST_vec,exponential,mean_ST_vec,rayleigh,mean_ST_vec,erlang);
grid on
axis([0 2 0 100]);
xlabel('Mean Sojourn-Time ($\overline{S_t}$)','interpreter','latex');
ylabel('$H_2$ quadratic guaranteed cost','interpreter','latex');
%%
mean_ST = 0.6;
Lambda_exp = 1/mean_ST;
sigma_ray = mean_ST*sqrt(2/pi);
k_erl = 3;
Lambda_erl = k_erl/mean_ST;
Nsim = 1e6;

cost_exponential = zeros(Nsim,1);
status_exponential = zeros(Nsim,1);

parfor i = 1:Nsim
    [cost_exponential(i),status_exponential(i)] = simulation_erlang(A,C,P,Lambda_exp*ones(2,1),ones(2,1),[1;1],1);
end


% cost_rayleigh = zeros(Nsim,1);
% status_rayleigh = zeros(Nsim,1);
% 
% parfor i = 1:Nsim
%     [cost_rayleigh(i),status_rayleigh(i)] = simulation_rayleigh(A,C,P,sigma_ray*ones(2,1),[1;1],1);
% end

% cost_erlang = zeros(Nsim,1);
% status_erlang = zeros(Nsim,1);
% 
% parfor i = 1:Nsim
%     [cost_erlang(i),status_erlang(i)] = simulation_erlang(A,C,P,Lambda_erl*ones(2,1),k_erl*ones(2,1),[1;1],1);
% end
%%
clear h;
syms h

    numQ1 = [vc(1)+h*vc(2)+h^2*vc(3)+h^3*vc(4)+h^4*vc(5), vc(6)+h*vc(7)+h^2*vc(8)+h^3*vc(9)+h^4*vc(10);
             vc(6)+h*vc(7)+h^2*vc(8)+h^3*vc(9)+h^4*vc(10), vc(11)+h*vc(12)+h^2*vc(13)+h^3*vc(14)+h^4*vc(15)];

    numQ2 = [vc(46)+h*vc(47)+h^2*vc(48)+h^3*vc(49)+h^4*vc(50), vc(51)+h*vc(52)+h^2*vc(53)+h^3*vc(54)+h^4*vc(55);
             vc(51)+h*vc(52)+h^2*vc(53)+h^3*vc(54)+h^4*vc(55), vc(56)+h*vc(57)+h^2*vc(58)+h^3*vc(59)+h^4*vc(60)];

%% Example 2 article Rational
clear, close, clc

A = cell(1,2);
A{1} = [0.7 -4; 0 -7];
A{2} = [-7 4; 0 0.7];

C = cell(1,2);
C{1} = [eye(2);zeros(1,2)];
C{2} = [eye(2);zeros(1,2)];

E = cell(1,2);
E{1} = ones(2,1);
E{2} = ones(2,1);

B = cell(1,2);
B{1} = [0;1];
B{2} = [0;1];

D = cell(1,2);
D{1} = [0;0;0.5];
D{2} = [0;0;0.5];
P = [0 1; 1 0];

mean_ST = 0.6;
sigma_ray = mean_ST*sqrt(2/pi);

mu = [1 0];

order = 4;
[norm,sol,v,Qq,res] = sf_rayleigh3(sigma_ray*ones(2,1),A,B,E,C,D,P,mu,order);

%%
% vc = value(c);
% syms t
% 
% vX1 = [vc(1)+t*vc(2)+t^2*vc(3)+t^3*vc(4)+t^4*vc(5), vc(6)+t*vc(7)+t^2*vc(8)+t^3*vc(9)+t^4*vc(10);
%        vc(6)+t*vc(7)+t^2*vc(8)+t^3*vc(9)+t^4*vc(10), vc(11)+t*vc(12)+t^2*vc(13)+t^3*vc(14)+t^4*vc(15)];
%    
% vX2 = [vc(332)+t*vc(333)+t^2*vc(334)+t^3*vc(335)+t^4*vc(336), vc(337)+t*vc(338)+t^2*vc(339)+t^3*vc(340)+t^4*vc(341);
%        vc(337)+t*vc(338)+t^2*vc(339)+t^3*vc(340)+t^4*vc(341), vc(342)+t*vc(343)+t^2*vc(344)+t^3*vc(345)+t^4*vc(346)];
%    
% vY1 = [vc(16)+t*vc(17)+t^2*vc(18)+t^3*vc(19), vc(20)+t*vc(21)+t^2*vc(22)+t^3*vc(23)];
% vY2 = [vc(347)+t*vc(348)+t^2*vc(349)+t^3*vc(350), vc(351)+t*vc(352)+t^2*vc(353)+t^3*vc(354)];

K = cell(1,2);
K{1} = @(t) [ (0.0625*(1.621885e+34*t^7 - 5.724701e+34*t^6 + 1.050687e+35*t^5 - 1.220608e+35*t^4 + 8.918467e+34*t^3 - 3.137756e+34*t^2 + 1.976129e+32*t + 2.796155e+33))/(2.603853e+31*t^8 - 2.998711e+32*t^7 + 1.096179e+33*t^6 - 1.441013e+33*t^5 + 1.082775e+33*t^4 - 2.67714e+32*t^3 + 7.119712e+31*t^2 - 5.894894e+31*t + 9.14449e+31), -(0.125*(2.177302e+34*t^7 - 5.007344e+34*t^6 + 5.510761e+34*t^5 - 3.526651e+34*t^4 + 1.827408e+34*t^3 - 9.445499e+33*t^2 + 2.213259e+33*t + 8.797739e+32))/(2.603853e+31*t^8 - 2.998711e+32*t^7 + 1.096179e+33*t^6 - 1.441013e+33*t^5 + 1.082775e+33*t^4 - 2.67714e+32*t^3 + 7.119712e+31*t^2 - 5.894894e+31*t + 9.14449e+31)];
K{2} = @(t) [ -(1.0*(- 1.363362e+34*t^7 + 3.314473e+34*t^6 - 2.601512e+34*t^5 + 2.468506e+33*t^4 + 5.095799e+33*t^3 + 1.737892e+33*t^2 - 2.054468e+33*t + 2.392531e+32))/(2.786361e+32*t^8 - 3.02114e+33*t^7 + 1.021047e+34*t^6 - 1.143719e+34*t^5 + 5.962093e+33*t^4 - 3.172942e+32*t^3 + 1.741551e+33*t^2 - 2.571381e+33*t + 1.30206e+33), (- 1.141806e+34*t^7 + 5.463893e+34*t^6 - 1.116549e+35*t^5 + 1.115335e+35*t^4 - 5.680413e+34*t^3 + 8.173146e+33*t^2 + 6.176738e+33*t - 3.270224e+33)/(2.786361e+32*t^8 - 3.02114e+33*t^7 + 1.021047e+34*t^6 - 1.143719e+34*t^5 + 5.962093e+33*t^4 - 3.172942e+32*t^3 + 1.741551e+33*t^2 - 2.571381e+33*t + 1.30206e+33)];


t = 0:0.01:2;
vK1 = zeros(length(t),2);
vK2 = zeros(length(t),2);

for i = 1:length(t)
    h = t(i);
    vK1(i,:) = K{1}(h);
    vK2(i,:) = K{2}(h);
end

subplot(1,2,1);
plot(t,vK1,'LineWidth',2);
xlabel('h');
ylabel('K_1(1,1) [blue] and K_1(1,2) [red]');
axis([0 2 -15 5])

subplot(1,2,2);
plot(t,vK2,'LineWidth',2);
xlabel('h');
ylabel('K_2(1,1) [blue] and K_2(1,2) [red]');
axis([0 2 -15 5])

%%
Nsim = 3e4;

cost_rayleigh = zeros(Nsim,1);
status_rayleigh = zeros(Nsim,1);

for i = 1:Nsim
    [cost_rayleigh(i),status_rayleigh(i)] = simulation_rayleigh_sf(A,B,C,D,K,P,sigma_ray*ones(2,1),[1;1],1);
end

%% Exemplo 2
clear
close
clc

sigma = sqrt(2/pi)./[0.3 1.5 2];
    
Ahat = [ -11.4540  2.7185 -19.4399 0;
           0.5068 -2.9875  23.3434 0;
           0.0922 -0.9957  -0.4680 0.3256;
           1       0.0926   0      0];
       
Bhat = [ 78.4002 -3.4690 0 0;
         -2.7282 13.9685 0 0]';
     
Jhat = eye(4);

Chat = [eye(2) zeros(2);zeros(2,4)];
Dhat = [zeros(2);eye(2)];

A = cell(1,3);
B = cell(1,3);
E = cell(1,3);
C = cell(1,3);
D = cell(1,3);

for i = 1:3
    A{i} = Ahat;
    B{i} = Bhat;
    E{i} = Jhat;
    C{i} = Chat;
    D{i} = Dhat;
end

B{2} = Bhat*[0 0; 0 1];
B{3} = B{2};
B{3} = Bhat*[0 0; 0 1/2];

Pi = [   0   2/3 1/3;
       11/15  0  4/15;
        1/2  1/2  0   ];
    
mu = [1 0 0];

custo = zeros(1,7);
tempo_parser = zeros(1,7);
tempo_solver = zeros(1,7);
for idx = 1:7
    order = 2*idx;
    [norm,sol,v,Qq,res] = sf_rayleigh3(sigma,A,B,E,C,D,Pi,mu,order);
    tempo_parser(idx) = sol.yalmiptime;
    tempo_solver(idx) = sol.solvertime;
    custo(idx) = norm;
end

%%
Lambda = [-0.3 0.2 0.1; 1.1 -1.5 0.4; 1 1 -2];
[norma,vK] = markov_h2_sf(Lambda, A, B, E, C, D, mu);

Nsim = 3e4;

cost_rayleigh1 = zeros(Nsim,1);
status_rayleigh1 = zeros(Nsim,1);
cost_rayleigh2 = zeros(Nsim,1);
status_rayleigh2 = zeros(Nsim,1);
cost_rayleigh3 = zeros(Nsim,1);
status_rayleigh3 = zeros(Nsim,1);
cost_rayleigh4 = zeros(Nsim,1);
status_rayleigh4 = zeros(Nsim,1);

K = cell(1,3);
K{1} = @(t) vK{1};
K{2} = @(t) vK{2};
K{3} = @(t) vK{3};

parfor i = 1:Nsim
    [cost_rayleigh1(i),status_rayleigh1(i)] = simulation_rayleigh_sf(A,B,C,D,K,Pi,sigma,[1;0;0;0],1);
    [cost_rayleigh2(i),status_rayleigh2(i)] = simulation_rayleigh_sf(A,B,C,D,K,Pi,sigma,[0;1;0;0],1);
    [cost_rayleigh3(i),status_rayleigh3(i)] = simulation_rayleigh_sf(A,B,C,D,K,Pi,sigma,[0;0;1;0],1);
    [cost_rayleigh4(i),status_rayleigh4(i)] = simulation_rayleigh_sf(A,B,C,D,K,Pi,sigma,[0;0;0;1],1);
end
