clear all;
clc;

Q_val = [0;50000;100000];
    %Costs associated with acutely symptomatic population
thetaz_val = [1;500;5000;50000000000];
    %Costs associated with vaccines
C_dth = [0; 3000; 6000; 9000; 12000; 15000; 20000];
    %Costs associated with deceased population

Rho = 5.08;
    %Basic reproduction number for D-variant
H_in = 0.06925;
    %Percentage of hospitalized - range between 5% and 12%
a_d = 0.0066/H_in;
    %so the infection mortality rate is 0.66%

gamma_i = 1/14;
    %Recovery rate from infected undetected
gamma_d = 1/14;
    %Recovery rate from infected detected
gamma_a = 1/12.39;
    %Recovery rate from hospitalized
gamma_w = 1/14;
    %Recovery rate from Vaccinated Infected undetected
gamma_y = 1/12.39;
    %Recovery rate from Vaccinated Acutely Symptomatic

ksi_i = H_in/(1-H_in)*gamma_i;
    %Transition rate from infected undetected to acutely symptomatic
ksi_d = H_in/(1-H_in)*gamma_d;
    %Transition rate from infected detected to acutely symptomatic       
ksi_w = 0.000265;
    %Transition rate from Vaccinated infected Undetected
    %to Vaccinated acutely symptomatic

mu = a_d/(1-a_d)*gamma_a;
    %Transition rate from acutely symptomatic to Deceased
muy = 0.0085;
    %Transition rate from Vaccinated acutely symptomatic to Deceased

beta = Rho*(gamma_i + ksi_i);
    %Definition of R0 in SIDARE, proven in our paper          
beta1 = 0.19495;
    %Rate which a vaccinated person infect an unvaccinated person
beta2 = 0.05849;
    %Rate which a vaccinated person infect another vaccinated person
beta3 = 0.11697;
    %Rate which an unvaccinated person infect a vaccinated person

v_val=0.05;                         %Testing rate value adopted

dt = 1;                     %time increments  
N = length(C_dth);          %number of iterations

               
%***Q = 0 or 50000 or 100000***
Q = diag([0;0;0;Q_val(3,1);0;0]);  %Cost associated with states

%***theta_z = 1 or 500 or 5000 or 50000000000***
theta_z = thetaz_val(4,1);

f=1;
parfor i=1 + (f-1)*N:N + (f-1)*N
[x{i}, u(i,:),zeta(i,:), C(:,i), C1(:,i), C2(:,i), C3(:,i), C4(:,i)] = Sim_simple(dt, beta, beta1, beta2, beta3, gamma_i, gamma_d, gamma_a, gamma_w, gamma_y, ksi_i, ksi_d, ksi_w, mu, muy, C_dth(i-(f-1)*N,1), Q, v_val, theta_z);
end
