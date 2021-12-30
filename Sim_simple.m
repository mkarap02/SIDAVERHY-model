function [x,u,zeta,C,C1,C2,C3,C4] = Sim_simple(dt, beta, beta1, beta2, beta3, gamma_i, gamma_d, gamma_a, gamma_w, gamma_y, ksi_i, ksi_d, ksi_w, mu, muy, C_dth, Q, v_val, theta_z)

T_days = 365;   %Number of days

R = 1;              %Cost associated with government strategy (used as basis)
theta_a = Q(4,4);   %Cost associated with the threatened population 
z_max = 0.05;

%Initial conditions
r = 0.00001;
x(1,1) = 1 - r; %S : Susceptible
x(2,1) = r;     %I : Infected Undetected
x(3,1) = 0;     %D : Infected Detected
x(4,1) = 0;     %A : Acutely Symptomatic
x(5,1) = 0;     %R : Recovered
x(6,1) = 0;     %E : Extinct
x(7,1) = 0;     %V : Vaccinated
x(8,1) = 0;     %Y : Vaccinated Infected Undetected
x(9,1) = 0;     %H : Vaccinated Acutely Symptomatic 

T = T_days/dt; 
l(1:length(x(:,1)),T) = 0;  %Lambda boundary conditions
l(6,T) = C_dth;             %Cost attributed to number of deaths
u_max = 0.8;                %Maximum value for u

u(1:T,1) = 0.4;             %Initialisation of u
v(1:T,1) = v_val;           %Constant value of testing rate, ν
zeta(1:T,1) = 0.1;          %Initialisation of z


%Initialization of states and costs
for k=2:T
x(:,k) = epidem(dt, x(:,k-1), beta(1,1), beta1(1,1), beta2(1,1), beta3(1,1), u(k-1,1),v(k-1,1), zeta(k-1,1) , gamma_i, gamma_d, gamma_a, gamma_w, gamma_y, ksi_i, ksi_d, ksi_w, mu, muy);
end

for k=T-1:-1:1
[l(:,k), dl(:,k)] = pontr(dt, l(:,k+1), x(:,k+1), u(k+1,1), v(k+1,1), zeta(k+1,1) , beta(1,1), beta1(1,1), beta2(1,1), beta3(1,1), gamma_i, gamma_d, gamma_a, gamma_w, gamma_y, ksi_i, ksi_d, ksi_w, mu, theta_a, muy);
end


%Cost function - aggregate and components----------------------------------
C1(1,1) = 0.5*dt*(R(1,1)*(u.'*u));                              %cost associated with government strategy u
C2(1,1) = 0.5*dt*theta_a*(x(4,:)+x(9,:))*((x(4,:)+x(9,:)).');   %cost associated with the acutely symptomatic population
C3(1,1) = x(length(x(:,1))-3,T)*C_dth;                          %cost associated with number of deaths
C4(1,1) = 0.5*dt*(zeta.'*zeta)*theta_z;                         %cost associated with vaccination
C(1,1) = C1(1,1)+C2(1,1)+C3(1,1)+C4(1,1);                       %Total Cost


N_iter = 100000;               %number of iterations for the convergence of the algorithm

for j=1:N_iter

    %Calculation of the new value for u and zeta
    u0 = u;
    zeta0 = zeta;
    for k=1:T
        u1(k,1) = min(max(inv(R(1,1))*(beta(1,1)*x(1,k)*x(2,k)*(l(2,k)-l(1,k)) + x(2,k)*x(7,k)*beta3*(l(8,k)-l(7,k))) ,0),u_max);
        zeta1(k,1) = min(max((l(1,k)*x(1,k) + l(2,k)*x(2,k) - l(7,k)*x(1,k) - l(8,k)*x(2,k)) / theta_z , 0),z_max);
        
    end
   
    a = 0.99;                      %coefficient used to update the current u 
    u = a*u0 + (1-a)*u1;            %new strategy u
    zeta = a*zeta0 + (1-a)*zeta1;
    
    %Update the SIDAVERHY model trajectory based on current u, zeta
    for k=2:T
        %Controlled SIDAVERHY epidemic model
        x(:,k) = epidem(dt, x(:,k-1), beta(1,1), beta1(1,1), beta2(1,1), beta3(1,1), u(k-1,1), v(k-1,1), zeta(k-1,1), gamma_i, gamma_d, gamma_a,  gamma_w, gamma_y, ksi_i, ksi_d, ksi_w, mu,muy);
    end
    
    %Update the costate variables
    for k=T-1:-1:1
        %Pontryagin equations
        [l(:,k), dl(:,k)] = pontr(dt, l(:,k+1), x(:,k+1), u(k+1,1), v(k+1,1),zeta(k+1,1), beta(1,1), beta1(1,1), beta2(1,1), beta3(1,1), gamma_i, gamma_d, gamma_a,  gamma_w, gamma_y, ksi_i, ksi_d, ksi_w, mu, theta_a, muy);
    end


    %Cost function - aggregate and components associated with iteration j
    C1(j,1) = 0.5*dt*(R(1,1)*(u.'*u));                                  %cost associated with government strategy u
    C2(j,1) = 0.5*dt*theta_a*((x(4,:)+x(9,:))*((x(4,:)+x(9,:)).'));     %cost associated with the acutely symptomatic population
    C3(j,1) = x(length(x(:,1))-3,T)*C_dth;                              %cost associated with number of deaths
    C4(j,1) = 0.5*dt*(zeta.'*zeta)*theta_z;                             %cost associated with vaccination
    C (j,1) = C1(j,1) + C2(j,1) + C3(j,1) + C4(j,1);                    %Total Cost
    
end