function [y,dy] = epidem(dt, x, beta, beta1, beta2, beta3, u, v, zeta, gamma_i, gamma_d, gamma_a,  gamma_w, gamma_y, ksi_i, ksi_d, ksi_w, mu, muy)
 
    
    %Controlled SIDAVERHY model
    dy(1,1) = -beta*(1 - u)*x(1,1)*x(2,1) - zeta*x(1,1) - x(8,1)*x(1,1)*beta1;                                              %Susceptible State
    dy(2,1) = beta*(1 - u)*x(1,1)*x(2,1) - gamma_i*x(2,1) - ksi_i*x(2,1) - v*x(2,1) - zeta*x(2,1) + x(8,1)*x(1,1)*beta1;    %Infected undetected State
    dy(3,1) = v*x(2,1) - gamma_d*x(3,1) - ksi_d*x(3,1);                                                                     %Detected infected State
    dy(4,1) = ksi_i*x(2,1) +  ksi_d*x(3,1) - gamma_a*x(4,1) - mu*x(4,1);                                                    %Acutely symptomatic State
    dy(5,1) = gamma_i*x(2,1) +  gamma_d*x(3,1) + gamma_a*x(4,1) + gamma_w*x(8,1) + gamma_y*x(9,1);                          %Recovered State
    dy(6,1) = mu*x(4,1) + muy*x(9,1);                                                                                       %Extinct (Deceased) State
    dy(7,1) = zeta*x(1,1) - x(8,1)*x(7,1)*beta2 - x(2,1)*x(7,1)*beta3*(1-u);                                                %Vaccination State
    dy(8,1) = x(2,1)*x(7,1)*beta3*(1-u) + x(8,1)*x(7,1)*beta2 + zeta*x(2,1) - gamma_w*x(8,1) - ksi_w*x(8,1);                %Vacc. Infected
    dy(9,1) = ksi_w*x(8,1) - gamma_y*x(9,1) - muy*x(9,1) ;
   y = max(x + dt*dy,0);                                                                                                    %State update
   
end

