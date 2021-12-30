%function of costate variables update based on pontryagin's minimum principle

function [y,dy] = pontr(dt, l, x, u,v,zeta, beta, beta1, beta2, beta3,gamma_i, gamma_d, gamma_a,gamma_w,gamma_y, ksi_i, ksi_d, ksi_w, mu, theta_a, muy)

    
    %Equations based on Pontryagin's minimum principle
     dy = -[l(1,1)*(-beta*x(2,1)+beta*x(2,1)*u-zeta-x(8,1)*beta1) + l(2,1)*(beta*x(2,1)-beta*x(2,1)*u+x(8,1)*beta1) + l(7,1)*zeta;
         l(1,1)*(-beta*x(1,1)+beta*x(1,1)*u) + l(2,1)*(beta*x(1,1)-gamma_i-ksi_i-v-beta*x(1,1)*u-zeta) + l(3,1)*v + l(4,1)*ksi_i + l(5,1)*gamma_i + l(7,1)*(-x(7,1)*beta3+x(7,1)*beta3*u) + l(8,1)*(x(7,1)*beta3-x(7,1)*beta3*u+zeta);
         l(3,1)*(-gamma_d-ksi_d) + l(4,1)*ksi_d + l(5,1)*gamma_d;
         l(4,1)*(-gamma_a-mu) + l(5,1)*gamma_a + l(6,1)*mu + theta_a*(x(4,1)+x(9,1));
         0;
         0;
         l(7,1)*(-x(8,1)*beta2-x(2,1)*beta3+x(2,1)*beta3*u) + l(8,1)*(x(2,1)*beta3+x(8,1)*beta2-x(2,1)*beta3*u) ;
         l(1,1)*(-x(1,1)*beta1) + l(2,1)*(x(1,1)*beta1) + l(5,1)*gamma_w - l(7,1)*x(7,1)*beta2 + l(8,1)*(x(7,1)*beta2-gamma_w-ksi_w) + l(9,1)*ksi_w ;
         l(5,1)*gamma_y + l(9,1)*(-gamma_y-muy) + l(6,1)*muy + theta_a*(x(4,1)+x(9,1));
         ];
    y = l - dy*dt;  %backwards in time, since boundary condition for costate variables is at t = T
    
end

