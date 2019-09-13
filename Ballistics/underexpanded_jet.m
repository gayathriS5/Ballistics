
M_0 = input("Enter exit Mach number :");
P_stag_0 = input("Enter initial stagnation pressure conditions:");
gamma = 1.4;
p_1 = 1;
if M_0<1
    disp("Not a probem");
elseif M_0>1
    p_0 = P_stag_0/((gamma-1)*(M_0)^2/2 + 1)^(gamma/(gamma-1))
    if p_0>p_1
        % expansion fan occurs
        M_1 = ((((P_stag_0/p_1)^((gamma-1)/gamma)-1))*(2/(gamma-1)))^0.5   
        mu_0 = asind(1/M_0);
        mu_1 = asind(1/M_1);
        nu_0 = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_0^2)-1))^0.5))-atand(((M_0^2)-1)^0.5);
        nu_1 = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_1^2)-1))^0.5))-atand(((M_1^2)-1)^0.5);
        def_1 = nu_1 - nu_0
        nu_2 = nu_1 + def_1;
        x = fzero(@(x) myfun(x,gamma,nu_2),2);
        M_2 = x
        p_2 = P_stag_0/((gamma-1)*(M_2)^2/2 + 1)^(gamma/(gamma-1))
        p_3 = 1;
        M_2_n = (((p_3/p_2)*(gamma+1) + (gamma-1))/(2*gamma))^0.5;
        beta_1 = asind(M_2_n/M_2);
        num =2*cotd(beta_1)*((M_2^2)*(sind(beta_1)^2)-1);
        den = (M_2^2)*(gamma + cosd(2*beta_1)) + 2;
        req = num/den;
        theta_1 = atand(req);
        M_3_n = ((gamma-1)*(M_2_n^2) + 2)/((2*gamma*(M_2_n^2))-(gamma-1));
        M_3 = M_3_n/sind(beta_1 - theta_1);
    end  
end
function f = myfun(x,gamma,nu_2)
f = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((x^2)-1))^0.5))-atand(((x^2)-1)^0.5) - nu_2;
end
