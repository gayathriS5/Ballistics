gamma = 1.4;
M_blast = 1.3;
p_ratios = (2*gamma*(M_blast^2)/(gamma+1))-((gamma-1)/(gamma+1));
M_blast = 1.3;
p_inf = 1;
p_b = p_ratios * p_inf;
M_exp = ((gamma-1)*(M_blast^2) + 2)/((2*gamma*(M_blast^2))-(gamma-1));
M_in = 1.95;
p_e = (2*gamma*(M_in^2)/(gamma+1))-((gamma-1)/(gamma+1))*p_inf;
M_tube = ((gamma-1)*(M_in^2) + 2)/((2*gamma*(M_in^2))-(gamma-1));
M_tube_r = 1.1789;
p_stag = p_e*(((gamma-1)*(M_tube_r)^2/2 + 1)^(gamma/(gamma-1)));
M_1 = ((((p_stag/p_b)^((gamma-1)/gamma)-1))*(2/(gamma-1)))^0.5;
mu_0 = asind(1/M_tube_r);
mu_1 = asind(1/M_1);
nu_0 = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_tube_r^2)-1))^0.5))-atand(((M_tube_r^2)-1)^0.5);
nu_1 = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_1^2)-1))^0.5))-atand(((M_1^2)-1)^0.5);
def_1 = nu_1 - nu_0;
nu_2 = nu_1 + def_1;
x = fzero(@(x) myfun(x,gamma,nu_2),2);
M_2 = x;
p_2 = p_stag/((gamma-1)*(M_2)^2/2 + 1)^(gamma/(gamma-1));
p_3 = 1;
M_2_n = (((p_3/p_2)*(gamma+1) + (gamma-1))/(2*gamma))^0.5;
beta_1 = asind(M_2_n/M_2);
num =2*cotd(beta_1)*((M_2^2)*(sind(beta_1)^2)-1);
den = (M_2^2)*(gamma + cosd(2*beta_1)) + 2;
req = num/den;
theta_1 = atand(req);
M_3_n = ((gamma-1)*(M_2_n^2) + 2)/((2*gamma*(M_2_n^2))-(gamma-1));
M_3 = M_3_n/sind(beta_1 - theta_1);


function f = myfun(x,gamma,nu_2)
f = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((x^2)-1))^0.5))-atand(((x^2)-1)^0.5) - nu_2;
end