clc;
clear all;
global M_blast P_atm P_e  theta_f exp_fan_f
global p_e p_b gamma theta D x_co_f
M_blast = input('Enter the Mach number of the blast wave initially coming out of the tube:');
P_atm = input('Enter the ambient pressure into which the blast wave is entering into:');
P_e = input('Enter the pressure at the exit of tube from experiment(in bar)in the front case:');
gamma = input('Enter gamma:');
D = input('Enter Diameter of the tube :');
n_f = input('Enter number of expansion fans excluding first exp_fan: ');
P_b = (2*gamma*(M_blast^2)/(gamma+1))-((gamma-1)/(gamma+1))*P_atm;
%before expansion fan starts just at the end of tube
%lets consider that the flow is just choked that is M=1
i = 1;
M_f(i) = 1;
pr_f(i) = P_e;
p_stag_f = P_e*(((gamma-1)*(M_f(1))^2/2 + 1)^(gamma/(gamma-1)));
M_b_1 = ((((p_stag_f/P_b)^((gamma-1)/gamma)-1))*(2/(gamma-1)))^0.5;
mu0 = asind(1/M_f(1));
mu1 = asind(1/M_b_1);
nu_0 = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_f(1)^2)-1))^0.5))-atand(((M_f(1)^2)-1)^0.5);
nu_1 = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_b_1^2)-1))^0.5))-atand(((M_b_1^2)-1)^0.5);
theta_f = nu_1 - nu_0;
exp_fan_f = mu0+(theta_f)-mu1;
def_par_f = (theta_f)/n_f;

tot_points_f = n_f*(n_f+1)/2 + n_f;
table_f = ones(tot_points_f,7);
header = {'POINT','THETA','NU','THETA+NU','THETA-NU','MACH ANGLE','MACH NO'};
table_f = [header; num2cell(table_f)];
%as we took initial M_f = 1 so theta == nu as nu = 0 for M=1
i = n_f;
a = 0;
while i>=1
    for k = 1:i+1
        a = a+1;
        if a == 1
            table_f(a,1) = a;
            table_f(a,2) = 0;
            table_f(a,3) = 0;
            table_f(a,4) = table_f(a,2)+table_f(a,3);
            table_f(a,5) = table_f(a,2)-table_f(a,3);
            table_f(a,7) = 1;
            table_f(a,6) = asind(1/table_f(a,7)); 
        end
        table_f(a,1) = a;
        table_f(a,4) = table_f(a,2)+table_f(a,3);
        table_f(a,5) = table_f(a,2)-table_f(a,3);
        table_f(a,7) = 1;
        table_f(a,6) = asind(1/table_f(a,7));
    end
    i = i-1;
end
function f = prandtl(x,nu_2)
global gamma 
f = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((x^2)-1))^0.5))-atand(((x^2)-1)^0.5) - nu_2;
end