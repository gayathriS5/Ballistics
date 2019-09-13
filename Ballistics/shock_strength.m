function f = shock_strength(x)
global p_e p_b gamma theta 
f(1) = x(2) - (((((p_e/p_b)^((gamma-1)/gamma)-1))*(2/(gamma-1)))^0.5)*x(1) ;

f(2) = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((x(2)^2)-1))^0.5))-atand(((x(2)^2)-1)^0.5)...
    -((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((x(1)^2)-1))^0.5))-atand(((x(1)^2)-1)^0.5) - theta;
end