gamma = 1.4;
M_blast = 1.3;
p_ratios = (2*gamma*(M_blast^2)/(gamma+1))-((gamma-1)/(gamma+1));
M_blast = 1.3;
p_atm = 1;
p_b = p_ratios * p_atm; %due to blast wave
Temp = 308;
M_tube_r = 1;
M_2_tube = 2.1;
p_ratios_2 = (2*gamma*(M_2_tube^2)/(gamma+1))-((gamma-1)/(gamma+1));
p_e = p_ratios_2 * p_atm;  % pressure at the end of tube
p_stag = p_e*(((gamma-1)*(M_tube_r)^2/2 + 1)^(gamma/(gamma-1)));
M_1 = ((((p_stag/p_b)^((gamma-1)/gamma)-1))*(2/(gamma-1)))^0.5; %after first expansion fan
mu_0 = asind(1/M_tube_r);
mu_1 = asind(1/M_1);
nu_0 = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_tube_r^2)-1))^0.5))-atand(((M_tube_r^2)-1)^0.5);
nu_1 = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_1^2)-1))^0.5))-atand(((M_1^2)-1)^0.5);
def_1 = nu_1 - nu_0;
nu_2 = nu_1 + def_1;
x = fzero(@(x) myfun(x,gamma,nu_2),2);
M_2 = x; %after second expansion fan
p_2 = p_stag/((gamma-1)*(M_2)^2/2 + 1)^(gamma/(gamma-1));
M_2_n = (((p_b/p_2)*(gamma+1) + (gamma-1))/(2*gamma))^0.5;
beta_1 = asind(M_2_n/M_2);
beta_1 = asind(M_2_n/M_2);
num =2*cotd(beta_1)*((M_2^2)*(sind(beta_1)^2)-1);
den = (M_2^2)*(gamma + cosd(2*beta_1)) + 2;
req = num/den;
theta_1 = atand(req);
M_3_n = ((gamma-1)*(M_2_n^2) + 2)/((2*gamma*(M_2_n^2))-(gamma-1));
M_3 = M_3_n/sind(beta_1 - theta_1); %after oblique incident shock
% a mach reflection occurs
M_4 = ((gamma-1)*(M_2^2) + 2)/((2*gamma*(M_2^2))-(gamma-1));%Normal shock
p_ratios_3 = (2*gamma*(M_2^2)/(gamma+1))-((gamma-1)/(gamma+1));
p_4 = p_ratios_3*p_2;
p_pro_b = 12;

%symmetric body of length :50 width :24 and angle at the front is 30deg and
%thickness 't'
%initially when barel is between p_pro_b and exit
%t = input("Enter thickness of projectile:")
%area = 24*t;%in mm^2
%load_1 = ((p_pro_b)*(area) - (p_e*cosd(60)*2*area))/10;
%load_2 = ((p_pro_b)*(area) - (p_b*cosd(60)*2*area))/10;
%load_3 = ((p_pro_b)*(area) - (p_2*cosd(60)*2*area))/10;
%load_4 = ((p_pro_b)*(area) - (p_4*cosd(60)*2*area))/10;
%p = [p_e,p_b,p_2,p_4];
%for m = 1:length(p)
%load = ((p_pro_b)*(area) - (p(m)*cosd(60)*2*area))/10;
%x = [1,2,3,4];
%plot(x(m),load,'r'); hold on;
%end
%hold off
p_stag_2 = p_pro_b*(((gamma-1)*(M_tube_r)^2/2 + 1)^(gamma/(gamma-1)));
%consider M=1 at the exit again of the tube
M_1_b = ((((p_stag_2/p_b)^((gamma-1)/gamma)-1))*(2/(gamma-1)))^0.5;
mu_0_b = asind(1/M_tube_r);
mu_1_b =  asind(1/M_1);
nu_0_b = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_tube_r^2)-1))^0.5))-atand(((M_tube_r^2)-1)^0.5);
nu_1_b = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_1_b^2)-1))^0.5))-atand(((M_1_b^2)-1)^0.5);
def_1_b = nu_1_b - nu_0_b;
nu_2_b = nu_1_b + def_1_b;
fan_angle = mu_0_b + def_1_b - mu_1_b;
def_par = def_1_b/3;
y = fzero(@(y) myfun(y,gamma,def_par),2);
M_expf_1 = y;
p_expf_1 = p_stag_2/(((gamma-1)*(M_expf_1)^2/2 + 1)^(gamma/(gamma-1)));
z = fzero(@(z) myfun(z,gamma,2*def_par),2);
M_expf_2 = z;
p_expf_2 = p_stag_2/(((gamma-1)*(M_expf_2)^2/2 + 1)^(gamma/(gamma-1)));
w = fzero(@(w) myfun(w,gamma,3*def_par),2);
M_expf_3 = w;
p_expf_3 = p_stag_2/(((gamma-1)*(M_expf_3)^2/2 + 1)^(gamma/(gamma-1)));
%after first expansion reflection with symmetric first expansion fan
v = fzero(@(v) myfun(v,gamma,2*def_par),2);
M_expf1_1 = v;
p_expf1_1 = p_stag_2/(((gamma-1)*(M_expf1_1)^2/2 + 1)^(gamma/(gamma-1)));
g = fzero(@(g) myfun(g,gamma,3*def_par),2);
M_expf2_1 = g;
p_expf2_1 = p_stag_2/(((gamma-1)*(M_expf2_1)^2/2 + 1)^(gamma/(gamma-1)));
h = fzero(@(h) myfun(h,gamma,4*def_par),2);
M_expf3_1 = h;
p_expf3_1 = p_stag_2/(((gamma-1)*(M_expf3_1)^2/2 + 1)^(gamma/(gamma-1)));
i = fzero(@(g) myfun(g,gamma,4*def_par),3);
M_expf2_1_2 = i;
p_expf2_1_2 = p_stag_2/(((gamma-1)*(M_expf2_1_2)^2/2 + 1)^(gamma/(gamma-1)));
j = fzero(@(j) myfun(j,gamma,5*def_par),3);
M_expf2_1_2_3 = j;
p_expf2_1_2_3 = p_stag_2/(((gamma-1)*(M_expf2_1_2_3)^2/2 + 1)^(gamma/(gamma-1)));
k = fzero(@(k) myfun(k,gamma,6*def_par),3);
M_expf2_1_2_3_3 = k;
p_expf2_1_2_3_3 = p_stag_2/(((gamma-1)*(M_expf2_1_2_3_3)^2/2 + 1)^(gamma/(gamma-1)));
% we will consider 4 expansion fans including the first expansion fan
% we will get so many points where u have to calculate pressure at the back
% of the projectile it is divided into 10 regions and distribution is
% symmetric about centre. we will consider points where the expansion fans
% intersection is on the projectile so we get such 3 points
Reg_1_p = 12;
Reg_2_p = p_expf_1;
Reg_3_p = p_expf_2;
Reg_7_p = p_expf_3;
Reg_4_p = p_expf1_1;
Reg_5_p = p_expf2_1;
Reg_8_p = p_expf3_1;
Reg_6_p = p_expf2_1_2;
Reg_9_p = p_expf2_1_2_3;
Reg_10_p = p_expf2_1_2_3_3;



%Point_1(Reg_1,Reg_2,Reg_3,Reg_7)
reg2h = 12*(12*tand(2*def_par)-12*tand(def_par))/(12*tand(2*def_par));%mm
reg3h = 12*(12*tand(3*def_par)-12*tand(def_par))/(12*tand(3*def_par))-reg2h;%mm
reg7h = 12 - reg2h -reg3h;
t = input("Enter thickness of projectile:");
load_b_p1 = ((Reg_2_p*reg2h*t)+(Reg_3_p*reg3h*t)+(Reg_7_p*reg7h*t))*2/10;%N

%Point_2(Reg_4,3,7)
reg4h2 = ((6*tand(2*def_par))-(6*tand(def_par)))*tand(90-def_par);
reg3h2 = 12*(12*tand(3*def_par)-12*tand(def_par))/(12*tand(3*def_par))-reg4h2;%mm
reg7h2 = 12 - reg4h2 - reg3h2;
load_b_p2 = ((Reg_4_p*reg4h2*t)+(Reg_3_p*reg3h2*t)+(Reg_7_p*reg7h2*t))*2/10;%N

%Point_3(Reg_8,5)
reg5h3 = 12*(12*tand(3*def_par)-12*tand(2*def_par))/(12*tand(2*def_par));%mm
reg8h3 = tand(90-def_par)*(12*tand(2*def_par))-12-reg5h3+(12-reg5h3-reg8h3);%mm
load_b_p3 = ((Reg_5_p*reg5h3*t)+(Reg_8_p*reg8h3*t))*2/10;%N

%point_4(Reg_6,8)
reg6h4 = ((6*tand(3*def_par))-(6*tand(2*def_par)))*tand(90-(2*def_par));
reg8h4 = 12 - reg6h4;
load_b_p4 = ((Reg_6_p*reg6h4*t)+(Reg_8_p*reg8h4*t))*2/10;%N

%point_5(Reg_9,8)
reg9h5 = (tand(3*def_par)*12)*tand(90-(2*def_par))-12;
reg8h5 = 12-reg9h5;
load_b_p5 = ((Reg_9_p*reg9h5*t)+(Reg_8_p*reg8h5*t))*2/10;%N

%point_6(reg_9,10,8)
reg10h6 = (tand(90-(3*def_par)))*(1.5+(12*tand(3*def_par)))-12;
reg9h6 = 12 - reg10h6;
load_b_p6 = ((Reg_9_p*reg9h6*t)+(Reg_10_p*reg10h6*t))*2/10;%N

%Now consider the front parts of the projectile consider different points
area = 24*t;%mm^2
load_f_1 = (p_e*cosd(60)*2*area)/10;
load_f_2 = (p_b*cosd(60)*2*area)/10;
load_f_3 = (p_2*cosd(60)*2*area)/10;
load_f_4 = (p_4*cosd(60)*2*area)/10;

%Now calculate load acting at differnt points and combinations
load_init = (p_pro_b*area)/10 - load_f_1;
load_just = (p_pro_b*area)/10 - load_f_2;
load_f = [load_f_1,load_f_2,load_f_3,load_f_4];
load_b = [load_b_p1,load_b_p2,load_b_p3,load_b_p4,load_b_p5,load_b_p6];
figure(1)
for m = 1:length(load_b)
    load1 = load_b(m) - load_f_1;
    load2 = load_b(m) - load_f_2;
    load3 = load_b(m) - load_f_3;
    load4 = load_b(m) - load_f_4;
    x = [1,2,3,4,5,6];
    plot(x(m),load1,'o'); hold on;
    plot(x(m),load2,'o'); hold on;
    plot(x(m),load3,'o'); hold on;
    plot(x(m),load4,'o'); hold on;
end
hold off
grid
figure(2)
%plot x vs p alongcental line 
x = [1,2,3,4];
Reg_p = [Reg_1_p,Reg_4_p,Reg_6_p,Reg_10_p]
for m=1:length(x)
    p = Reg_p(m);
    x = [1,2,3,4];
    plot(x(m),p,'o'); hold 'on';
    xlabel('X along central line');
    ylabel('Pressure variations');
end
hold off
grid

figure(3)
%for x = 1
y = [1,2,3,4];
Reg1 = [Reg_1_p,Reg_2_p,Reg_3_p,Reg_7_p];
for m=1:length(y)
    p = Reg1(m);
    plot(y(m),p,'o'); hold 'on';
    xlabel('Y-locations for point_1');
    ylabel('Pressure variations');
    
end
hold off
grid
figure(4)
%for x = 3
y = [1,2,3];
Reg3 = [Reg_4_p,Reg_5_p,Reg_7_p];
for m=1:length(y)
    p = Reg3(m);
    plot(y(m),p,'o'); hold 'on';
    xlabel('Y-locations for point_3');
    ylabel('Pressure variations');
    
end
hold off
grid

figure(5)
%for x = 5
y = [1,2];
Reg5 = [Reg_5_p,Reg_8_p];
for m=1:length(y)
    p = Reg5(m);
    plot(y(m),p,'o'); hold 'on';
    xlabel('Y-locations for point_5');
    ylabel('Pressure variations');
    
end
hold off
grid

figure(6)
%for x = 6
y = [1,2];
Reg6 = [Reg_10_p,Reg_9_p];
for m=1:length(y)
    p = Reg6(m);
    plot(y(m),p,'o'); hold 'on';
    xlabel('Y-locations for point_6');
    ylabel('Pressure variations');
    
end
hold off
grid
function f = myfun(x,gamma,nu_2)
f = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((x^2)-1))^0.5))-atand(((x^2)-1)^0.5) - nu_2;
end


