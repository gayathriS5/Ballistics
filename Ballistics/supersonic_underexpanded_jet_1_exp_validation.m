clc;
clear all;
global p_e p_b theta D x_co_f ang gamma
global M_blast P_atm P_e theta_f

%M_blast = input('Enter the Mach number of the blast wave initially coming out of the tube:');
P_atm = input('Enter the ambient pressure into which the blast wave is entering into:');
%P_e = input('Enter the pressure at the exit of tube from experiment(in bar)in the front case:');
gamma = input('Enter gamma:');
D = input('Enter Diameter of the tube :');
ang = input('Enter the angle of the projectile');
n_fp = input('Enter number of expansion fans excluding first exp_fan: ');
p_e = input('Enter the pressure at the exit of tube from experiment just behind projectile(in bar):');
np = input('Enter number of expansion fans u need to proceed excluding the first expansion fan:');
%M_e = 1;
M_e = input('Enter the mach number at the exit of tube backside:');
t = input('Enter thickness of projectile:');
length = input('Enter length of the projectile:');
mass = input('Enter the mass of projectile in kg:');
init_velocity = input('Enter initial velocity of the projectile in m/s:');
W_s = 725.6+583;
T_1 = 300;
v_1 = 0;
R = 8314/28.97;
M_f_1 = W_s/(gamma*R*T_1)^0.5;
T_2 = T_1*(2*gamma*(M_f_1)^2-(gamma-1))*(((gamma-1)*(M_f_1)^2)+2)/(((gamma+1)*M_f_1)^2);
M_f_2 = ((gamma-1)*(M_f_1^2) + 2)/((2*gamma*(M_f_1^2))-(gamma-1));
V_2 = W_s-(M_f_2*(gamma*R*T_2)^0.5);
M_f(1) = V_2/(gamma*R*T_2)^0.5;
P_e =(2*gamma*(M_f_1^2)/(gamma+1))-((gamma-1)/(gamma+1))*P_atm;
x_co_f = 0;
time = 0;
M_blast = M_f_1 - 0.5;
trend = 0;
while M_blast>=1.123
    n_f = n_fp;
    n = np;
    trend = trend+1;
    y_f = [];
    h_f = [];
    w_f = [];
    ld_f = [];
    P_f = [];
    P = [];
    points_x_f = [];
    points_y_f = [];
    M_f = [];
    M = [];
    pr_f = [];
    pr = [];
    exp_fan = [];
    exp_fan_f = [];
    beta = [];
    beta_f = [];
    M_f(1) = V_2/(gamma*R*T_2)^0.5;
    P_b = (2*gamma*(M_blast^2)/(gamma+1))-((gamma-1)/(gamma+1))*P_atm;
    p_b = P_b;
    %before expansion fan starts just at the end of tube
    %lets consider that the flow is just choked that is M=1
    i = 1;
    pr_f(i) = P_e;
    p_stag_f = P_e*(((gamma-1)*(M_f(1))^2/2 + 1)^(gamma/(gamma-1)));
    M_b_1 = ((((p_stag_f/P_b)^((gamma-1)/gamma)-1))*(2/(gamma-1)))^0.5;
    mu0 = asind(1/M_f(1));
    mu1 = asind(1/M_b_1);
    nu_0 = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_f(1)^2)-1))^0.5))-atand(((M_f(1)^2)-1)^0.5);
    nu_1 = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_b_1^2)-1))^0.5))-atand(((M_b_1^2)-1)^0.5);
    theta_f = nu_1 - nu_0;
    exp_fan_ft = mu0+(theta_f)-mu1;
    def_par_f = (theta_f)/n_f;
    while 1
        def_par_f = theta_f/n_f;
        defmin = 1.25;
        if def_par_f>=defmin
            break;
        end
        n_f = n_f -1;
    end
    %after expansion fans positive are crossed what are the properties of gas
    for s_f=2:n_f+1
        nu_f(1) = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_f(1)^2)-1))^0.5))-atand(((M_f(1)^2)-1)^0.5);
        nu_f(s_f) = nu_f(s_f-1)+def_par_f;
        x = fzero(@(x) prandtl(x,nu_f(s_f)),M_f(s_f-1)+1);
        M_f(s_f) = x;
    end
    for e = 1:n_f+1
        beta_f(e) = asind(1/M_f(e));
    end
    %if exp_fan angle is greater than 90deg then it wont involve in mesh and it
    %will be above the projectile which effectively get cancelled so no need to
    %bother about it
    f_f = 0;
    for e=1:n_f
        exp_fan_f(e) = (90-beta_f(1))+beta_f(1)+(e*def_par_f)-beta_f(e+1);
        if exp_fan_f(e)>=90
            f_f = f_f+1;
        end
    end
    n_f = n_f-f_f;
    tot_regions_f = ((n_f+2)*(n_f+1))/2;
    tot_rows_f = n_f+1;
    tot_points_f =(n_f*(n_f+1)/2) + n_f + 2;
    P_f = ones(tot_points_f,2);
    r_f = ones(tot_points_f,1);
    k_f = 0;
    i = tot_rows_f;
    %Allocating mach numbers of different regions
    %mach series is not along x axis but along the lines of expfans
    %Now allocating pressure to those regions allotted to mach numbers
    % M denotes mach number and pr denotes pressure

    s_f = 1;
    i = tot_rows_f;
    while i>=1
        for j_f = 1:i-1
            s_f = s_f+1;
            pr_f(1) = P_e;
            nu_f(1) = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_f(1)^2)-1))^0.5))-atand(((M_f(1)^2)-1)^0.5);
            nu_f(s_f) = nu_f(s_f-1)+def_par_f;
            x = fzero(@(x) prandtl(x,nu_f(s_f)),M_f(s_f-1)+1);
            M_f(s_f) = x;
            pr_f(s_f) = p_stag_f/(((gamma-1)*(M_f(s_f))^2/2 + 1)^(gamma/(gamma-1)));
        end
        i = i-1;
        if i>=1
            s_f = s_f+1;
            r_f = s_f-i;
            nu_f(s_f) = nu_f(r_f) + def_par_f;
            x = fzero(@(x) prandtl(x,nu_f(s_f)),M_f(r_f)+1);
            M_f(s_f) = x;
            pr_f(s_f) = p_stag_f/(((gamma-1)*(M_f(s_f))^2/2 + 1)^(gamma/(gamma-1)));

        end
    end

    %Allocating an array of slopes for both negative and positive expansions
    for o = 1:n_f
        slopep_f(o) = tand(90+(exp_fan_f(o)));
        slopen_f(o) = tand(90-(exp_fan_f(o)));
    end
    %Finding all mesh points of intersections of expansion fans behind the
    %projectile and P() is a matrix with both x and y coordinates 
    % and r denotes an array containing the rows of the mesh

    i = n_f+1;
    while i>1
        i = i-1;
        for j_f = 1:i
            k_f = k_f+1;
            r_f(k_f,1) = i;
            d = abs(i-(tot_rows_f-1));
            if i == n_f
                P_f(k_f,1) = D/2*(tand(exp_fan_f(j_f)));
                P_f(k_f,2) = 0;
            else
                P_f(k_f,2) = (P_f(j_f+d,1)-P_f(j_f,1))/(cotd(90-(exp_fan_f(j_f)))+cotd(90-exp_fan_f(d+j_f)));
                P_f(k_f,1) = (D/2 + P_f(k_f,2))/tand(90-(exp_fan_f(j_f)));
            end

        end
    end
    i = n_f+1;
    for l_f = 1:i-1
        k_f = k_f+1;
        r_f(k_f,1) = i;
        P_f(k_f,1) = D/tand(90-(exp_fan_f(l_f)));
        P_f(k_f,2) = D/2;
    end
    r_f(k_f+1) = i;
    P_f(k_f+1,1) = 0;
    P_f(k_f+1,2) = D/2;
    r_f(k_f+2,1) = i-1;
    P_f(k_f+2,1) = 0;
    P_f(k_f+2,2) = 0;
    %include even oblique shocks backside the load only for last expansion fan
    %before u have to consider the expansion fans greater than 90deg
    %find their pressures and mach numbers in that region
    b_f(1) = asind(1/M_f(tot_regions_f));
    a_f(1) = P_f((n_f),1);
    c_f(1) = P_f((n_f),2);
    for q_f = 1:f_f
        angle_f(q_f) = exp_fan_f(n_f+q_f) - 90;
        x = fsolve(@(x) ntng(x,angle_f(q_f),b_f(q_f)+(q_f-1)*def_par_f,a_f(q_f),c_f(q_f)),[0 0]);
        P_f((tot_points_f+q_f),1) = x(1);
        P_f((tot_points_f+q_f),2) = x(2);
        a_f(q_f+1) = x(1);
        c_f(q_f+1) = x(2); 
        nu_f(tot_regions_f+q_f) = nu_f(tot_regions_f+q_f-1)+ def_par_f;
        x = fzero(@(x) prandtl(x,nu_f(tot_regions_f+q_f)),M_f(tot_regions_f+q_f-1)+1);
        M_f(tot_regions_f+q_f) = x;
        pr_f(tot_regions_f+q_f) = p_stag_f/((gamma-1)*(M_f(tot_regions_f+q_f))^2/2 + 1)^(gamma/(gamma-1));
        b_f(q_f+1) = asind(1/M_f(tot_regions_f+q_f));
    end
    if f_f==0
        q_f=0;
    end
    a_f = n_f+f_f ;
    hz_f = 0;
    g_f = 1;
    while a_f>f_f
        angles_f = a_f*def_par_f;
        for z_f = g_f:n_f-1
            if 90-exp_fan_f(z_f)>angles_f
                hz_f = hz_f+1;
                if a_f == n_f+f_f
                    x = fsolve(@(x) tem(x,angles_f,90-exp_fan_f(z_f)),[0 0]);
                    points_x_f(hz_f) = x(1);
                    points_y_f(hz_f) = x(2);
                else
                    x = fsolve(@(x) temp(x,angles_f,points_x_f(hz_f-1),points_y_f(hz_f-1),90-exp_fan_f(z_f)),[0 0]);
                    points_x_f(hz_f) = x(1);
                    points_y_f(hz_f) = x(2);
                end
                a_f = n_f+f_f-z_f;
                g_f = z_f+1;
            else
                d = abs(90-exp_fan_f(z_f)-angles_f);
                ins = ceil(d+90-exp_fan_f(z_f));
                hz_f = hz_f+1;
                if a_f == n_f+f_f
                    x = fsolve(@(x) tem(x,angles_f,ins),[0 0]);
                    points_x_f(hz_f) = x(1);
                    points_y_f(hz_f) = x(2);
                else
                    x = fsolve(@(x) temp(x,angles_f,points_x_f(hz_f-1),points_y_f(hz_f-1),ins),[0 0]);
                    points_x_f(hz_f) = x(1);
                    points_y_f(hz_f) = x(2);
                end
                a_f = n_f+f_f-z_f;
                g_f = z_f+1;
            end
            
        end
        a_f = n_f+f_f-z_f;
        g_f = z_f+1;
    end
    if n_f >1
        if f_f 
            x = fsolve(@(x) final(x,(n_f+f_f-hz_f)*def_par_f,points_x_f(hz_f),points_y_f(hz_f),b_f(f_f+1)+(f_f*def_par_f),P_f(tot_points_f+f_f,1),P_f(tot_points_f+f_f,2)),[0 0]);
            P_f((tot_points_f+f_f+1),1) = x(1);
            P_f((tot_points_f+f_f+1),2) = x(2);
        else
            x = fsolve(@(x) final(x,(n_f+f_f-hz_f)*def_par_f,points_x_f(hz_f),points_y_f(hz_f),asind(1/M_f(tot_regions_f)),P_f(n_f,1),P_f(n_f,2)),[0 0]);
            P_f((tot_points_f+f_f+1),1) = x(1);
            P_f((tot_points_f+f_f+1),2) = x(2);
        end
    end
    if n_f == 1
        x = fsolve(@(x) ntng(x,(n_f+f_f-hz_f)*def_par_f,b_f(f_f+1)+(f_f*def_par_f),P_f(tot_points_f+f_f,1),P_f(tot_points_f+f_f,2)),[0 0]);
        P_f((tot_points_f+f_f+1),1) = x(1);
        P_f((tot_points_f+f_f+1),2) = x(2);
    end
    M_obq_n_f = (((P_b/pr_f(tot_regions_f+q_f))*(gamma+1) + (gamma-1))/(2*gamma))^0.5;
    beta_obq_f = asind(M_obq_n_f/M_f(tot_regions_f+q_f));
    num =2*cotd(beta_obq_f)*((M_f(tot_regions_f+q_f)^2)*(sind(beta_obq_f)^2)-1);
    den = (M_f(tot_regions_f+q_f)^2)*(gamma + cosd(2*beta_obq_f)) + 2;
    req = num/den;
    theta_obq_f = atand(req);
    M_af_n_f = ((gamma-1)*(M_obq_n_f^2) + 2)/((2*gamma*(M_obq_n_f^2))-(gamma-1));
    M_af_f = M_af_n_f/sind(beta_obq_f - theta_obq_f); %after oblique incident shock
    x = fsolve(@(x) fun(x,abs(beta_obq_f-f_f*def_par_f),P_f(tot_points_f+f_f+1,1),P_f(tot_points_f+f_f+1,2)),[0 0]);
    P_f((tot_points_f+f_f+2),1) = x(1);
    P_f((tot_points_f+f_f+2),2) = x(2);
    %consider the height of mach stem is equal to diameter of pipe
    %Now consider the properies across the normal shock
    M_f(tot_regions_f+q_f+1) = ((gamma-1)*(M_f(tot_regions_f+q_f)^2) + 2)/((2*gamma*(M_f(tot_regions_f+q_f)^2))-(gamma-1));%Normal shock
    Prat_f = (2*gamma*(M_f(tot_regions_f+q_f)^2)/(gamma+1))-((gamma-1)/(gamma+1));
    pr_f(tot_regions_f+q_f+1) = Prat_f*pr_f(tot_regions_f);
    %Npw backside the projectile
    p_stag = p_e*(((gamma-1)*(M_e)^2/2 + 1)^(gamma/(gamma-1)));
    M_b = ((((p_stag/p_b)^((gamma-1)/gamma)-1))*(2/(gamma-1)))^0.5;
    mu_e = asind(1/M_e);
    mu_b = asind(1/M_b);
    nu_e = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_e^2)-1))^0.5))-atand(((M_e^2)-1)^0.5);
    nu_b = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M_b^2)-1))^0.5))-atand(((M_b^2)-1)^0.5);
    theta = nu_b - nu_e;
    def_par = theta/n;
    while 1
        def_par = theta/n;
        defmin = 1.25;
        if def_par>=defmin
            break;
        end
        n = n -1;
    end
    %Now calculate all possible regions in expansion fan jets and find
    %pressures in those regions neglect exp fan intersections only consider th
    %this is to find points and intersection of fans
    % nu means prandtl angle and mu is mach angle
    for s=2:n+1
        M(1) = M_e;
        nu(1) = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M(1)^2)-1))^0.5))-atand(((M(1)^2)-1)^0.5);
        nu(s) = nu(s-1)+def_par;
        x = fzero(@(x) prandtl(x,nu(s)),M(s-1)+1);
        M(s) = x;
    end

    for e = 1:n+1
        beta(e) = asind(1/M(e));
    end
    %if exp_fan angle is greater than 90deg then it wont involve in mesh and it
    %will be above the projectile which effectively get cancelled so no need to
    %bother about it
    f = 0;
    for e=1:n
        exp_fan(e) = (90-beta(1))+beta(1)+(e*def_par)-beta(e+1);
        if exp_fan(e)>=90
            f = f+1;
        end
    end
    n = n-f;
    tot_regions = ((n+2)*(n+1))/2;
    tot_rows = n+1;
    tot_points =(n*(n+1)/2) + n + 2;
    P = ones(tot_points,2);
    r = ones(tot_points,1);
    k = 0;
    i = tot_rows;
    %Allocating mach numbers of different regions
    %mach series is not along x axis but along the lines of expfans
    %Now allocating pressure to those regions allotted to mach numbers
    % M denotes mach number and pr denotes pressure
    p_stag = p_e*(((gamma-1)*(M_e)^2/2 + 1)^(gamma/(gamma-1)));
    s = 1;
    i = tot_rows;
    while i>=1
        for j = 1:i-1
            s = s+1;
            M(1) = M_e;
            pr(1) = p_e;
            nu(1) = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((M(1)^2)-1))^0.5))-atand(((M(1)^2)-1)^0.5);
            nu(s) = nu(s-1)+def_par;
            x = fzero(@(x) prandtl(x,nu(s)),M(s-1)+1);
            M(s) = x;
            pr(s) = p_stag/(((gamma-1)*(M(s))^2/2 + 1)^(gamma/(gamma-1)));
        end
        i = i-1;
        if i>=1
            s = s+1;
            r = s-i;
            nu(s) = nu(r) + def_par;
            x = fzero(@(x) prandtl(x,nu(s)),M(r)+1);
            M(s) = x;
            pr(s) = p_stag/(((gamma-1)*(M(s))^2/2 + 1)^(gamma/(gamma-1)));

        end
    end

    %Allocating an array of slopes for both negative and positive expansions
    for o = 1:n
        slopep(o) = tand(90+(exp_fan(o)));
        slopen(o) = tand(90-(exp_fan(o)));
    end
    %Finding all mesh points of intersections of expansion fans behind the
    %projectile and P() is a matrix with both x and y coordinates 
    % and r denotes an array containing the rows of the mesh

    i = n+1;
    while i>1
        i = i-1;
        for j = 1:i
            k = k+1;
            r(k,1) = i;
            d = abs(i-(tot_rows-1));
            if i == n
                P(k,1) = D/2*(tand(exp_fan(j)));
                P(k,2) = 0;
            else
                P(k,2) = (P(j+d,1)-P(j,1))/(cotd(90-(exp_fan(j)))+cotd(90-exp_fan(d+j)));
                P(k,1) = (D/2 + P(k,2))/tand(90-(exp_fan(j)));
            end

        end
    end
    i = n+1;
    for l = 1:i-1
        k = k+1;
        r(k,1) = i;
        P(k,1) = D/tand(90-(exp_fan(l)));
        P(k,2) = D/2;
    end
    r(k+1) = i;
    P(k+1,1) = 0;
    P(k+1,2) = D/2;
    r(k+2,1) = i-1;
    P(k+2,1) = 0;
    P(k+2,2) = 0;
    %include even oblique shocks backside the load only for last expansion fan
    %before u have to consider the expansion fans greater than 90deg
    %find their pressures and mach numbers in that region
    b(1) = asind(1/M(tot_regions));
    a(1) = P((n),1);
    c(1) = P((n),2);
    for q = 1:f
        angle(q) = exp_fan(n+q) - 90;
        x = fsolve(@(x) ntng(x,angle(q),b(q)+(q-1)*def_par,a(q),c(q)),[0 0]);
        P((tot_points+q),1) = x(1);
        P((tot_points+q),2) = x(2);
        a(q+1) = x(1);
        c(q+1) = x(2); 
        nu(tot_regions+q) = nu(tot_regions+q-1)+ def_par;
        x = fzero(@(x) prandtl(x,nu(tot_regions+q)),M(tot_regions+q-1)+1);
        M(tot_regions+q) = x;
        pr(tot_regions+q) = p_stag/((gamma-1)*(M(tot_regions+q))^2/2 + 1)^(gamma/(gamma-1));
        b(q+1) = asind(1/M(tot_regions+q));
    end
    if f==0
        q=0;
    end
    a = n+f ;
    hz = 0;
    g = 1;
    while a>f
        angles = a*def_par;
        for z = g:n-1
            if 90-exp_fan(z)>angles
                hz = hz+1;
                if a == n+f
                    x = fsolve(@(x) tem(x,angles,90-exp_fan(z)),[0 0]);
                    points_x(hz) = x(1);
                    points_y(hz) = x(2);
                else
                    x = fsolve(@(x) temp(x,angles,points_x(hz-1),points_y(hz-1),90-exp_fan(z)),[0 0]);
                    points_x(hz) = x(1);
                    points_y(hz) = x(2);
                end
                a = n+f-z;
                g = z+1;
            else
                d = abs(angles-(90-exp_fan(z)));
                ins_b = ceil(90-exp_fan(z)+d);
                hz = hz+1;
                if a == n+f
                    x = fsolve(@(x) tem(x,angles,ins_b),[0 0]);
                    points_x(hz) = x(1);
                    points_y(hz) = x(2);
                else
                    x = fsolve(@(x) temp(x,angles,points_x(hz-1),points_y(hz-1),ins_b),[0 0]);
                    points_x(hz) = x(1);
                    points_y(hz) = x(2);
                end
                a = n+f-z;
                g = z+1;
            end
        end
        a = n+f-z;
        g = z+1;
    end
    if f
        x = fsolve(@(x) final(x,(n+f-hz)*def_par,points_x(hz),points_y(hz),b(f+1)+(f*def_par),P(tot_points+f,1),P(tot_points+f,2)),[0 0]);
        P((tot_points+f+1),1) = x(1);
        P((tot_points+f+1),2) = x(2);
    else
        x = fsolve(@(x) final(x,(n+f-hz)*def_par,points_x(hz),points_y(hz),asind(1/M(tot_regions)),P(n,1),P(n,2)),[0 0]);
        P((tot_points+f+1),1) = x(1);
        P((tot_points+f+1),2) = x(2);
    end
    M_obq_n = (((P_b/pr(tot_regions+q))*(gamma+1) + (gamma-1))/(2*gamma))^0.5;
    beta_obq = asind(M_obq_n/M(tot_regions+q));
    num =2*cotd(beta_obq)*((M(tot_regions+q)^2)*(sind(beta_obq)^2)-1);
    den = (M(tot_regions+q)^2)*(gamma + cosd(2*beta_obq)) + 2;
    req = num/den;
    theta_obq = atand(req);
    M_af_n = ((gamma-1)*(M_obq_n^2) + 2)/((2*gamma*(M_obq_n^2))-(gamma-1));
    M_af = M_af_n/sind(beta_obq - theta_obq); %after oblique incident shock
    x = fsolve(@(x) fun(x,abs(beta_obq-f*def_par),P(tot_points+f+1,1),P(tot_points+f+1,2)),[0 0]);
    P((tot_points+f+2),1) = x(1);
    P((tot_points+f+2),2) = x(2);
    %consider the height of mach stem is equal to diameter of pipe
    %Now consider the properies across the normal shock
    M(tot_regions+q+1) = ((gamma-1)*(M(tot_regions+q)^2) + 2)/((2*gamma*(M(tot_regions+q)^2))-(gamma-1));%Normal shock
    Prat_b = (2*gamma*(M(tot_regions+q)^2)/(gamma+1))-((gamma-1)/(gamma+1));
    pr(tot_regions+q+1) = Prat_b*pr(tot_regions);

    %Now as points and regions are defined its time to find forces
    %x_co_f = input('Enter the x-position of front of projectile from the end of the tube:');
    %x_co_f DENOTES THE X POSITION of front most part of projectile ALONG central projectile
    %x_co denotes x position of backside of the projectile
    %projectile is 50mm long 24 mm diameter and 24 mm thickness with
    %symmetrical beak of angle 30deg
    %figure() 
    x_co = x_co_f - length ;
    dome = cotd(ang)*D/2;
    if x_co_f<=0
        tot_load_f = P_e*D*t/10;
    end
    flip = 0;
    if x_co_f<=P_f(1,1) && x_co_f>0
        for i=1:n_f
            x = fsolve(@(x) fan(x,slopep_f(i)),[0 0])
            if x(2)<=D/2 && x(2)>0
                y_f(i) = x(2);
                flip = 1
                w_f(1) = y_f(1);
                if i>1
                    w_f(i) = y_f(i)-y_f(i-1);
                end
                ld_f(i) = pr_f(i)*w_f(i)*t;
            else
                break
            end
        end
        if flip 
            w_f(i+1) = D/2-y_f(i);
            ld_f(i+1) = pr_f(i+1)*w_f(i+1)*t;
            tot_load_f = (sum(ld_f))*2/10;
        else
            tot_load_f = P_e*D*t/10;
        end

    end

    r_f = n_f+1;
    r0_f = n_f+1;
    j_f = 0;
    s_f = 0;
    for k_f=1:n_f
        v_f =0;
        y_f = [];
        h_f = [];
        w_f = [];
        ld_f = [];
        u_f = [];
        z_f = [];
        g_f = [];
        ht_f = [];
        if (P_f(k_f,1)<=x_co_f && x_co_f<P_f(k_f+1,1) && k_f<n_f && x_co_f<dome)
            for v_f = 1:k_f
                x = fsolve(@(x) fan(x,slopep_f(v_f)),[0 0]);
                u_f(v_f) = x(2);
                if u_f(v_f)<=D/2 && u_f(v_f)>=0
                    j_f = j_f+1;
                    y_f(j_f) = u_f(v_f)
                end
            end
            for m = 1:k_f
                x = fsolve(@(x) fluke(x,slopen_f(m)),[0 0]);
                g_f(m) = x(2);
                if g_f(m)<=D/2 && g_f(m)>=0
                    s_f = s_f+1;
                    z_f(s_f) = g_f(m)
                end
            end
            ht_f = [y_f,z_f];
            tot = s_f+j_f;
            h_f = sort(ht_f,'ascend');
            w_f(1) = h_f(1);
            ld_f(1) = pr_f(r_f+1)*w_f(1)*t;
            indi = r_f+1;
            row = n_f+1-k_f;
            for i = 2:tot
                if ismember(h_f(i-1),z_f)
                    w_f(i) = h_f(i)-h_f(i-1);
                    ld_f(i) = pr_f(indi-row)*w_f(i)*t;
                    indi = indi-row;
                    row = row+1;
                else
                    w_f(i) = h_f(i)-h_f(i-1);
                    ld_f(i) = pr_f(indi-1)*w_f(i)*t;
                    indi = indi-1;                  
                end
            end
            w_f(i+1) = D/2 - h_f(tot);
            ld_f(i+1) = pr_f(1)*w_f(i+1)*t;
            tot_load_f = sum(ld_f)*2/10;
            break
        end
            
        if (P_f(k_f,1)<=x_co_f && x_co_f<P_f(tot_points_f-2,1) && k_f==n_f && x_co_f<dome)
            for v_f = 1:k_f
                x = fsolve(@(x) fan(x,slopep_f(v_f)),[0 0]);
                u_f(v_f) = x(2);
                if u_f(v_f)<=D/2 && u_f(v_f)>=0
                    j_f = j_f+1;
                    y_f(j_f) = u_f(v_f)
                end
            end
            for m = 1:k_f
                x = fsolve(@(x) fluke(x,slopen_f(m)),[0 0]);
                g_f(m) = x(2);
                if g_f(m)<=D/2 && g_f(m)>=0
                    s_f = s_f+1;
                    z_f(s_f) = g_f(m)
                end
            end
            ht_f = [y_f,z_f];
            tot = s_f+j_f;
            h_f = sort(ht_f,'ascend');
            w_f(1) = h_f(1);
            ld_f(1) = pr_f(r_f+1)*w_f(1)*t;
            indi = r_f+1;
            row = n_f+1-k_f;
            for i = 2:tot
                if ismember(h_f(i-1),z_f)
                    w_f(i) = h_f(i)-h_f(i-1);
                    ld_f(i) = pr_f(indi-row)*w_f(i)*t;
                    indi = indi-row;
                    row = row+1;
                else
                    w_f(i) = h_f(i)-h_f(i-1);
                    ld_f(i) = pr_f(indi-1)*w_f(i)*t;
                    indi = indi-1;                  
                end
            end
            w_f(i+1) = D/2 - h_f(tot);
            ld_f(i+1) = pr_f(1)*w_f(i+1)*t;
            tot_load_f = sum(ld_f)*2/10;
            break
        end
        if (P_f(k_f,1)<=x_co_f && x_co_f<P_f(k_f+1,1) && k_f<n_f && x_co_f>=dome)
            for v_f = 1:n_f-k_f
                x = fsolve(@(x) fan(x,slopep_f(v_f+k_f)),[0 0]);
                u_f(v_f) = x(2);
                if u_f(v_f)<=D/2 && u_f(v_f)>=0
                    j_f = j_f+1;
                    y_f(j_f) = u_f(v_f);
                end
            end
            for m = 1:k_f
                x = fsolve(@(x) fluke(x,slopen_f(m)),[0 0]);
                g_f(m) = x(2);
                if g_f(m)<=D/2 && g_f(m)>=0
                    s_f = s_f+1;
                    z_f(s_f) = g_f(m)
                end
            end
            ht_f = [y_f,z_f];
            tot = s_f+j_f;
            h_f = sort(ht_f,'ascend');
            w_f(1) = h_f(1);
            ld_f(1) = pr_f(r_f+1)*w_f(1)*t;
            indi = r_f+1;
            row = n_f+1-k_f;
            i = 1;
            for i = 2:tot
                if ismember(h_f(i-1),z_f)
                    w_f(i) = h_f(i)-h_f(i-1);
                    ld_f(i) = pr_f(indi-row)*w_f(i)*t;
                    indi = indi-row;
                    row = row+1;
                else
                    w_f(i) = h_f(i)-h_f(i-1);
                    ld_f(i) = pr_f(indi+1)*w_f(i)*t;
                    indi = indi+1;                  
                end
            end
            w_f(i+1) = D/2 - h_f(tot);
            if s_f
                if h_f(tot) == z_f(s_f)
                    ld_f(i+1) = pr_f(indi-row)*w_f(i+1)*t;
                else
                    ld_f(i+1) = pr_f(indi+1)*w_f(i+1)*t;
                end
            else
                ld_f(i+1) = pr_f(indi+1)*w_f(i+1)*t;
            end
            tot_load_f = sum(ld_f)*2/10;            
            break
        end

        if (P_f(k_f,1)<=x_co_f && x_co_f<P_f(tot_points_f-2,1) && k_f==n_f && x_co_f>=dome)
            for m = 1:k_f
                x = fsolve(@(x) fluke(x,slopen_f(m)),[0 0]);
                g_f(m) = x(2);
                if g_f(m)<=D/2 && g_f(m)>=0
                    s_f = s_f+1;
                    z_f(s_f) = g_f(m);
                end
            end
            h_f = sort(z_f,'ascend');
            w_f(1) = h_f(1);
            ld_f(1) = pr_f(r_f+1)*w_f(1)*t;
            indi = r_f+1;
            row = n_f+1-k_f;
            for i = 2:s_f
                w_f(i) = h_f(i)-h_f(i-1);
                ld_f(i) = pr_f(indi-row)*w_f(i)*t;
                indi = indi-row;
                row = row+1;
            end
            w_f(i+1) = D/2 - h_f(s_f);
            ld_f(i+1) = pr_f(indi-row)*w_f(i+1)*t;
            tot_load_f = sum(ld_f)*2/10;
            break
        end
    r_f = r_f+r0_f-k_f;
    end

    if x_co_f >= P_f(tot_points_f-2,1) && x_co_f<=P_f(tot_points_f+f_f+2,1)
        tot_load_f = (pr_f(tot_regions_f)*D*t)/10;

    end
    if x_co_f>=P_f(tot_points_f+f_f+2)
        y_f = tand(180-ang)*(P_f(tot_points_f+f_f+2)-x_co_f);
        if y_f>=0 && y_f<= D/2
            w_f = y_f;
            load_f(1) = pr_f(tot_regions_f+q_f+1)*w_f*t;
            load_f(2) = pr_f(tot_regions_f)*(D/2-w_f)*t;
            tot_load_f = (load_f(1)+load_f(2))*2/10;
        else
            tot_load_f = (pr_f(tot_regions_f+q_f+1)*D*t)/10;
        end
    end
    if x_co<=0
        tot_load_b = p_e*D*t/10;
    end

    if x_co<=P(1,1) && x_co>0
        for i=1:n
            y(i) = (slopep(i)*x_co)+P(tot_points-1,2);
            w(1) = y(1);
            if i>1
                w(i) = y(i)-y(i-1);
            end
            ld(i) = pr(i)*w(i)*t;

        end
        w(i+1) = D/2-y(i);
        ld(i+1) = pr(i+1)*w(i+1)*t;
        tot_load_b = (sum(ld))*2/10;

    end

    r = n+1;
    r0 = n+1;
    l = 0;
    tup = 0;
    j = 0;
    for k=1:n
        v =0;
        y = [];
        h = [];
        w = [];
        ld = [];
        u = [];
        z = [];
        g = [];
        if (P(k,1)<=x_co && x_co<P(k+1,1) && k<n)
            for v = 1:n-k
                u(v) = (slopep(v+k)*x_co)+P(tot_points-1,2);
                if u(v)<=12 && u(v)>=0
                    j = j+1;
                    y(j) = u(v);
                end
            end
            u(v+1) = (slopen(k)*x_co)-P(tot_points-1,2);
            if u(v+1)<=12
                y(j+1) = u(v+1);
                tup = 1;
            end
            h = sort(y,'ascend');
            w(1) = h(1);
            ld(1) = pr(r+1)*w(1)*t;
            if tup
                result = find(h==y(j+1));
                for i = 1:result
                    if i>1
                        w(i) = h(i)-h(i-1);
                        ld(i) = pr(r+i)*w(i)*t;
                    end
                end
                for i=result+1:j+1
                    w(i) = h(i) - h(i-1);
                    ld(i) = pr(l+i)*w(i)*t;

                end
                s = 0;
                for m = 1:k-1
                    g(m) = (slopen(m)*x_co)-P(tot_points-1,2);
                    if g(m)<=12 && g(m)>=0
                        s = s+1;
                        z(s) = g(m);
                    end
                end
                p = 1;
                if s
                    i = j+1;
                    i = i+1;
                    ht = sort(z,'ascend');
                    w(i) = ht(p)-h(i-1);
                    ld(i) = pr(l+i)*w(i)*t;
                    flip = l+i;
                    for m=1:s-1
                        if s>1
                            i = i+1;
                            p = p+1;
                            w(i) = ht(p)-ht(p-1);
                            ld(i) = pr(flip-(n+1-k-(p-2)))*w(i)*t;
                            flip = flip-(n+1-k+(p-2));
                        end
                    end
                    i = i+1;
                    p = p+1;
                    w(i) = D/2-ht(s);
                    ld(i) = pr(flip-(n+1-k-(p-2)))*w(i)*t;
                    tot_load_b = (sum(ld));
                    tot_load_b = tot_load_b*2/10;

                else
                    i = j+1;
                    i = i+1;
                    w(i) = D/2-h(j+1);
                    ld(i) = pr(l+i)*w(i)*t;
                    tot_load_b = (sum(ld));
                    tot_load_b = tot_load_b*2/10;

                end
            else
                for i = 1:j
                    if i>1
                        w(i) = h(i)-h(i-1);
                        ld(i) = pr(r+i)*w(i)*t;
                    end
                    w(i+1) = D/2 - w(i);
                    ld(i+1) = pr(r+i+1)*w(i+1)*t;
                    tot_load_b = sum(ld)*2/10;
                end
            end
            break
        end

        if (P(k,1)<=x_co && x_co<P(tot_points-2,1) && k==n)
            u(v+1) = (slopen(k)*x_co)-P(tot_points-1,2);
            if u(v+1)<=12
                y(j+1) = u(v+1);
            end
            ld(1) = pr(r+1)*y(1)*t;
            s = 0;
            for m = 1:k-1
                g(m) = (slopen(m)*x_co)-P(tot_points-1,2);
                if g(m)<=12 && g(m)>=0
                    s = s+1;
                    z(s) = g(m);
                end
            end
            p = 1;
            if s
                i = j+1;
                i = i+1;
                ht = sort(z,'ascend');
                w(i) = ht(p)-y(i-1);
                ld(i) = pr(l+i)*w(i)*t;
                flip = l+i;
                for m=1:s-1
                    if s>1
                        i = i+1;
                        p = p+1;
                        w(i) = ht(p)-ht(p-1);
                        ld(i) = pr(flip-(n+1-k-(p-2)))*w(i)*t;

                        flip = flip-(n+1-k+(p-2));
                    end
                end
                i = i+1;
                w(i) = D/2-ht(s);
                ld(i) = pr(flip-(n+1-k-(p-2)))*w(i)*t;
                tot_load_b = (sum(ld));
                tot_load_b = tot_load_b*2/10;

            else
                i = j+1;
                i = i+1;
                w(i) = D/2-y(j+1);
                ld(i) = pr(l+i)*w(i)*t;
                tot_load_b = (sum(ld));
                tot_load_b = tot_load_b*2/10;
                break

            end
        end
    l = r;
    r = r+r0-k;
    end

    if x_co >= P(tot_points-2) && x_co<=P(tot_points+f+2)
        tot_load_b = (pr(tot_regions)*D*t)/10;
    end
    if x_co>=P(tot_points+f+2)
        tot_load_b = (pr(tot_regions+q+1)*D*t)/10;
    end 
    net_load = tot_load_b - tot_load_f;
    acceleration = net_load/mass;
    dt = 1/62626;%(in sec);
    velocity = acceleration*dt + init_velocity;
    distance = (acceleration*(dt^2)/2) + (init_velocity*dt); %(initial position = 0) in m
    figure(1)
    plot(x_co_f,acceleration,'-o'); hold on;
    xlabel('time in seconds');
    ylabel('acceleration of projectile');
    figure(2)
    plot(time,P_b,'o');hold on;
    xlabel('time in seconds');
    ylabel('Pressure into which projectile is entering into');
    x_co_f = x_co_f+distance*1000;
    init_velocity = velocity;
    time = time+dt;
    blast_speed = 725.6*exp(-28.97*1000*(time))+583*exp(-0.753*1000*(time));
    T_1 = 300;
    v_1 = 0;
    R = 8314/28.97;
    if M_blast >=1.4
        P_e = P_e - 0.8;
        p_e = p_e - 0.2;
    else 
        P_e = P_e - 0.1;
        p_e = p_e -0.65;
    end
    M_blast = blast_speed/(gamma*R*T_1)^0.5;
end
hold off;
function f = prandtl(x,nu_2)
global gamma
f = ((((gamma+1)/(gamma-1))^0.5)*atand((((gamma-1)/(gamma+1))*((x^2)-1))^0.5))-atand(((x^2)-1)^0.5) - nu_2;
end
function g = myfun(x,beta_f,a)
global theta_f D 
g(1) = x(2)-(tand(theta_f)*x(1)+ (D/2));
g(2) = x(2)-(tand(beta_f)*(x(1)-a));
end
function h = fun(x,beta_1,a,b)
global D
h(1) = x(2)-b-(tand(180-beta_1)*(x(1)-a));
h(2) = x(2)-3*D/4;
end
function z = fan(x,slope)
global x_co_f D ang
z(1) = x(2)-(((slope)*x(1))+(D/2));
z(2) = x(2)-(tand(180-ang)*(x(1)-x_co_f));
end
function m = fluke(x,neg_slope)
global x_co_f ang D
m(1) = x(2)-(((neg_slope)*x(1))-D/2);
m(2) = x(2)-(tand(180-ang)*(x(1)-x_co_f));
end
function s = ntng(x,angle,b,a,c)
global D 
s(1) = x(2)-(tand(angle)*x(1)+ (D/2));
s(2) = x(2)-(tand(b)*(x(1)-a)+c);
end
function  o = final(x,d,e,f,b,a,c)
o(1) = x(2)-(tand(d)*(x(1)-e)+f);
o(2) = x(2)-(tand(b)*(x(1)-a)+c);
end
function t = tem(x,angle1,angle2)
global D
t(1) = x(2)-(tand(angle1)*x(1)+D/2);
t(2) = x(2)-(tand(angle2)*x(1)-D/2);
end
function u = temp(x,angle1,points_x,points_y,angle2) 
global D
u(1) = x(2)-(tand(angle1)*(x(1)-points_x)+points_y);
u(2) = x(2)-(tand(angle2)*x(1)-D/2);
end