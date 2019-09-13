M0 = 1.1:0.1:2;
M1 = 2.2:0.2:3;
M2 = 3.5:0.5:7;
M3 = 8:2:14;
M = [M0,M1,M2,M3];
gamma = 1.4; 
figure
for m=1:length(M) 
    beta = asind(1/M(m)):90;
    num =2*cotd(beta).*((M(m)^2)*(sind(beta).^2)-1);
    den = (M(m)^2)*(gamma + cosd(2*beta)) + 2;
    req = num./den;
    theta = atand(req);
    maxi(m) = find(theta==max(theta))
    A(m,1) = theta(maxi(m));
    B(m,1) = beta(maxi(m));
    plot(theta,beta,A,B); hold on
    axis([0 90 0 90]);
end
hold off