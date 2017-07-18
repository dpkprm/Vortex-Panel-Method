
%% Data Acquisition and Re-arranging
data        = load('naca2410_160cos.txt'); %Load Airfoil Coordianted from http://airfoiltools.com/airfoil/naca4digit
XB          = flip(data(:,1));
YB          = flip(data(:,2));
YB(1,1)     = 0;
YB(end,1)    = 0;

M = size(XB,1);

z = 1;
%% Calculation of control points and other geometric parameters
alpha = 0;
alpha = (alpha) * (pi/180);

for i = 1:M-1
    X(i,1)        = 0.5*( XB(i) + XB(i+1) );
    Y(i,1)        = 0.5*( YB(i) + YB(i+1) );
    S(i,1)        = sqrt( (XB(i+1)-XB(i))^2 + (YB(i+1)-YB(i))^2 );
    theta(i,1)    = atan2( (YB(i+1)-YB(i)) , (XB(i+1)-XB(i)) );
    RHS(i,1)      = sin(theta(i) - alpha);
end

%% Calculation of Coefficients

for i = 1:M-1
    for j = 1:M-1
        if (i == j)
            CN1(i,j) = -1;
            CN2(i,j) = 1 ;
            CT1(i,j) = 0.5*pi;
            CT2(i,j) = 0.5*pi;
        else
            A = - (X(i) - XB(j))*(cos(theta(j))) - (Y(i) - YB(j))*(sin(theta(j)));
            B = (X(i) - XB(j))^2 + (Y(i) - YB(j))^2;
            C = sin(theta(i) - theta(j));
            D = cos(theta(i) - theta(j));
            E = (X(i) - XB(j))*sin(theta(j)) - (Y(i) - YB(j))*cos(theta(j));
            F = log(1 + ((S(j))^2 + (2*A*S(j))) / B);
            G = atan2((E*S(j)) , (B + A*S(j)));
            P = ((X(i) - XB(j)) * sin(theta(i) - 2*theta(j))) + ((Y(i) - YB(j)) * cos(theta(i) - 2*theta(j)));
            Q = ((X(i) - XB(j)) * cos(theta(i) - 2*theta(j))) - ((Y(i) - YB(j)) * sin(theta(i) - 2*theta(j)));
            
            CN2(i,j) = D + ((0.5*Q*F)/S(j)) - ((A*C + D*E)*(G/S(j)));
            CN1(i,j) = 0.5*D*F + C*G - CN2(i,j);
            CT2(i,j) = C + ((0.5*P*F)/S(j)) + ((A*D - C*E)*(G/S(j)));
            CT1(i,j) = 0.5*C*F - D*G - CT2(i,j);
        end
    end
end

%% Computation of Influence Coefficients
for i = 1:M-1
    AN(i,1) = CN1(i,1);
    AN(i,M) = CN2(1,M-1);
    AT(i,1) = CT1(i,1);
    AT(i,M) = CT2(i,M-1);
    for j = 2:M-1
        AN(i,j) = CN1(i,j) + CN2(i,j-1);
        AT(i,j) = CT1(i,j) + CT2(i,j-1);
    end
end
AN(M,1) = 1;
AN(M,M) = 1;
for j = 2:M-1
    AN(M,j) = 0;
end
RHS(M) = 0;

%% Solve for Gamma and velocity/pressure
Gama = AN\RHS;                                  % Solving for a syetem of linear equations

for i = 1:M-1
    V(i) = cos(theta(i)-alpha);
    for j = 1:M
        V(i) = V(i) + AT(i,j)*Gama(j);
        CP(i) = 1 - (V(i))^2;
    end
end

figure(1)
plot(X,CP);
set(gca,'Ydir','reverse');
xlabel('x/c');
ylabel('Coefficient of Pressure');
grid on;
grid minor;

%% Calculation of Lift Coefficient

CPl = CP(1:((M-1)/2));
CPl = flip(CPl);
CPu = CP((((M-1)/2)+1):end);

dCP = CPl - CPu;
dx = X((((M-1)/2)+1):end);

Cl(z) = trapz(dx,dCP);
z = z+1;