%           circuit 1 to solve nodal
%    ----------------R22/1000/--------------
%   |                                       |
%   |                                       |
%   |                                       |
%   X---C21/47*10^-9F/---Y----L21/0.01H/----Z
%   |                    |                  |
%   |                    |                  |
%   R21/1000/       C22/10^-7F/        C23/10^-7F/
%   |                    |                  |
%   |                    A                  |
%   |                    |                  |
%   |                   R50                 |
%   |                    |                  |
%   |                 U/0.386V              |
%   |                    |                  |
%  GND                  GND                GND

%looking for voltages XYZ

%Frequencies 4kHz, 6kHz, 8kHz
freq = [4,6,8];
Zr21=[];
Zr22=[];
Zc21=[];
Zc22=[];
Zc23=[];
Zl21=[];

detG = [];
detX = [];
detY = [];
detZ = [];

Vx = [];
Vy = [];
Vz = [];

%ans
Il21 = [];

%changing capp and inductance into complex domain

%for 4kHz
U = [0.386 0.34 0.345];
i=1; %ctr
flag=0; %there is resistor b4 voltage sorce with 50 Ohm, flag changes if we use it or not

for i = 1:3
   omega = [2 * pi() * freq(1)*1000 ,2 * pi() * freq(2)*1000,2 * pi() * freq(3)*1000];
Zr21 = [complex(1000)];
Zr22 = [complex(1000)];

Zc21 = [Zc21,Zc(omega(i),47*10^-9)];

%flat out +50Ohms from diagram
Zc22 = [Zc22,Zc(omega(i),10^-7)+flag*50];
Zc23 = [Zc23,Zc(omega(i),10^-7)];

Zl21 = [Zl21,Zl(omega(i),0.01)];

%at X: Vx(1/Zr21+1/Zc21+1/Zr22)-Vy(1/Zc21)-Vz(1/Zr22)=0
%at Y: -Vx(1/Zc21)+Vy(1/(Zc22+50) + 1/Zl21 +1/Zc21)-Vx(1/Zc21)-Vz(1/Zl21)= U/Zc22
%at Z: -Vx(1/Zr22)-Vy(1/Zl21)+Vz(1/Zl21+1/Zc23+1/Zr22)=0

G_matrix=[(1/Zr21+1/Zc21(i)+1/Zr22) (-1/Zc21(i)) (-1/Zr22);
    (-1/Zc21(i)) (1/(Zc22(i)+flag*50)+1/Zl21(i)+1/Zc21(i)) (-1/Zc21(i));
    (-1/Zr22) (-1/Zl21(i)) (1/Zl21(i)+1/Zc23(i)+1/Zr22)];

%V_matrix = [Vx;Vy;Vz];

I_matrix = [0;(U(i)/(Zc22(i)+flag*50));0];

detG = [detG, det(G_matrix)];

Vx_matrix = G_matrix';
Vx_matrix(:,1) = I_matrix;
detX = [detX,det(Vx_matrix)];

Vy_matrix = G_matrix';
Vy_matrix(:,2) = I_matrix;
detY = [detY, det(Vy_matrix)];

Vz_matrix = G_matrix';
Vz_matrix(:,3) = I_matrix;
detZ = [detZ, det(Vz_matrix)];

Vx = [Vx,detX(i)/detG(i)];
Vy = [Vy,detY(i)/detG(i)];
Vz = [Vz,detZ(i)/detG(i)];

Il21 = [Il21,((Vz(i)-Vy(i))/Zl21(i))];
end;
%changing to polar form so we can compare results

polarVx = [0 0; 0 0; 0 0];
polarVy = [0 0; 0 0; 0 0];
polarVz = [0 0; 0 0; 0 0];
for i = 1:3
    polarVx(i,:) = [abs(Vx(i)) angle(Vx(i))];
    polarVy(i,:) = [abs(Vy(i)) angle(Vy(i))];
    polarVz(i,:) = [abs(Vz(i)) angle(Vz(i))];
end
%Current across L21 for 
for i = 1:3
    fprintf("Current going through Incuctor L21 = %f%+fj, at frequency %dkHz\n",real(Il21(i)),imag(Il21(i)),freq(i))
end
