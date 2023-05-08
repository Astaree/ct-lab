%           circuit 2 to solve nodal
%   -----R41/1000/ ----B----- R42/4700/------
%   |                  |                     |
%   |                  |                     |
%   |              C41/10nF/                 |
%   |                  |                     |
%   A----C42/100nF/----C------L41/10mH/------D
%   |                                        |
%  R/50/                                     |
%   U//                                  C43/22nF/
%   |                                        |
%  GND                                      GND

%looking for voltages XYZ

%Frequencies 4kHz, 6kHz, 8kHz
freq = [4,6,8];
Zr41=[];
Zr42=[];
Zc41=[];
Zc42=[];
Zc43=[];
Zl41=[];

detG = [];
detX = [];
detY = [];
detZ = [];

Vx = [];
Vy = [];
Vz = [];

%ans
Il41 = [];

%changing capp and inductance into complex domain

%for 4kHz
U = [0.392 0.391 0.389];
flag=0; %there is resistor b4 voltage sorce with 50 Ohm, flag changes if we use it or not

for i = 1:3
   omega = [2 * pi() * freq(1)*1000 ,2 * pi() * freq(2)*1000,2 * pi() * freq(3)*1000];
Zr41 = [complex(1000)];
Zr42 = [complex(4700)];

Zc41 = [Zc41,Zc(omega(i),10*10^-9)];

%flat out +50Ohms from diagram
Zc42 = [Zc42,Zc(omega(i),100*10^-9)];
Zc43 = [Zc43,Zc(omega(i),22*10^-9)];

Zl41 = [Zl41,Zl(omega(i),0.01)];

%at A: Va(flag*(1/50)+1/Zc42(i)+1/Zr41)-Vb(1/Zr41 + 1/Zc41(i) +
%1/Zr42)-Vc(1/Zc42(i)+1/Zc41(i)+Zl41(i))-0Vd=U(i)/R50
%at B: -Va(1/Zr41)+Vb(1/Zr41+1/Zc41(i)+1/Zr42)-Vc(1/Zc41)-Vd(Zr42)=0
%at C:-Va(1/Zc42(i))-Vb(1/Zc41(i))+Vc(1/Zc42(i)+1/Zc41(i)+1/Zl41(i))-Vd(1/Zl41(i))=0
%at D:-0Va - Vb(1/Zr42)-Vc(1/Zl41(i))+Vd(1/Zc43(i)+1/Zl41(i)+1/Zr42)=0




G_matrix=[(flag*(1/50)+1/Zc42(i)+1/Zr41) -(1/Zr41)  -(1/Zc42(i)) 0;


    ];

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
