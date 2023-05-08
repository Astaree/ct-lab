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

% looking for voltages A, B, C, D

% eq A: U_a(1/R50 + 1/Xc42 + 1/R41) - U_b(1/r41) - U_c(1/Xc42) -U_d(0) =
% = U_src/R50
% eq B: -U_a(1/R41) + U_b (1/R41 + 1/Xc41 + 1/R42) - U_c (1/Xc41) -
% U_d(1/R42) = 0
% eq C: -U_a(1/Xc42) - U_b (1/Xc41) + U_c (1/Xc42 + 1/Xc41 +1/ Xl41) -
% U_d(1/Xl41) = 0
% eq D: -U_a(0) - U_b (1/R42) - U_c (1/Xl41) + U_d(1/R42 +1/Xl41 + 1/Xc43) = 0

%Frequencies 4kHz, 6kHz, 8kHz
% voltages = .392, .391, .389 V

% we can form matrix equation G*U=A

% preperation:

matrix_G = complex(zeros(4,4,3));
matrix_I = complex(zeros(4,1,3));
matrix_V = complex(zeros(4,1,3));

detG = [0 0 0];
detA = [0 0 0];
detB = [0 0 0];
detC = [0 0 0];
detD = [0 0 0];

U_a = [0 0 0];
U_b = [0 0 0];
U_c = [0 0 0];
U_d = [0 0 0];

X_r_50 = [50 50 50 50];
X_r_41 = [1000 1000 1000];
X_r_42 = [4700 4700 4700];

X_c_41 = [10*10^-9 10*10^-9 10*10^-9];
X_c_42 = [10^-7 10^-7 10^-7];
X_c_43 = [22*10^-9 22*10^-9 22*10^-9];

X_l_41 = [0.01 0.01 0.01];
I_C_41 = [0 0 0];

freq = [4000 6000 8000];
U_src = [.392 .391 .389];
omega = [calc_omega(freq(1)) calc_omega(freq(2)) calc_omega(freq(3))];

for i = 1:3
    X_c_41 (:,i) = Zc(omega(:,i),X_c_41(:,i));
    X_c_42 (:,i) = Zc(omega(:,i),X_c_42(:,i));
    X_c_43 (:,i) = Zc(omega(:,i),X_c_43(:,i));
    X_l_41 (:,i) = Zl(omega(:,i),X_l_41(:,i));
end

% eq A: U_a(1/R50 + 1/Xc42 + 1/R41) - U_b(1/r41) - U_c(1/Xc42) -U_d(0) =
% = U_src/R50
% eq B: -U_a(1/R41) + U_b (1/R41 + 1/Xc41 + 1/R42) - U_c (1/Xc41) -
% U_d(1/R42) = 0
% eq C: -U_a(1/Xc42) - U_b (1/Xc41) + U_c (1/Xc42 + 1/Xc41 +1/ Xl41) -
% U_d(1/Xl41) = 0
% eq D: -U_a(0) - U_b (1/R42) - U_c (1/Xl41) + U_d(1/R42 +1/Xl41 + 1/Xc43) = 0

for i = 1:3
matrix_G(:,:,i) = [
    (1/X_r_50(i) + 1/X_c_42(i) + 1/X_r_41(i)) -(1/X_r_41(i)) -(1/X_c_42(i)) -(0);
    -(1/X_r_41(i)) (1/X_r_41(i) + 1/X_c_41(i) + 1/X_r_42(i)) -(1/X_c_41(i)) -(1/X_r_42(i));
    -(1/X_c_42(i)) -(1/X_c_41(i)) (1/X_c_42(i) + 1/X_c_41(i) +1/X_l_41(i)) -(1/X_l_41(i));
    -(0) -(1/X_r_42(i)) -(1/X_l_41(i)) (1/X_r_42(i) +1/X_l_41(i) + 1/X_c_43(i))    ];
end

for i=1:3
matrix_I(:,:,i) = [
    U_src(i)/X_r_50(i);
    0;
    0;
    0];
end

matr_A = matrix_G;
matr_A(:,1,:) = matrix_I(:,1,:);
matr_B = matrix_G;
matr_B(:,2,:) = matrix_I(:,1,:);
matr_C = matrix_G;
matr_C(:,3,:) = matrix_I(:,1,:);
matr_D = matrix_G;
matr_D(:,4,:) = matrix_I(:,1,:);

for i = 1:3
detG(i) = det(matrix_G(:,:,i));
detA(i) = det(matr_A(:,:,i));
detB(i) = det(matr_B(:,:,i));
detC(i) = det(matr_C(:,:,i));
detD(i) = det(matr_D(:,:,i));
end

for i =1:3
U_a(i) = detA(i)/detG(i);
U_b(i) = detB(i)/detG(i);
U_c(i) = detC(i)/detG(i);
U_d(i) = detD(i)/detG(i);
end

%calculating voltage drop accros C41
rmsA = [abs(U_a(1)) rad2deg(angle(U_a(1)));
    abs(U_a(2)) rad2deg(angle(U_a(2)));
    abs(U_a(3)) rad2deg(angle(U_a(3)));
];
rmsB = [abs(U_b(1)) rad2deg(angle(U_b(1)));
    abs(U_b(2)) rad2deg(angle(U_b(2)));
    abs(U_b(3)) rad2deg(angle(U_b(3)));
];
rmsC = [abs(U_c(1)) rad2deg(angle(U_c(1)));
    abs(U_c(2)) rad2deg(angle(U_c(2)));
    abs(U_c(3)) rad2deg(angle(U_c(3)));
];
rmsD = [abs(U_d(1)) rad2deg(angle(U_d(1)));
    abs(U_d(2)) rad2deg(angle(U_d(2)));
    abs(U_d(3)) rad2deg(angle(U_d(3)));
];

%expected values
for i =1:3
I_C_41(i) = (U_c(i)-U_b(i))/(X_c_41(i));
end

%measured values:
U_m_B = [.384*exp(deg2rad(-2.3)*j);
    .390*exp(deg2rad(-4)*j);
    .396*exp(deg2rad(-8)*j)    ];

U_m_C = [.315*exp( deg2rad(-.5)*j);
    .302*exp( deg2rad(.5)*j);
    .273*exp( deg2rad(3)*j)    ];

expected_c_41 = (U_m_C(:,:) - U_m_B(:,:));
expected_c_41 = expected_c_41(:)/X_c_41(:);
