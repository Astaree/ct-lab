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
%   |                U/0.386V               |
%   |                    |                  |
%   |                   R50                 |
%   |                    |                  |
%  GND                  GND                GND

% eq A: U_a=U_src, skiped as it provides nothing
% eq X: U_x(1/R21+1/XC21+1/R22) - U_y(1/XC21) - U_z(1/R22) = 0
% eq Y: - U_x(1/XC21) + U_y(1/XC21 + 1/(XC22 + R50) + 1/XL21) - U_z(1/XL21) =
% = U/(XC22 + R50)
% eq Z: - U_x(1/R22) - U_y(1/XL21) + U_z(1/R22 + 1/XL21 + 1/XC23) = 0

% frequencies 4k, 6k, 8k
% voltages = .386, .340, .345 V

% we can form matrix equation G*U=A

% preperation:

matrix_G = complex(zeros(3,3,3));
matrix_I = complex(zeros(3,1,3));
matrix_V = complex(zeros(3,1,3));

detG = [0 0 0];
detX = [0 0 0];
detY = [0 0 0];
detZ = [0 0 0];

U_x = [0 0 0];
U_y = [0 0 0];
U_z = [0 0 0];


X_r_50 = [50 50 50];
X_r_22 = [1000 1000 1000];
X_r_21 = [1000 1000 1000];

X_c_21 = [47*10^-9 47*10^-9 47*10^-9];
X_c_22 = [10^-7 10^-7 10^-7];
X_c_23 = [10^-7 10^-7 10^-7];

X_l_21 = [0.01 0.01 0.01];
I_l_21 = [0 0 0];

freq = [4000 6000 8000];
U_src = [.386 .340 .345];

% calc omega from frequencies
omega = [calc_omega(freq(1)) calc_omega(freq(2)) calc_omega(freq(3))];

% calc complex resistance from incuctance and cappacitance
for i = 1:3
    X_c_21 (:,i) = Zc(omega(:,i),X_c_21(:,i));
    X_c_22 (:,i) = Zc(omega(:,i),X_c_22(:,i));
    X_c_23 (:,i) = Zc(omega(:,i),X_c_23(:,i));
    X_l_21 (:,i) = Zl(omega(:,i),X_l_21(:,i));
end

% eq A: U_a=U_src, skiped as it provides nothing
% eq X: U_x(1/R21+1/XC21+1/R22) - U_y(1/XC21) - U_z(1/R22) = 0
% eq Y: - U_x(1/XC21) + U_y(1/XC21 + 1/(XC22 + R50) + 1/XL21) - U_z(1/XL21) =
% = U/(XC22 + R50)
% eq Z: - U_x(1/R22) - U_y(1/XL21) + U_z(1/R22 + 1/XL21 + 1/XC23) = 0

% constructing matrix G that is 3x3 because A contributes "nothing":

for i = 1:3
matrix_G(:,:,i) = [
    ((1/X_r_21(i)) + (1/X_c_21(i)) + (1/X_r_22(i))) -(1/X_c_21(i)) -(1/X_r_22(i));
    -(1/X_c_21(i)) ((1/X_c_21(i)) + (1/(X_c_22(i) + X_r_50(i))) + (1/X_l_21(i))) -(1/X_l_21(i));
    -(1/X_r_22(i)) -(1/X_l_21(i)) ((1/X_r_22(i))+(1/X_l_21(i))+ (1/X_c_23(i)))
    ];
end



for i=1:3
matrix_I(:,:,i) = [
    0;
    (U_src(i)/(X_c_22(i)+X_r_50(i)));
    0];
end

%prep for det's, remember we only move in frequency domain!!!!
matr_X = matrix_G;
matr_X(:,1,:) = matrix_I(:,1,:);
matr_Y = matrix_G;
matr_Y(:,2,:) = matrix_I(:,1,:);
matr_Z = matrix_G;
matr_Z(:,3,:) = matrix_I(:,1,:);

for i = 1:3
detG(i) = det(matrix_G(:,:,i));
detX(i) = det(matr_X(:,:,i));
detY(i) = det(matr_Y(:,:,i));
detZ(i) = det(matr_Z(:,:,i));
end

%calculate voltages in nodes:

for i =1:3
U_x(i) = detX(i)/detG(i);
U_y(i) = detY(i)/detG(i);
U_z(i) = detZ(i)/detG(i);
end

%calculating voltage drop accros L21
rmsX = [abs(U_x(1)) rad2deg(angle(U_x(1)));
    abs(U_x(2)) rad2deg(angle(U_x(2)));
    abs(U_x(3)) rad2deg(angle(U_x(3)))
];
rmsY = [abs(U_y(1)) rad2deg(angle(U_y(1)));
    abs(U_y(2)) rad2deg(angle(U_y(2)));
    abs(U_y(3)) rad2deg(angle(U_y(3)))
];
rmsZ = [abs(U_z(1)) rad2deg(angle(U_z(1)));
    abs(U_z(2)) rad2deg(angle(U_z(2)));
    abs(U_z(3)) rad2deg(angle(U_z(3)))
];

%expected values
for i =1:3
I_l_21(i) = (U_z(i)-U_y(i))/(X_l_21(i));
end

%measured values:
U_m_Y = [.127*exp(deg2rad(24)*j);
    .290*exp(deg2rad(85)*j);
    .425*exp(deg2rad(29.5)*j)    ];

U_m_Z = [.271*exp( deg2rad(-.5)*j);
    .465*exp( deg2rad(-43.5)*j);
    .228*exp( deg2rad(-111)*j)    ];

expected_l_21 = (U_m_Z(:,:) - U_m_Y(:,:));
expected_l_21 = expected_l_21(:)/X_l_21(:);



