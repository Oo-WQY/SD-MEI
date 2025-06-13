function M_data = SDMEI_BLOCH_SSFP(param)

gamma = param.gamma;

% param of pluse 
tp = param.pulseparam.tp;  % [ms] Single cycle excitation time
trep = param.pulseparam.trep;  % [ms] Periodic repetition interval

% param of H
B1_data = param.B1;
Delta = param.Delta;
theta = param.theta;
T1 = param.T1;
T2 = param.T2;
R1 = 1/T1;
R2 = 1/T2;

% Initial State M
M0 = 1;

% Steady state M
M_xss_B1 = [];
M_yss_B1 = [];
M_zss_B1 = [];


for bn = 1:length(B1_data)
Mx_0 = 0;
My_0 = 0;
Mz_0 = M0;
M = [Mx_0; My_0; Mz_0];

B1 = B1_data(bn);
omega_1 = B1*gamma;

% Find the coefficients of the Bloch equation
% Spectral Diagonalization of the Excitation Process
H = [R2, -Delta, omega_1*sin(theta);
     Delta, R2, -omega_1*cos(theta);
     -omega_1*sin(theta), omega_1*cos(theta), R1];
[V,D] = eig(H);
D1 = diag(D);

% Adjust the eigenvalues and eigenvectors to make them 
% consistent with the spin transient dynamics process
if numel(find(sign(imag(D1))==0))>1
    [~, ind] = sort(D1);
    indc = ind(3);
    ind(3) = ind(1);
    ind(1) = indc;
else
    ind(1) = find(sign(imag(D1))==0);
    ind(2) = find(sign(imag(D1))==1);
    ind(3) = find(sign(imag(D1))==-1);
end
Ds = D(ind,ind);
Vs = V(:,ind);

% Excitation process equation coefficients
expD11=exp(-Ds(1,1)*tp);
expD22=exp(-Ds(2,2)*tp);
expD33=exp(-Ds(3,3)*tp);
expD = diag([expD11, expD22, expD33]);

A_p = Vs * expD * inv(Vs); 

A11_p = real(A_p(1,1));
A22_p = real(A_p(2,2));
A12_p = real(A_p(1,2));
A21_p = real(A_p(2,1));
A23_p = real(A_p(2,3));
A32_p = real(A_p(3,2));
A13_p = real(A_p(1,3));
A31_p = real(A_p(3,1));
A33_p = real(A_p(3,3)); 

Mx_ss_p = R1*Delta*omega_1*M0/(R1*R2^2+R1*Delta^2+R2*omega_1^2);
My_ss_p = R1*R2*omega_1*M0/(R1*R2^2+R1*Delta^2+R2*omega_1^2);
Mz_ss_p = R1*(R2^2+Delta^2)*M0/(R1*R2^2+R1*Delta^2+R2*omega_1^2);


% Spectral Diagonalization of the Relaxation Process
H_r = [R2, -Delta, 0;
     Delta, R2, 0;
     0, 0, R1];
[V_r,D_r] = eig(H_r);
D1_r = diag(D_r);

% Adjust the eigenvalues and eigenvectors to make them 
% consistent with the spin transient dynamics process
if numel(find(sign(imag(D1_r))==0))>1
    [~, ind] = sort(D1_r);
    indc = ind(3);
    ind(3) = ind(1);
    ind(1) = indc;
else
    ind(1) = find(sign(imag(D1_r))==0);
    ind(2) = find(sign(imag(D1_r))==1);
    ind(3) = find(sign(imag(D1_r))==-1);
end
Ds_r = D_r(ind,ind);
Vs_r = V_r(:,ind);

% Excitation process equation coefficients
tr = trep-tp;
expD11_r=exp(-Ds_r(1,1)*tr);
expD22_r=exp(-Ds_r(2,2)*tr);
expD33_r=exp(-Ds_r(3,3)*tr);
expD_r = diag([expD11_r, expD22_r, expD33_r]);

A_r = Vs_r * expD_r * inv(Vs_r); 

A11_r = real(A_r(1,1));
A22_r = real(A_r(2,2));
A12_r = real(A_r(1,2));
A21_r = real(A_r(2,1));
A23_r = real(A_r(2,3));
A32_r = real(A_r(3,2));
A13_r = real(A_r(1,3));
A31_r = real(A_r(3,1));
A33_r = real(A_r(3,3)); 

Mx_ss_r = 0;
My_ss_r = 0;
Mz_ss_r = M0;


% Steady-state solution process
max_iter = 100;     % Set a maximum number of iterations to prevent infinite loops
tol = 1e-12;        % Absolute error threshold
i = 0;
abs_error = inf;

while abs_error > tol && i < max_iter
    i = i+1;

    Mx = A11_p*(Mx_0 - Mx_ss_p) + A12_p*(My_0 - My_ss_p) + A13_p*(Mz_0 - Mz_ss_p) + Mx_ss_p;
    My = A21_p*(Mx_0 - Mx_ss_p) + A22_p*(My_0 - My_ss_p) + A23_p*(Mz_0 - Mz_ss_p) + My_ss_p;
    Mz = A31_p*(Mx_0 - Mx_ss_p) + A32_p*(My_0 - My_ss_p) + A33_p*(Mz_0 - Mz_ss_p) + Mz_ss_p;
    
    Mp = [Mx; My; Mz];
    M = [M, Mp];
    if i > 1
        % Load the magnetization vector after the last excitation
        Mx_prev = M(1, end-2);
        My_prev = M(2, end-2);   
        Mz_prev = M(3, end-2);
        % Calculate the maximum absolute error of the magnetization vector
        abs_error = max([abs(Mx - Mx_prev), abs(My - My_prev), abs(Mz - Mz_prev)]);
    end
    
    Mx_0 = M(1, end);
    My_0 = M(2, end);
    Mz_0 = M(3, end);
    
    Mx = A11_r*(Mx_0 - Mx_ss_r) + A12_r*(My_0 - My_ss_r) + A13_r*(Mz_0 - Mz_ss_r) + Mx_ss_r;
    My = A21_r*(Mx_0 - Mx_ss_r) + A22_r*(My_0 - My_ss_r) + A23_r*(Mz_0 - Mz_ss_r) + My_ss_r;
    Mz = A31_r*(Mx_0 - Mx_ss_r) + A32_r*(My_0 - My_ss_r) + A33_r*(Mz_0 - Mz_ss_r) + Mz_ss_r;
    
    Mr = [Mx; My; Mz];
    M = [M, Mr];
    
    Mx_0 = M(1, end);
    My_0 = M(2, end);
    Mz_0 = M(3, end);
end

% Save the magnetization vector in steady state
M_xss_B1 = [M_xss_B1, Mp(1)];
M_yss_B1 = [M_yss_B1, Mp(2)];
M_zss_B1 = [M_zss_B1, Mp(3)];
end

M_data = [B1_data; M_xss_B1; M_yss_B1; M_zss_B1];