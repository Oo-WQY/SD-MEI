function M_data = SDMEI_BLOCH_SSFP_Mt(param)

gamma = param.gamma;

% param of pluse 
tp = param.pulseparam.tp;  % [ms] Single cycle excitation time
trep = param.pulseparam.trep;  % [ms] Periodic repetition interval
td = trep - tp;

fs = param.pulseparam.fs; 

% param of H
B1 = param.B1;
Delta = param.Delta;
theta = param.theta;
T1 = param.T1;
T2 = param.T2;
R1 = 1/T1;
R2 = 1/T2;

% Initial State M
M0 = 1;
Mx_0 = 0;
My_0 = 0;
Mz_0 = M0;
M = [Mx_0; My_0; Mz_0];


omega_1 = B1*gamma;

% Find the coefficients of the Bloch equation
% Excitation Process
% Spectral Diagonalization
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

% Calculate the matrix A at the sampling time point
ts_p = 1/fs : 1/fs : tp;
% If tp is not in ts_p, add
if abs(ts_p(end) - tp) > eps(tp)
    t_p = [ts_p, tp];
else
    t_p = ts_p;
end

expD11=exp(-Ds(1,1)*t_p);
expD22=exp(-Ds(2,2)*t_p);
expD33=exp(-Ds(3,3)*t_p);
expD = zeros(3, 3, length(t_p));
for k = 1:length(t_p)
    expD(:,:,k) = diag([expD11(k), expD22(k), expD33(k)]);
end

A_p = zeros(3, 3, length(t_p));
for k = 1:length(t_p)
    A_p(:,:,k) = Vs * expD(:,:,k) * inv(Vs); 
end
A11_p = real(squeeze(A_p(1,1,:))');
A22_p = real(squeeze(A_p(2,2,:))');
A12_p = real(squeeze(A_p(1,2,:))');
A21_p = real(squeeze(A_p(2,1,:))');
A23_p = real(squeeze(A_p(2,3,:))');
A32_p = real(squeeze(A_p(3,2,:))');
A13_p = real(squeeze(A_p(1,3,:))');
A31_p = real(squeeze(A_p(3,1,:))');
A33_p = real(squeeze(A_p(3,3,:))'); 

% The magnetization vector in steady state
Mx_ss_p = R1*Delta*omega_1*M0/(R1*R2^2+R1*Delta^2+R2*omega_1^2);
My_ss_p = R1*R2*omega_1*M0/(R1*R2^2+R1*Delta^2+R2*omega_1^2);
Mz_ss_p = R1*(R2^2+Delta^2)*M0/(R1*R2^2+R1*Delta^2+R2*omega_1^2);


% Relaxation Process
% Spectral Diagonalization 
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

% Relaxation Process equation coefficients
% Calculate sampling time points
ts_r = 1/fs : 1/fs : td;
% If tp is not in ts_p, add
if abs(ts_r(end) - td) > eps(td)
    t_r = [ts_r, td];
else
    t_r = ts_r;
end

% Calculate the matrix A at the sampling time point
expD11_r=exp(-Ds_r(1,1)*t_r);
expD22_r=exp(-Ds_r(2,2)*t_r);
expD33_r=exp(-Ds_r(3,3)*t_r);
expD_r = zeros(3, 3, length(t_r));
for k = 1:length(t_r)
    expD_r(:,:,k) = diag([expD11_r(k), expD22_r(k), expD33_r(k)]);
end

A_r = zeros(3, 3, length(t_r));
for k = 1:length(t_r)
    A_r(:,:,k) = Vs_r * expD_r(:,:,k) * inv(Vs_r); 
end

A11_r = real(squeeze(A_r(1,1,:))');
A22_r = real(squeeze(A_r(2,2,:))');
A12_r = real(squeeze(A_r(1,2,:))');
A21_r = real(squeeze(A_r(2,1,:))');
A23_r = real(squeeze(A_r(2,3,:))');
A32_r = real(squeeze(A_r(3,2,:))');
A13_r = real(squeeze(A_r(1,3,:))');
A31_r = real(squeeze(A_r(3,1,:))');
A33_r = real(squeeze(A_r(3,3,:))');

Mx_ss_r = 0;
My_ss_r = 0;
Mz_ss_r = M0;


% Steady-state solution process
max_iter = 100;     % Set a maximum number of iterations to prevent infinite loops
tol = 1e-12;        % Absolute error threshold
i = 0;
abs_error = inf;

t = [0];
tadd = 0;

while abs_error > tol && i < max_iter  
    i = i+1;

    % Excitation process
    Mx = A11_p*(Mx_0 - Mx_ss_p) + A12_p*(My_0 - My_ss_p) + A13_p*(Mz_0 - Mz_ss_p) + Mx_ss_p;
    My = A21_p*(Mx_0 - Mx_ss_p) + A22_p*(My_0 - My_ss_p) + A23_p*(Mz_0 - Mz_ss_p) + My_ss_p;
    Mz = A31_p*(Mx_0 - Mx_ss_p) + A32_p*(My_0 - My_ss_p) + A33_p*(Mz_0 - Mz_ss_p) + Mz_ss_p;
           
    if i > 1
        % Load the magnetization vector after the last excitation
        Mx_prev = M_prev(1, end);
        My_prev = M_prev(2, end);   
        Mz_prev = M_prev(3, end);
        % Calculate the maximum absolute error of the magnetization vector
        abs_error = max([abs(Mx(end) - Mx_prev), abs(My(end) - My_prev), abs(Mz(end) - Mz_prev)]);
    end

    Mp = [Mx; My; Mz];
    M = [M, Mp];
    M_prev = Mp;

    t = [t, tadd+t_p];
    tadd = tadd + tp;
    
    % Relaxation Process
    Mx_0 = M(1, end);
    My_0 = M(2, end);
    Mz_0 = M(3, end);
    
    Mx = A11_r*(Mx_0 - Mx_ss_r) + A12_r*(My_0 - My_ss_r) + A13_r*(Mz_0 - Mz_ss_r) + Mx_ss_r;
    My = A21_r*(Mx_0 - Mx_ss_r) + A22_r*(My_0 - My_ss_r) + A23_r*(Mz_0 - Mz_ss_r) + My_ss_r;
    Mz = A31_r*(Mx_0 - Mx_ss_r) + A32_r*(My_0 - My_ss_r) + A33_r*(Mz_0 - Mz_ss_r) + Mz_ss_r;
    
    Mr = [Mx; My; Mz];
    M = [M, Mr];
    
    t = [t, tadd+t_r];
    tadd = tadd + td;

    Mx_0 = M(1, end);
    My_0 = M(2, end);
    Mz_0 = M(3, end);
end


M_data = [t; M];