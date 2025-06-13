function M_data = SDMEI_BLOCH_SE_Mt(param)

gamma = param.gamma;

% param of pluse 
tp = param.pulseparam.tp;
trep = param.pulseparam.trep;
td = trep - tp;

% param of H
B1 = param.B1;
Delta_abs = param.Delta;
theta = param.theta;
T1 = param.T1;
T2 = param.T2;
R1 = 1/T1;
R2 = 1/T2;

% Delta
Delta_data = -Delta_abs : 1 : Delta_abs;

M_all = cell(1, length(Delta_data));

for Dn = 1:length(Delta_data)
Delta = Delta_data(Dn);

% Initial State M
M0 = 1;
Mx_0 = 0;
My_0 = 0;
Mz_0 = M0;
M = [Mx_0; My_0; Mz_0];


% Spectral Diagonalization of the Excitation Process
omega_1 = B1*gamma;
K = [R2, -Delta, omega_1*sin(theta);
     Delta, R2, -omega_1*cos(theta);
     -omega_1*sin(theta), omega_1*cos(theta), R1];
[V,D] = eig(K);
D1 = diag(D);
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

% Calculate sampling time points
ts_p = 1/fs : 1/fs : tp;
% If tp is not in ts_p, add
if abs(ts_p(end) - tp) > eps(tp)
    t_p = [ts_p, tp];
else
    t_p = ts_p;
end

% save the value of t
t = [0, t_p];

% Calculate the matrix A at the sampling time point
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

% Calculate the magnetization vector at the sampling time point
Mx = A11_p*(Mx_0 - Mx_ss_p) + A12_p*(My_0 - My_ss_p) + A13_p*(Mz_0 - Mz_ss_p) + Mx_ss_p;
My = A21_p*(Mx_0 - Mx_ss_p) + A22_p*(My_0 - My_ss_p) + A23_p*(Mz_0 - Mz_ss_p) + My_ss_p;
Mz = A31_p*(Mx_0 - Mx_ss_p) + A32_p*(My_0 - My_ss_p) + A33_p*(Mz_0 - Mz_ss_p) + Mz_ss_p;

Mp = [Mx;My;Mz];
M = [M,Mp];

Mx_0 = M(1, end);
My_0 = M(2, end);
Mz_0 = M(3, end);

% td
% Spectral Diagonalization of the Relaxation Process
omega_1 = 0;
K = [R2, -Delta, omega_1*sin(theta);
     Delta, R2, -omega_1*cos(theta);
     -omega_1*sin(theta), omega_1*cos(theta), R1];
[V,D] = eig(K);
D1 = diag(D);
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

% Calculate sampling time points
ts_r = 1/fs : 1/fs : td;
% If tp is not in ts_p, add
if abs(ts_r(end) - td) > eps(td)
    t_r = [ts_r, td];
else
    t_r = ts_r;
end

% save the value of t
t = [t, td+t_r];

% Calculate the matrix A at the sampling time point
expD11=exp(-Ds(1,1)*t_r);
expD22=exp(-Ds(2,2)*t_r);
expD33=exp(-Ds(3,3)*t_r);
expD = zeros(3, 3, length(t_r));
for k = 1:length(t_r)
    expD(:,:,k) = diag([expD11(k), expD22(k), expD33(k)]);
end

A_r = zeros(3, 3, length(t_r));
for k = 1:length(t_r)
    A_r(:,:,k) = Vs * expD(:,:,k) * inv(Vs); 
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

Mx_0 = M(1, end);
My_0 = M(2, end);
Mz_0 = M(3, end);

Mx = A11_r*(Mx_0 - Mx_ss_r) + A12_r*(My_0 - My_ss_r) + A13_r*(Mz_0 - Mz_ss_r) + Mx_ss_r;
My = A21_r*(Mx_0 - Mx_ss_r) + A22_r*(My_0 - My_ss_r) + A23_r*(Mz_0 - Mz_ss_r) + My_ss_r;
Mz = A31_r*(Mx_0 - Mx_ss_r) + A32_r*(My_0 - My_ss_r) + A33_r*(Mz_0 - Mz_ss_r) + Mz_ss_r;

Mr = [Mx; My; Mz];
M = [M, Mr];

% 90 tp
omega_1 = B1*gamma;
theta = pi/2;
K = [R2, -Delta, omega_1*sin(theta);
     Delta, R2, -omega_1*cos(theta);
     -omega_1*sin(theta), omega_1*cos(theta), R1];
[V,D] = eig(K);
D1 = diag(D);
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

% Calculate sampling time points
ts_p = 1/fs : 1/fs : 2*tp;
% If tp is not in ts_p, add
if abs(ts_p(end) - 2*tp) > eps(2*tp)
    t_p = [ts_p, 2*tp];
else
    t_p = ts_p;
end

% save the value of t
t = [t, trep+t_p];

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

Mx_ss_p = -R1*R2*omega_1*M0/(R1*R2^2+R1*Delta^2+R2*omega_1^2);
My_ss_p = R1*Delta*omega_1*M0/(R1*R2^2+R1*Delta^2+R2*omega_1^2);
Mz_ss_p = R1*(R2^2+Delta^2)*M0/(R1*R2^2+R1*Delta^2+R2*omega_1^2);

Mx_0 = M(1, end);
My_0 = M(2, end);
Mz_0 = M(3, end);

Mx = A11_p*(Mx_0 - Mx_ss_p) + A12_p*(My_0 - My_ss_p) + A13_p*(Mz_0 - Mz_ss_p) + Mx_ss_p;
My = A21_p*(Mx_0 - Mx_ss_p) + A22_p*(My_0 - My_ss_p) + A23_p*(Mz_0 - Mz_ss_p) + My_ss_p;
Mz = A31_p*(Mx_0 - Mx_ss_p) + A32_p*(My_0 - My_ss_p) + A33_p*(Mz_0 - Mz_ss_p) + Mz_ss_p;

Mp = [Mx;My;Mz];
M = [M,Mp];

% 2*td
omega_1 = 0;
K = [R2, -Delta, omega_1*sin(theta);
     Delta, R2, -omega_1*cos(theta);
     -omega_1*sin(theta), omega_1*cos(theta), R1];
[V,D] = eig(K);
D1 = diag(D);
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

% Calculate sampling time points
ts_r = 1/fs : 1/fs : td;
% If tp is not in ts_p, add
if abs(ts_r(end) - td) > eps(td)
    t_r = [ts_r, td];
else
    t_r = ts_r;
end

% save the value of t
t = [t, td+t_r];

t_r = 1/fs : 1/fs : 2*td;
expD11=exp(-Ds(1,1)*t_r);
expD22=exp(-Ds(2,2)*t_r);
expD33=exp(-Ds(3,3)*t_r);
expD = zeros(3, 3, length(t_r));
for k = 1:length(t_r)
    expD(:,:,k) = diag([expD11(k), expD22(k), expD33(k)]);
end

A_r = zeros(3, 3, length(t_r));
for k = 1:length(t_r)
    A_r(:,:,k) = Vs * expD(:,:,k) * inv(Vs); 
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

Mx_0 = M(1, end);
My_0 = M(2, end);
Mz_0 = M(3, end);

Mx = A11_r*(Mx_0 - Mx_ss_r) + A12_r*(My_0 - My_ss_r) + A13_r*(Mz_0 - Mz_ss_r) + Mx_ss_r;
My = A21_r*(Mx_0 - Mx_ss_r) + A22_r*(My_0 - My_ss_r) + A23_r*(Mz_0 - Mz_ss_r) + My_ss_r;
Mz = A31_r*(Mx_0 - Mx_ss_r) + A32_r*(My_0 - My_ss_r) + A33_r*(Mz_0 - Mz_ss_r) + Mz_ss_r;

Mr = [Mx; My; Mz];
M = [M, Mr];

t = 0 : 1/fs : (length(M)-1)*1e-2;
M_data = [t;M];

M_data = [B1_data; M_xss_B1; M_yss_B1; M_zss_B1];