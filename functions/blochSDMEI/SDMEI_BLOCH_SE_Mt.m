function M_data = SDMEI_BLOCH_SE_Mt(param)

gamma = param.gamma;

% param of pluse 
fs = param.pulseparam.fs; 

% param of H
B1 = param.B1;
Delta_abs = param.Delta;
theta = param.theta;
T1 = param.T1;
T2 = param.T2;
T2_1 = param.T2_eff;
R1 = 1/T1;
R2 = 1/T2;

% Delta sampling
Delta_data = -Delta_abs : 1 : Delta_abs;

% Initialize M_all to save the magnetization vector M corresponding to each Delta
M_all = cell(1, length(Delta_data));

for Dn = 1:length(Delta_data)
    Delta = Delta_data(Dn);
    
    % Initial State M
    M0 = 1;
    Mx_0 = 0;
    My_0 = 0;
    Mz_0 = M0;
    M = [Mx_0; My_0; Mz_0];
    t = [0];
    tadd = 0;
    
    % Excitation relaxation time corresponding to SE sequence
    tp = param.pulseparam.tp;
    trep = param.pulseparam.trep;
    td = trep - tp;
    t_int = [tp, td, 2*tp, 2*td];

    % Parameter setting of H
    B1 = param.B1;
    B1_Data = [B1, 0, B1, 0];
    theta_Data = [0, 0, pi/2, pi/2];

    for tnum = 1:length(t_int)
        tp = t_int(tnum);
        B1 = B1_Data(tnum);
        theta = theta_Data(tnum);
        
        % Calculate sampling time points
        t_s = 1/fs : 1/fs : tp;
        % If tp is not in t_s, add
        if abs(t_s(end) - tp) > eps(tp)
            t_s = [t_s, tp];
        else
            t_s = t_s;
        end
        
        % save the value of t
        t = [t, tadd + t_s];
        tadd = tadd + tp;

        % Spectral Diagonalization
        omega_1 = B1*gamma;
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
        
        % Solve the coefficient matrix A of the Bloch equation
        expD11=exp(-Ds(1,1)*t_s);
        expD22=exp(-Ds(2,2)*t_s);
        expD33=exp(-Ds(3,3)*t_s);
        expD = zeros(3, 3, length(t_s));
        for k = 1:length(t_s)
            expD(:,:,k) = diag([expD11(k), expD22(k), expD33(k)]);
        end
        
        A_p = zeros(3, 3, length(t_s));
        for k = 1:length(t_s)
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
        
        % Bloch Equation
        % Solving for the magnetization vector
        Mx = A11_p*(Mx_0 - Mx_ss_p) + A12_p*(My_0 - My_ss_p) + A13_p*(Mz_0 - Mz_ss_p) + Mx_ss_p;
        My = A21_p*(Mx_0 - Mx_ss_p) + A22_p*(My_0 - My_ss_p) + A23_p*(Mz_0 - Mz_ss_p) + My_ss_p;
        Mz = A31_p*(Mx_0 - Mx_ss_p) + A32_p*(My_0 - My_ss_p) + A33_p*(Mz_0 - Mz_ss_p) + Mz_ss_p;
        
        Mp = [Mx;My;Mz];
        M = [M,Mp];
        
        Mx_0 = M(1, end);
        My_0 = M(2, end);
        Mz_0 = M(3, end);
    end
    M_all{Dn} = M;  
end

% Weighted solution
T2_ih = 1/T2_1-1/T2;
y=1./(Delta_data.^2+T2_ih^2);
W = y./sum(y);
ln = length(t);
Mx = zeros(1,ln);
My = zeros(1,ln);
Mz = zeros(1,ln);
for dn = 1:length(Delta_data)
    MData = M_all{dn};
    Mx = Mx + W(dn)*MData(1,:);
    My = My + W(dn)*MData(2,:);
    Mz = Mz + W(dn)*MData(3,:);
end
M_data = [t; Mx; My; Mz];