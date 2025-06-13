function [T_single,N_now,Flag_on] = getsingletime(t,Trep,Tp)
N_rep = floor(t/Trep);
N_now = N_rep + 1;
T_single = t-(N_rep) * Trep;
if T_single > Tp
    Flag_on = 0;
else
    Flag_on = 1;
end

