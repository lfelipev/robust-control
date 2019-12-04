
function [Er, a_max] = d_swing_up(eta_m, eta_g, Kg, kt, Rm, Mr, Lr)
    % Reference Energy (J)
    Er = 0.42;
    % Maximum torque for 5 V
    tau_max = eta_m*eta_g*Kg*kt*5/Rm;
    % Maximum acceleration of pivot (m/s^2)
    a_max = tau_max / (Mr*Lr);
end
