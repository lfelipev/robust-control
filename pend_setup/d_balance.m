function K = d_balance(A,B,zeta,wn,p3,p4)
    % Location of dominant poles along real-axis
    sigma = zeta*wn;
    % Location of dominant poles along img axis (damped natural freqency)
    wd = wn*sqrt(1-zeta^2);
    % Desired poles (-30 and -40 are given)
    DP = [-sigma+j*wd, -sigma-j*wd, p3, p4];
    % Find control gain using Matlab pole-placement command
    K = acker(A,B,DP);
end
