function Xdot = propagateOrbit(X, t, P)
    % Set local variables
    Sun = 'SUN';
    Mars = '499';
    Deimos = '402';
    rf = 'ECLIPJ2000';
    abcorr = 'NONE';
    R_1to2 = cspice_pxform(rf, 'IAU_DEIMOS', t);
    R_2to1 = cspice_pxform('IAU_DEIMOS', rf, t);
   
    % Unpack state vector
    rVec = X(1:3);
    norm_rVec = norm(rVec);
    vVec = X(4:6);

    % Find Deimos 2 body acceleration term
    acc_D2 = -P.Mu_deimos*rVec/norm_rVec^3;

    % Find Deimos nonspherical acceleration term
    [acc_Dns, ~] = SHacc4x4only(R_1to2*rVec, P.R_deimos, P.Mu_deimos, P.J, P.C, P.S);
    acc_Dns = R_2to1*acc_Dns;

    % Find position of Mars wrt Deimos
    M = mice_spkezr(Mars, t, rf, abcorr, Deimos);
    M_pos = M.state(1:3);
    norm_M_pos = norm(M_pos);
    % Find position of Mars wrt Spacecraft
    M_pos2 = M_pos - rVec;
    norm_M_pos2 = norm(M_pos2);
    % Find Mars three body acceleration term
    acc_M3 = P.Mu_mars*(M_pos2/norm_M_pos2^3-M_pos/norm_M_pos^3);

    % Find position of Sun wrt to Deimos
    S = mice_spkezr(Sun, t, rf, abcorr, Deimos);
    S_pos = S.state(1:3);
    norm_S_pos = norm(S_pos);
    % Find position of Sun wrt Spacecraft
    S_pos2 = S_pos - rVec;
    norm_S_pos2 = norm(S_pos2);
    % Find Sun three body acceleration term
    acc_S3 = P.Mu_sun*(S_pos2/norm_S_pos2^3-S_pos/norm_S_pos^3);

    % Choose total acceleration
    if P.mode == 1
        acc = acc_D2 + acc_Dns + acc_M3 + acc_S3;
    elseif P.mode == 2
        acc = acc_D2;
    elseif P.mode == 3
        acc = acc_D2 + acc_Dns;
    elseif P.mode == 4
        acc = acc_D2 + acc_M3 + acc_S3;
    elseif P.mode == 5
        acc = acc_D2 + acc_M3;
    end

    % Calculate first derivatives with respect to time of each variable
    rVecDot = vVec;
    vVecDot = acc;

    % Pack Xdot vector
    Xdot = [rVecDot; vVecDot];
end