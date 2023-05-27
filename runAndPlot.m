function [] = runAndPlot(TOF, OE_initial, P)
    Sun = 'SUN';
    Mars = '499';
    Deimos = '402';
    frame1 = 'ECLIPJ2000';
    frame2 = 'IAU_DEIMOS';
    abcorr = 'NONE';
    R_1to3 = cspice_pxform(frame1, frame2, P.ET_initial);
    R_3to1 = cspice_pxform(frame2, frame1, P.ET_initial);
    
    X_initial = oe2vec(OE_initial, P.Mu_deimos);
    X_initial(1:3) = R_3to1*X_initial(1:3);
    X_initial(4:6) = R_3to1*X_initial(4:6);
    
    tspan = linspace(P.ET_initial,(P.ET_initial+TOF*24*60*60),TOF*200);
    eps = 1e-10;
    options = odeset('RelTol',eps,'AbsTol',eps,'Event',@meanRadiusFct);
    
    Xdot = @(t, X)propagateOrbit(X, t, P);
    [tVec, X, ~, ~, ~] = ode45(Xdot, tspan, X_initial, options);
    X_length = length(X);
    tDays = (tVec - tVec(1))/(24*60*60);
    
    for i = 1:X_length
        ET = P.ET_initial + i*24*60*60/200;
    
        R_1to2 = cspice_pxform(frame1, frame2, ET);
        R_2to1 = cspice_pxform(frame2, frame1, ET);
    
        OE(i,:) = vec2oe(P.Mu_deimos, R_1to2*X(i,1:3)', R_1to2*X(i,4:6)');
    
        rVec = X(i,1:3);
        norm_rVec = norm(rVec);
        r(i) = abs(norm_rVec - P.R_deimos);
        X_f2(i,1:3) = R_1to2*X(i,1:3)';
        X_f3(i,1:3) = R_1to3*X(i,1:3)';
    
        % Find Deimos 2 body acceleration term
        acc_D2(i) = norm(-P.Mu_deimos*rVec/norm_rVec^3);
    
        % Find Deimos nonspherical acceleration term
        [acc_Dns_a, ~] = SHacc4x4only(R_1to2*rVec', P.R_deimos, P.Mu_deimos, P.J, P.C, P.S);
        acc_Dns(i) = norm(R_2to1*acc_Dns_a);
    
        % Find position of Mars wrt Deimos
        M = mice_spkezr(Mars, ET, frame1, abcorr, Deimos);
        M_pos = M.state(1:3);
        norm_M_pos = norm(M_pos);
        % Find position of Mars wrt Spacecraft
        M_pos2 = M_pos - rVec';
        norm_M_pos2 = norm(M_pos2);
        % Find Mars three body acceleration term
        acc_M3(i) = norm(P.Mu_mars*(M_pos2/norm_M_pos2^3-M_pos/norm_M_pos^3));
    
        % Find position of Sun wrt to Deimos
        S = mice_spkezr(Sun, ET, frame1, abcorr, Deimos);
        S_pos = S.state(1:3);
        norm_S_pos = norm(S_pos);
        % Find position of Sun wrt Spacecraft
        S_pos2 = S_pos - rVec';
        norm_S_pos2 = norm(S_pos2);
        % Find Sun three body acceleration term
        acc_S3(i) = norm(P.Mu_sun*(S_pos2/norm_S_pos2^3-S_pos/norm_S_pos^3));
    end
    
    frame1 = figure;
    tiledlayout(2,2);
    nexttile
    plot(X(1,1),X(1,2),'>','MarkerEdgeColor','green','MarkerFaceColor','green')
    hold on
    plot(X(end,1),X(end,2),'x','MarkerEdgeColor','red','MarkerFaceColor','red')
    ellipsoid(0,0,0,P.R_deimos,P.R_deimos,P.R_deimos)
    plot(X(:,1),X(:,2))
    xlabel('x (km)')
    ylabel('y (km)')
    title('Frame 1: ECLIPJ2000 frame')
    axis equal
    grid on
    nexttile
    plot(X(1,2),X(1,3),'>','MarkerEdgeColor','green','MarkerFaceColor','green')
    hold on
    plot(X(end,2),X(end,3),'x','MarkerEdgeColor','red','MarkerFaceColor','red')
    ellipsoid(0,0,0,P.R_deimos,P.R_deimos,P.R_deimos)
    plot(X(:,2),X(:,3))
    xlabel('y (km)')
    ylabel('z (km)')
    axis equal
    grid on
    nexttile
    plot(X(1,1),X(1,3),'>','MarkerEdgeColor','green','MarkerFaceColor','green')
    hold on
    plot(X(end,1),X(end,3),'x','MarkerEdgeColor','red','MarkerFaceColor','red')
    ellipsoid(0,0,0,P.R_deimos,P.R_deimos,P.R_deimos)
    plot(X(:,1),X(:,3))
    xlabel('x (km)')
    ylabel('z (km)')
    axis equal
    grid on
    nexttile
    plot3(X(1,1),X(1,2),X(1,3),'>','MarkerEdgeColor','green','MarkerFaceColor','green')
    hold on
    plot3(X(end,1),X(end,2),X(end,3),'x','MarkerEdgeColor','red','MarkerFaceColor','red')
    ellipsoid(0,0,0,P.R_deimos,P.R_deimos,P.R_deimos)
    plot3(X(:,1),X(:,2),X(:,3))
    xlabel('x (km)')
    ylabel('y (km)')
    zlabel('z (km)')
    legend('start', 'end', 'location', 'bestoutside')
    axis equal
    grid on
    
    frame2 = figure;
    tiledlayout(2,2);
    nexttile
    plot(X_f2(1,1),X_f2(1,2),'>','MarkerEdgeColor','green','MarkerFaceColor','green')
    hold on
    plot(X_f2(end,1),X_f2(end,2),'x','MarkerEdgeColor','red','MarkerFaceColor','red')
    ellipsoid(0,0,0,P.R_deimos,P.R_deimos,P.R_deimos)
    plot(X_f2(:,1),X_f2(:,2))
    xlabel('x (km)')
    ylabel('y (km)')
    title('Frame 2: IAU frame, rotating')
    axis equal
    grid on
    nexttile
    plot(X_f2(1,2),X_f2(1,3),'>','MarkerEdgeColor','green','MarkerFaceColor','green')
    hold on
    plot(X_f2(end,2),X_f2(end,3),'x','MarkerEdgeColor','red','MarkerFaceColor','red')
    ellipsoid(0,0,0,P.R_deimos,P.R_deimos,P.R_deimos)
    plot(X_f2(:,2),X_f2(:,3))
    xlabel('y (km)')
    ylabel('z (km)')
    axis equal
    grid on
    nexttile
    plot(X_f2(1,1),X_f2(1,3),'>','MarkerEdgeColor','green','MarkerFaceColor','green')
    hold on
    plot(X_f2(end,1),X_f2(end,3),'x','MarkerEdgeColor','red','MarkerFaceColor','red')
    ellipsoid(0,0,0,P.R_deimos,P.R_deimos,P.R_deimos)
    plot(X_f2(:,1),X_f2(:,3))
    xlabel('x (km)')
    ylabel('z (km)')
    axis equal
    grid on
    nexttile
    plot3(X_f2(1,1),X_f2(1,2),X_f2(1,3),'>','MarkerEdgeColor','green','MarkerFaceColor','green')
    hold on
    plot3(X_f2(end,1),X_f2(end,2),X_f2(end,3),'x','MarkerEdgeColor','red','MarkerFaceColor','red')
    ellipsoid(0,0,0,P.R_deimos,P.R_deimos,P.R_deimos)
    plot3(X_f2(:,1),X_f2(:,2),X_f2(:,3))
    xlabel('x (km)')
    ylabel('y (km)')
    zlabel('z (km)')
    legend('start', 'end', 'location', 'best')
    axis equal
    grid on
    
    frame3 = figure;
    tiledlayout(2,2);
    nexttile
    plot(X_f3(1,1),X_f3(1,2),'>','MarkerEdgeColor','green','MarkerFaceColor','green')
    hold on
    plot(X_f3(end,1),X_f3(end,2),'x','MarkerEdgeColor','red','MarkerFaceColor','red')
    ellipsoid(0,0,0,P.R_deimos,P.R_deimos,P.R_deimos)
    plot(X_f3(:,1),X_f3(:,2))
    xlabel('x (km)')
    ylabel('y (km)')
    title('Frame 3: IAU frame, fixed')
    axis equal
    grid on
    nexttile
    plot(X_f3(1,2),X_f3(1,3),'>','MarkerEdgeColor','green','MarkerFaceColor','green')
    hold on
    plot(X_f3(end,2),X_f3(end,3),'x','MarkerEdgeColor','red','MarkerFaceColor','red')
    ellipsoid(0,0,0,P.R_deimos,P.R_deimos,P.R_deimos)
    plot(X_f3(:,2),X_f3(:,3))
    xlabel('y (km)')
    ylabel('z (km)')
    axis equal
    grid on
    nexttile
    plot(X_f3(1,1),X_f3(1,3),'>','MarkerEdgeColor','green','MarkerFaceColor','green')
    hold on
    plot(X_f3(end,1),X_f3(end,3),'x','MarkerEdgeColor','red','MarkerFaceColor','red')
    ellipsoid(0,0,0,P.R_deimos,P.R_deimos,P.R_deimos)
    plot(X_f3(:,1),X_f3(:,3))
    xlabel('x (km)')
    ylabel('z (km)')
    axis equal
    grid on
    nexttile
    plot3(X_f3(1,1),X_f3(1,2),X_f3(1,3),'>','MarkerEdgeColor','green','MarkerFaceColor','green')
    hold on
    plot3(X_f3(end,1),X_f3(end,2),X_f3(end,3),'x','MarkerEdgeColor','red','MarkerFaceColor','red')
    ellipsoid(0,0,0,P.R_deimos,P.R_deimos,P.R_deimos)
    plot3(X_f3(:,1),X_f3(:,2),X_f3(:,3))
    xlabel('x (km)')
    ylabel('y (km)')
    zlabel('z (km)')
    legend('start', 'end', 'location', 'best')
    axis equal
    grid on
    
    OrbitalElements = figure;
    t = tiledlayout(7,1);
    nexttile
    plot(tDays, OE(:,1))
    ylabel('a')
    grid on
    nexttile
    plot(tDays, OE(:,2))
    ylabel('e')
    grid on
    nexttile
    plot(tDays, OE(:,3))
    ylabel('i')
    grid on
    nexttile
    plot(tDays, OE(:,4))
    ylabel('\omega')
    grid on
    nexttile
    plot(tDays, OE(:,5))
    ylabel('\Omega')
    grid on
    nexttile
    plot(tDays, real(OE(:,6)))
    ylabel('\nu')
    grid on
    nexttile
    plot(tDays, r)
    xlabel('t - t_{0} (days)')
    ylabel('r - R')
    grid on
    title(t, 'Orbital Elements (from frame 3, units of km and deg)')
    
    eccentricity = figure;
    plot(OE(1,2).*cosd(OE(1,4)), OE(1,2).*sind(OE(1,4)),'>','MarkerEdgeColor','green','MarkerFaceColor','green')
    hold on
    plot(OE(end,2).*cosd(OE(end,4)), OE(end,2).*sind(OE(end,4)),'x','MarkerEdgeColor','red','MarkerFaceColor','red')
    plot(OE(:,2).*cosd(OE(:,4)), OE(:,2).*sind(OE(:,4)))
    xlabel('e cos(\omega)')
    ylabel('e sin(\omega)')
    title('eccentricity vector evolution (t_{0} = May 01, 2023 03:00 PM CDT)')
    grid on
    
    acc = figure;
    semilogy(tDays, acc_D2)
    if P.mode == 1
        hold on
        semilogy(tDays, acc_Dns)
        semilogy(tDays, acc_M3)
        semilogy(tDays, acc_S3)
        legend('Deimos 2-body', 'Deimos oblate', 'Mars 3rd body', 'Sun 3rd body', 'Location', 'bestoutside')
    elseif P.mode == 2
        legend('Deimos 2-body', 'Location', 'bestoutside')
    elseif P.mode == 3
        hold on
        semilogy(tDays, acc_Dns)
        legend('Deimos 2-body', 'Deimos oblate', 'Location', 'bestoutside')
    elseif P.mode == 4
        hold on
        semilogy(tDays, acc_M3)
        semilogy(tDays, acc_S3)
        legend('Deimos 2-body', 'Mars 3rd body', 'Sun 3rd body', 'Location', 'bestoutside')
    elseif P.mode == 5
        hold on
        semilogy(tDays, acc_M3)
        legend('Deimos 2-body', 'Mars 3rd body', 'Location', 'bestoutside')
    end
    xlabel('t - t_{0} (days)')
    ylabel('acceleration (km^{2}/s)')
    title(['OE_{0} (km, deg) = ' num2str(OE_initial)])
    grid on
end