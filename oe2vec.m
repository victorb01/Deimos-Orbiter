function X = oe2vec(OE, mu)
    R1 = [cosd(OE(5)), -sind(OE(5)), 0;
         sind(OE(5)), cosd(OE(5)), 0;
         0, 0, 1];
    R2 = [1, 0, 0;
          0, cosd(OE(3)), -sind(OE(3));
          0, sind(OE(3)), cosd(OE(3))];
    R3 = [cosd(OE(4)), -sind(OE(4)), 0;
         sind(OE(4)), cosd(OE(4)), 0;
         0, 0, 1];
    R = R1*R2*R3;

    p = OE(1)*(1-OE(2)^2);
    r = p/(1+OE(2)*cosd(OE(6)));

    r_pqw = r*[cosd(OE(6)); sind(OE(6)); 0];
    v_pqw = [-sind(OE(6))*sqrt(mu/p); sqrt(mu/p)*(OE(2) + cosd(OE(6))); 0];

    X = [R*r_pqw; R*v_pqw];
end