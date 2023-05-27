%Function that returns orbital elements of an orbit based on the inputs of
%the position and velocity vectors
function OE = vec2oe(mu, rVec, vVec)

    %Calculate semi-major axis
    r = norm(rVec);
    v = norm(vVec);
    a = -mu*r/(v^2*r - 2*mu);
    OE(1) = a;
    %Calculate angular moment vector
    hVec = cross(rVec,vVec);
    %Calculate e vector
    eVec = cross(vVec, hVec)/mu - rVec/r;
    %Calculate magnitude of e
    OE(2) = norm(eVec);
    %Calculate inclination
    OE(3) = acosd(dot([0; 0; 1], hVec)/norm(hVec));
    %Calculate vector n
    n = cross([0; 0; 1], hVec)/norm(cross([0; 0; 1], hVec));
    %Calculate argument of periapsis
    if eVec(3) > 0
        OE(4) = acosd(dot(n, eVec)/norm(eVec));
    else
        OE(4) = -acosd(dot(n, eVec)/norm(eVec));
    end
    %Calculate longitude of ascending node
    OE(5) = atan2d(n(2), n(1));
    %Calculate true anomaly
    if dot(rVec, vVec) > 0
        OE(6) = acosd(dot(rVec, eVec)/(norm(rVec)*norm(eVec)));
    else
        OE(6) = -acosd(dot(rVec, eVec)/(norm(rVec)*norm(eVec)));
    end
end