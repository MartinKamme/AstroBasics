classdef BasicAstrodynamics
    properties
        statein     % [6x1] column [r,v]
        coes        % [6x1] column [h; inc (rad) ; RAAN (rad) ; ecc ; per (rad) ; theta (rad)]
        mu
    end
    methods
        
        function [R,V] = coes2rv(obj)
            h       = obj.coes(1);
            inc     = obj.coes(2);
            RAAN    = obj.coes(3);
            e       = obj.coes(4);
            per     = obj.coes(5);
            theta   = obj.coes(6);
            
            % h     [km^2/s] Specific angular momentum
            % i     [rad] Inclination
            % RAAN  [rad] Right ascension (RA) of the ascending node
            % e     Eccentricity
            % per   [rad] Argument of perigee
            % theta [rad] True anomaly
            
            % State Vectors in Perifocal coordinates
            rx = h^2/obj.mu*(1/(1 + e*cos(theta)))*[cos(theta);sin(theta);0];
            vx = obj.mu/h*[-sin(theta); (e +cos(theta));0];
            
            % Direction cosine matrix
            DCM = [cos(per), sin(per),0;-sin(per),cos(per),0;0,0,1]*...
                [1,0,0;0,cos(inc),sin(inc);0,-sin(inc),cos(inc)]*...
                [cos(RAAN), sin(RAAN),0;-sin(RAAN),cos(RAAN),0;0,0,1];
            
            % Transformation Matrix
            Dcm = inv(DCM);
            
            % ECI R
            R = Dcm*rx;
            
            % ECI V
            V = Dcm*vx;
            
        end
        
        function [a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(obj)
            %----------------------- Begin Code Sequence -----------------------------%
            % Purpose:                                                                %
            % Convert a given set of state vectors in ECI reference frame to orbital  %
            % elements.                                                               %
            %-------------------------------------------------------------------------%
            %                                                                         %
            % Inputs:                                                                 %
            %--------
            %r_ECI                  [3 x N]                         Position Vector in
            %                                                       ECI coordinate
            %                                                       frame of reference
            %
            %v_ECI                  [3 x N]                         Velocity vector in
            %                                                       ECI coordinate
            %                                                       frame of reference
            %
            %mu                     double                          Gravitational Constant
            %                                                       Defaults to Earth if
            %                                                       not specified
            % Outputs:
            %---------                                                                %
            %a                      [1 x N]                         Semi-Major Axis
            %                                                       (km)
            %
            %eMag                   [1 x N]                         Eccentricity
            %                                                       (unitless)
            %
            %i                      [1 x N]                         inclination
            %                                                       (radians)
            %
            %O                      [1 x N]                         Right Ascention of
            %                                                       the ascending node
            %                                                       (radians)
            %
            %o                      [1 x N]                         Argument of perigee
            %                                                       (radians)
            %
            %nu                      [1 x N]                         Mean Anomaly
            %                                                       (radians)
            %
            %truLon                 [1 x N]                         True Longitude
            %                                                       (radians)
            %
            %argLat                 [1 x N]                         Argument of Latitude
            %                                                       (radians)
            %
            %lonPer                 [1 x N]                         Longitude of Periapse
            %                                                       (radians)
            %
            %p                      [1 x N]                         Semilatus Rectum
            %                                                       (km)
            %
            % References:
            %-------------
            %Vallado,D. Fundamentals of Astrodynamics and Applications. 2007.
            %
            % Function Dependencies:
            %------------------
            %None
            %------------------------------------------------------------------       %
            % Programed by Darin Koblick  03-04-2012                                  %
            % Updated to address circular equatorial orbits       12/12/2013          %
            %------------------------------------------------------------------       %
            
            r = obj.statein(1:3,1);
            v = obj.statein(4:6,1);
            %if ~exist('mu','var');  t = getConst(); mu = t.Earth.Mu; end
            %Specific angular momentum
            h = cross(r,v);
            %h = h'; % mac does some whack stuff with cross
            n = cross(repmat([0;0;1],[1,size(r,2)]),h);
            %n = n';  % mac does some whack stuff with cross
            nMag = sqrt(sum(n.^2,1));
            vMag = sqrt(sum(v.^2,1));
            rMag = sqrt(sum(r.^2,1));
            hMag = sqrt(sum(h.^2,1));
            e = (1./obj.mu).*(bsxfun(@times,(vMag.^2 - obj.mu./rMag),r) - bsxfun(@times,dot(r,v),v));
            eMag = sqrt(sum(e.^2,1));
            zeta = (vMag.^2)./2 - obj.mu./rMag;
            
            %Special Procedure when we have a parabolic orbit
            idx = eMag ~= 1;
            a = NaN(size(eMag));
            p = NaN(size(eMag));
            if any(idx)
                a(idx) = -obj.mu./(2.*zeta(idx));
                p = a(idx).*(1-eMag(idx).^2);
            else
                a(idx) = Inf;
                p(idx) = (hMag(idx).^2)./obj.mu;
            end
            
            %Compute the angles
            i = acos(h(3,:)./hMag);
            O = acos(n(1,:)./nMag);
            o = acos(dot(n,e)./(nMag.*eMag));
            nu = acos(dot(e,r)./(eMag.*rMag));
            lonPer = acos(e(1,:)./eMag);
            argLat = acos(dot(n,r)./(nMag.*rMag));
            truLon = acos(r(1,:)./rMag);
            
            %Account for those cases where satellite is in circular orbit
            O(n(1,:) == 0) = 0;
            o(dot(n,e) == 0) = 0;
            lonPer(e(1,:) == 0) = 0;
            nu(dot(e,r) == 0) = 0;
            argLat(dot(n,r) == 0) = 0;
            
            %Apply Quadrant Checks to All Determined Angles
            idx = n(2,:) < 0; if any(idx);  O(idx) = 2*pi - O(idx);  end
            idx = e(3,:) < 0; if any(idx); o(idx) = 2*pi - o(idx); end
            idx = dot(r,v) < 0; if any(idx); nu(idx) = 2*pi - nu(idx); end
            idx = e(2,:) < 0; if any(idx); lonPer(idx) = 2*pi-lonPer(idx);  end
            idx = r(3,:) < 0; if any(idx); argLat(idx) = 2*pi - argLat(idx); end
            idx = r(2,:) < 0; if any(idx); truLon(idx) = 2*pi - truLon(idx); end
        end
        
         function [statedot] = orbitpropagator(obj,t,state)
            % All vectors must be input as COLUMNS
            
            % state includes the following:
            % state(1) = rAx eci
            % state(2) = rAy eci
            % state(3) = rAz eci
            % state(4) = vAx eci
            % state(5) = vAy eci
            % state(6) = vAz eci
            % state(7) = deltarx
            % state(8) = deltary
            % state(9) = deltarz
            % state(10) = deltavx
            % state(11) = deltavy
            % state(12) = deltavz
            
            R = state(1:3);
            V = state(4:6);
            
            dx = V(1);
            dy = V(2);
            dz = V(3);
            
            r = norm(R);
            ddx = -obj.mu*R(1)/r^3;
            ddy = -obj.mu*R(2)/r^3;
            ddz = -obj.mu*R(3)/r^3;
            
            statedot =  [dx;dy;dz;ddx;ddy;ddz];  % Propagation of R and V vectors
            
            
        end
        
        function [tnew,stateout] = PropagateOrbit(obj,time,options)
            tspan = [0 time];
            state = obj.statein;
            [tnew,stateout] = ode45(@obj.orbitpropagator,tspan,state,options);
        end
        
       
    end
    
end
