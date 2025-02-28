function [u, tau_yield] = Schoof2006_icestream( A, H, tantheta, L, m, y)

ice_density = 910;
grav = 9.81;

% Calculate the gravitational driving stress f
f = -ice_density * grav * H * tantheta;

% Calculate the ice hardness factor B
B = A^(-1/3);

% Calculate the "ice stream half-width" W
W = L * (m+1)^(1/m);

% Calculate the till yield stress across the stream
tau_yield = f * abs( y/L).^m;

% Calculate the analytical solution for u
ua = -2 * f.^ 3 * L.^ 4 / (B.^ 3 * H.^ 3);
ub = ( 1 / 4                    ) * (   ( y/L).^      4  - (m+1).^ (   4/m) );
uc = (-3 / ((m+1)     * (  m+4))) * (abs( y/L).^ (  m+4) - (m+1).^ (1+(4/m)));
ud = ( 3 / ((m+1).^ 2 * (2*m+4))) * (abs( y/L).^ (2*m+4) - (m+1).^ (2+(4/m)));
ue = (-1 / ((m+1).^ 3 * (3*m+4))) * (abs( y/L).^ (3*m+4) - (m+1).^ (3+(4/m)));
u = ua * (ub + uc + ud + ue);

% Outside the ice-stream, velocity is zero
u( abs(y) > W) = 0;

end