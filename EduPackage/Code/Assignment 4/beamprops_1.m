%% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% define a function beamprops_1
% −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

function [k_iz, k_Gz] = beamprops_1(energy, mass, inAngle, G)
    % PRECONDTION: − energy in meV
    %              − mass in atomic mass units
    %              − inAngle is known and in radians
    %              − G is a vector of (moduli of) allowed reciprocal lattice vectors
    %              − outAngle is variable
    %              − Scattering is elastic and the Lauecondition is satisfied
    % POSTCONDITION: − k_iz is a scalar
    %                − k_Gz is a vector where each element k_Gz(i) is an allowed component corresponding to a reciprocal vector of modulus G(i)
    
    % Reduced Planck's Constant in units such that k_i is in A^{−1}
    hbar = 2.05;
    
    % Incident beam wavevector /A^{−1}
    k_i = sqrt(2*mass*energy)/hbar;
    
    % Laue Condition || surface + elastic collision => k_i(sin(outAngle)−sin(inAngle)) = G
    outAngle = asin(G/k_i + sin(inAngle));
    
    % Perpendicular component
    k_iz = -k_i*cos(inAngle);
    % Elastic collision => |k_i|=|k_G| for all G.
    k_Gz = k_i*cos(outAngle);
end