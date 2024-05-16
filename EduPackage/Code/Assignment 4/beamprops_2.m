%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% define a function beamprops_2 to solve the problem
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
function [k_iz, k_Gz] = beamprops_2(energy, mass, sctAngle, G)
    % PRECONDTION: − energy in meV
    %              − mass in atomic mass units
    %             − sctAngle is known and in radians
    %             − G is a vector of (moduli of) allowed reciprocal lattice vectors
    %             − inAngle is variable
    %             − Scattering is elastic and the Laue condition is satisfied
    % POSTCONDITION: − k_iz is a scalar
    %                − k_Gz is a vector where each element k_Gz( i) is an allowed component corresponding to reciprocal vector G(i)

    % Reduced Planck's Constant in units such that k_i is in A ^{−1}
    hbar = 2.05;

    % Incident beam wavevector /A^{−1}
    k_i = sqrt(2*mass*energy)/hbar;

    % Laue Condition || surface + elastic collision
    theta_i = acos(G/(1.34085*k_i)) - 0.73478;

    % Perpendicular component
    k_iz = -k_i*cos(theta_i);
    % Elastic collision => |k_i|=|k_G| for all G.
    k_Gz = k_i*cos(sctAngle-theta_i);
    end