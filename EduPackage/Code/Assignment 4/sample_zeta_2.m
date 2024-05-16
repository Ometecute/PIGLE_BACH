function [zeta] = sample_zeta_2(R,h,G1,G2,G3)
    %Function that produces corrugation function with three Fourier components
    %Inputs:
    %   - R: positions to sample (Nx2 matrix, where N is number of positions)
    %   - h: corrugation height (scalar)
    %   - G1, G2, G3, : recpiprocal lattice vectors (column vectors)
    
    %Output:
    %   - zeta: corrugation function at each position (Nx1 matrix)
    
    zeta_0 = cos(R*G1) + cos(R*G2) + cos (R*G3); %Unscaled Fourier series
    zeta_0 = zeta_0 - min(min(zeta_0)); %So that zeta is always greater than 0
    zeta_0 = zeta_0/max(max(zeta_0)); %Giving that zeta has max value of 1
    
    zeta = h*zeta_0; %Giving zeta a max value of h
end