function [ Lfun ] = totlincom( taug, bases, epsd, coeffs )
% totlincomb.m - calculates the functional that when minimised in terms of
% the coefficients coeffs, provides the linear combination of the basis
% functions bases, that when combined with the gravitational stress field
% taug, provides the closest fit to the direction of the stress field eps
%
% Inputs:
% taug - gravitational stress field
% bases - array of basis functions
% eps - the strain rate field
% coeffs - vector of coefficients with same number as third dimension of
% bases
%
% Output:
% Lfun - integral value to be minimised
%
% Total stress is sum of gravitational stress and boundary stress expressed
% as linear combination of basis functions
%
% Written by Hamish Hirschberg

sz=size(taug);
taut=taug;     % set up total stresses
for ii=1:length(coeffs)
    taut=taut+squeeze(bases(:,:,ii,:)*coeffs(ii));   % boundary stresses as sum of basis functions
end

% remove points from the boundary from consideration
epsd=epsd(2:sz(1)-1,2:sz(2)-1,:);
taut=taut(2:sz(1)-1,2:sz(2)-1,:);
% E, T, etau - indetermediary functions
E=sqrt((epsd(:,:,1).^2+epsd(:,:,1).*epsd(:,:,2)+epsd(:,:,2).^2+epsd(:,:,3).^2)*2);
T=sqrt((taut(:,:,1).^2+taut(:,:,1).*taut(:,:,2)+taut(:,:,2).^2+taut(:,:,3).^2)*2);
etau=2*epsd(:,:,1).*taut(:,:,1)+taut(:,:,1).*epsd(:,:,2)+taut(:,:,2).*epsd(:,:,1)+...
    2*epsd(:,:,2).*taut(:,:,2)+2*taut(:,:,3).*epsd(:,:,3);

Lmat=(E.*T-etau);        % matrix of the integral

Lfun=sum(Lmat(:));     % total sum

end
