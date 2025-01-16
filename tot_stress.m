function [varargout]=tot_stress(taug,bases,vels,varargin)

% tot_stress.m - program to find the total deviatoric stress when given the
% gravitational stress, stress basis functions, and velocity.
%
% Based on the method of Flesch et al. (2001) but solving stress balance 
% equation in matrix form.
% Uses cartesian coordinates, assumes flat Earth but takes into account
% smaller spacing of lines of longitude at higher latitudes.
%
% Requires gravitational stress, basis functions, and velocities to be
% gridded according to the same rectangular grid.
%
% Written by Hamish Hirschberg 2016-17
%
% Input
% 1. m x n x 3 matlab variable of gravitationally-induced stresses
% 2. m x n x k x 3 matlab variable of basis stress functions
% AND EITHER
% 3. xyz-based file with East and North velocities in mm/yr
% OR
% 3. m x n x 3 matlab variable of strain rates in 10^-12/s
% WITH OPTIONAL (use [] to skip if using later options;
%   requires xyz-based velocity file)
% 4. name of output xyz-based file of total deviatoric stress
% 5. name of output xyz file of misfit
% 6. name of output xyz file of effective viscosity
% 7. name of output xyz-based file of boundary stress
% 8. name of output xyz-based file of strain rates
% AND OPTIONAL (use [] to skip if using later options)
% 9. minimum strain rate second invariant in 10^-12/s considered reliable
%   (default: 0, i.e. no minimum; ignored by matlab output)
% 10. maximum misfit considered reliable (default: 1, i.e. no maximum)
%
% Output
% 1. m x n x 3 matlab variable of total deviatoric stress
% 2. m x n matrix of misfit
% 3. m x n matrix of effective viscosity
% 4. m x n x 3 matlab variable of boundary stress
% 5. m x n x 3 matlab variable of strain rates
%
% Output text file format for stresses and strain rates:
% Column 1: Longitude
% Column 2: Latitude
% Column 3: Most extensional stress/strain rate (SHmin)
% Column 4: Most compressive stress/strain rate (SHmax)
% Column 5: Azimuth of SHmax
% Column 6: Magnitude of stress/strain rate - sqrt(sum(s_ij*s_ij)/2)

[ny,nx,~]=size(taug);               % dimensions of region
mine=[];                % initialise variables
tol=[];

% if velocities input as xyz-based file
if ischar(vels)
    rE=6371;                            % radius of the Earth
    yearsec=31556736;                   % number of seconds in a year
    data=dlmread(vels);                         % read in file
    dlo=data(2,1)-data(1,1);               % longitude spacing
    dla=data(nx+1,2)-data(1,2);                  % latitude spacing
    
    dy=rE*dla*pi/180;                           % colatitude spacing in km
    dx=zeros(ny,nx);                             % azimuthal spacing in km
    u=zeros(ny,nx);                    % matrix of x-velocity
    v=zeros(ny,nx);                     % matrix of y-velocity
    
    % create matrix of GPE and dx
    for ii=1:ny
        for j=1:nx
            k=(ii-1)*nx+j;          % index of point
            % account for smaller longitude spacing at higher latitudes
            dx(ii,j)=rE*cosd(data(k,2))*dlo*pi/180;
            u(ii,j)=data(k,3)/yearsec*1e9;      % read in x-velocity
            v(ii,j)=data(k,4)/yearsec*1e9;      % read in y-velocity
        end
    end

    % Calculate strain rates from gradients of velocities
    epsd=zeros(ny,nx,3);
    [ux,uy]=gradient(u);
    [vx,vy]=gradient(v);
    epsd(:,:,1)=ux./dx;
    epsd(:,:,2)=vy./dy;
    epsd(:,:,3)=(vx./dx+uy./dy)/2;

% else strain rates given as matlab variables
else
    epsd=vels;
end

% optional output files
if nargin > 3
    ntot=nx*ny;                         % total number of points
    otaut=varargin{1};                  % total stress output file
    mtaut=zeros(ntot,6);
    if nargin > 4
        omis=varargin{2};               % misfit output file
        mmis=zeros(ntot,3);
        if nargin > 5
            ovisc=varargin{3};          % viscosity output file
            mvisc=zeros(ntot,3);
            if nargin > 6
                otaub=varargin{4};      % boundary stress output file
                mtaub=zeros(ntot,6);
                if nargin > 7
                    oepsd=varargin{5};  % strain rate output file
                    mepsd=zeros(ntot,6);
                    if nargin > 8
                        mine=varargin{6};       % set min. strain rate
                        if nargin > 9
                            tol=varargin{7};    % set max. misfit
                        end
                    end
                end
            end
        end
    end
end

% default values
if isempty(mine)
    mine=0;             % minimum strain rate magnitude (i.e. all strain rates)
end
if isempty(tol)
    tol=1;              % maximum misfit (i.e. all misfits)
end

coeffs=ones(size(bases,3),1);      % create vector of coefficients for basis functions
%options=optimset('MaxIter',1e5,'MaxFunEvals',1e6,'TolFun',1e-9,'TolX',1e-8);
options=optimset('MaxIter',1e5,'MaxFunEvals',1e6,'TolFun',1e-7,'TolX',1e-6);
func=@(coeffs)totlincom(taug,bases,epsd,coeffs);        % set up function
coeffs=fminunc(func,coeffs,options);            % minimise the function

taub=zeros(size(taug));     % set up boundary stresses
for ii=1:length(coeffs)
    taub=taub+squeeze(bases(:,:,ii,:)*coeffs(ii));   % boundary stresses as sum of basis functions
end
taut=taug+taub;      % total stress as sum of gravitational and boundary stresses

% calculate misfit between stress and strain rate
% E, T, etau - indetermediary functions
E=sqrt(epsd(:,:,1).^2+epsd(:,:,1).*epsd(:,:,2)+epsd(:,:,2).^2+epsd(:,:,3).^2);
T=sqrt(taut(:,:,1).^2+taut(:,:,1).*taut(:,:,2)+taut(:,:,2).^2+taut(:,:,3).^2);
etau=2*epsd(:,:,1).*taut(:,:,1)+taut(:,:,1).*epsd(:,:,2)+taut(:,:,2).*epsd(:,:,1)+...
    2*epsd(:,:,2).*taut(:,:,2)+2*taut(:,:,3).*epsd(:,:,3);

mis=(1-etau./(2*E.*T))/2;     % calculate misfit
eta=0.5*T./E;               % calculate effective viscosity

if nargout>=1
    varargout{1}=taut;                      % output total stress
    if nargout>=2
        varargout{2}=mis;                   % output misfit
        if nargout>=3
            varargout{3}=eta;              % output viscosity
            if nargout>=4
                varargout{4}=taub;          % output boundary stress
                if nargout>=5
                    varargout{5}=epsd;      % output strain rate
                    if nargout>=6
                        varargout{6}=u;
                        varargout{7}=v;
                    end
                end
            end
        end
    end
end

bad=(isnan(mis)|mis>tol|E<mine);     % ignore points that have a misfit greater than tol

% calculate principal axes of boundary stress
paazb=princaxes(taub);
paazb(repmat(bad,[1,1,3]))=0;
B=sqrt(taub(:,:,1).^2+taub(:,:,1).*taub(:,:,2)+taub(:,:,2).^2+taub(:,:,3).^2);

B(bad)=NaN;       % ignore points if misfit too large
eta(bad)=NaN;
T(bad)=NaN;

% calculate principal axes of strain rate
paaze=princaxes(epsd);
paaze(repmat(bad,[1,1,3]))=0;      % convert unreliable values to 0

% calculate principal axes of total deviatoric stress
paazt=princaxes(taut);
paazt(repmat(bad,[1,1,3]))=0;      % convert NaN to 0
 
% output total stress file
if ~isempty(otaut)
    mtaut(:,1:2)=data(:,1:2);
    for ii=1:ny
        for j=1:nx
            k=(ii-1)*nx+j;
            mtaut(k,3)=paazt(ii,j,1);       % most extensional stress
            mtaut(k,4)=paazt(ii,j,2);       % most compressional stress
            mtaut(k,5)=paazt(ii,j,3);       % SHmax azimuth
            mtaut(k,6)=T(ii,j);     % second invariant
        end
    end
    dlmwrite(otaut,mtaut,'\t');             % write file to otaut
end
 
% output misfit file
if ~isempty(omis)
    mmis(:,1:2)=data(:,1:2);
    for ii=1:ny
        for j=1:nx
            k=(ii-1)*nx+j;
            mmis(k,3)=mis(ii,j);            % misfit
        end
    end
    dlmwrite(omis,mmis,'\t');
end

% output viscosity file
if ~isempty(ovisc)
    mvisc(:,1:2)=data(:,1:2);
    for ii=1:ny
        for j=1:nx
            k=(ii-1)*nx+j;
            mvisc(k,3)=eta(ii,j);           % effective viscosity
        end
    end
    dlmwrite(ovisc,mvisc,'\t');
end

% output boundary stress file
if ~isempty(otaub)
    mtaub(:,1:2)=data(:,1:2);
    for ii=1:ny
        for j=1:nx
            k=(ii-1)*nx+j;
            mtaub(k,3)=paazb(ii,j,1);       % same order as total stress
            mtaub(k,4)=paazb(ii,j,2);
            mtaub(k,5)=paazb(ii,j,3);
            mtaub(k,6)=B(ii,j);
        end
    end
    dlmwrite(otaub,mtaub,'\t');
end

% output strain rate file
if ~isempty(oepsd)
    mepsd(:,1:2)=data(:,1:2);
    for ii=1:ny
        for j=1:nx
            k=(ii-1)*nx+j;
            mepsd(k,3)=paaze(ii,j,1);       % same order as total stress
            mepsd(k,4)=paaze(ii,j,2);
            mepsd(k,5)=paaze(ii,j,3);
            mepsd(k,6)=E(ii,j);
        end
    end
    dlmwrite(oepsd,mepsd,'\t');
end

end
