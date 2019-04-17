function [varargout]=moho2gpe(elev,varargin)

% moho2gpe.m - program to determine the GPE from Moho depths/crustal
% thicknesses and topography for given crustal and mantle densities.
%
% Requires input Moho and crustal density to be at the same grid points as
% topography or to be scalars.
%
% Written by Hamish Hirschberg 2016-17
% 
% Input variables: (use [] to skip if using later options; 1-6 optional)
% 0. topography in m as xyz file name or matlab variable (required)
% 1. output name of output xyz file (requires >=1 input xyz file)
% 2. depth in km GPE averaged down to as a scalar
% 3. crustal density in g/cm3 as scalar, xyz file name, or matlab variable
% 4. Moho depth in km as xyz file name or matlab variable
% 5. mantle density in g/cm3 as scalar
% 6. points at sea (1) vs points on land (0) as xyz file or matlab variable
%    (use if region contains land below sea level)
% 
% Output variables (optional):
% 1. GPE in MPa as matlab variable

% intialise some variables
data=[];
output=[];
depth=[];
rhoc=[];
moho=[];
rhom=[];
water=[];

% read in variables
% read in topography
if ischar(elev)
    data=dlmread(elev);      % if topography from a file
    topo=data(:,3)/1000;
else
    topo=elev/1000;           % if topography from matlab variable
end
if nargin>1
    output=varargin{1};     % output file name
    if nargin>2
        depth=varargin{2};      % read in depth
        if nargin>3
            % read in crustal density
            if ischar(varargin{3})
                data=dlmread(varargin{3});      % if crustal density from file
                rhoc=data(:,3);
            else
                rhoc=varargin{3};  % if crustal density from scalar or matlab variable
            end
            if nargin>4
                % read in Moho depth
                if ischar(varargin{4})
                    data=dlmread(varargin{4});      % if Moho depth from file
                    moho=data(:,3);
                else
                    moho=varargin{4};  % if Moho depth from matlab variable
                end
                if nargin>5
                    rhom=varargin{5};       % set mantle density
                    if nargin>6
                        % read in submarine points
                        if ischar(varargin{6})
                            data=dlmread(varargin{6});      % if points from file
                            water=data(:,3);
                        else
                            water=varargin{6};  % if points from matlab variable
                        end
                    end
                end
            end
        end
    end
end

% insert defaults where necessary
if isempty(depth)
    depth=25;       % 25 km
end
if isempty(rhoc)
    rhoc=2.76;      % 2.76 g/cm3
end
if isempty(moho)
    moho=depth;     % equivalent to considering only crust
elseif max(abs(moho))>1000
    % warn if Moho depth appears to have values too large to be in km
    warning('Moho2GPE:mohoM',...
        'Moho depth appears to be in m. This program assumes Moho depth in km.')
end
if isempty(rhom)
    rhom=3.2;       % 3.2 g/cm3
end
if isempty(water)
    water=topo<0;   % assume no land below sea level
end

if max(abs(topo))<20
    % warn if topography appears to have values too small to be in m
    warning('Moho2GPE:topoKm',...
        'Topography appears to be in km. This program assumes topography in m.')
end

% set Earthly constants
g=9.81;         % gravity
rhow=1;         % density of water

% additional constants for datum
if length(rhoc)==1
    rhod=rhoc;
else
    rhod=2.76;
end
if length(moho)==1
    mohod=moho;
else
    mohod=min(25,depth);
end

% account for moho depths greater than averaging depth
moho(moho>depth)=depth;

% calculate datum
gped=g*(rhod*mohod^2+2*rhod*mohod*(depth-mohod)+rhom*(depth-mohod)^2)/(2*depth);

% calculate GPE
gpei=g*(rhoc.*(moho+topo).^2/2+rhoc.*(moho+topo).*(depth-moho)+...
    rhom*(depth-moho).^2/2)/depth;      % basic calculation
gpew=water.*rhow*g.*(topo*depth+topo.^2/2)/depth;   % account for water
gpe=gpei-gpew-gped;           % remove datum

% write out data to file named output only if xyz file inputted
if ~isempty(output) && ~isempty(data)
    dlmwrite(output,[data(:,1:2),gpe],'\t');
end

if nargout>=1
    varargout{1}=gpe;       % output GPE to matlab if requested
end

end