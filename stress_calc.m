% stress_calc.m - script to calculate total deviatoric stress from crustal
% structure and geodetic velocities.
%
% Written by Hamish Hirschberg 2017
% 
% This scripts calls a sequence of programs to calculate total deviatoric
% stress, including gravitational and boundary stresses, from topography,
% crustal structure and geodetic velocities. The codes were written as part
% of my Master's project and are based on the method of Flesch et al.
% (2001).
% 
% The programs called by this script are:
% moho2gpe.m
% grav_stress.m
% basis_fns.m
% tot_stress.m
%
% Set unwanted output file names to []

model='_100';            % model suffix; set to '' for no suffix
% directory for output; set to '' or './' for this directory
wdir='./';

% inputs to calculate GPE
% if GPE has already been calculated, set topo=[] and gpe below to file name
topo='topo.xyz';        % input topography file (m)
moho='moho.xyz';        % input moho depth file (km)
depth=100;              % depth averaged over
rhoc=2.81;              % crustal density
rhom=3.3;               % mantle density
sea=[];                 % offshore points file. Only necessary if region
                        % contains dry land below sea level

% gravitational stress files
gpe=strcat(wdir,'gpe',model,'.xyz');          % GPE file (MPa)
gravs=strcat(wdir,'grav',model,'.xyz');       % gravitational stress file

% total stress input and output
% if only calculating gravitational stress, set vel=[]
vel='vel.xyz';          % velocity file (mm/yr)
tot=strcat(wdir,'tot',model,'.xyz');          % total deviatoric stress file
mis=strcat(wdir,'mis',model,'.xyz');          % misfit file
visc=strcat(wdir,'visc',model,'.xyz');        % effective viscosity (10^21 Pa.s)
bound=strcat(wdir,'bound',model,'.xyz');      % boundary stress
sr=strcat(wdir,'sr',model,'.xyz');            % strain rate (10^-12/s)
mine=0.2;                 % minimum strain rate used (10^-12/s)
maxm=1;                 % maximum misfit used

if ~isempty(topo)
    % use moho2gpe.m to find gpe down to 25 km with a crustal density of 2.8
    moho2gpe(topo,gpe,depth,rhoc,moho,rhom,sea);
end

% use grav_stress.m to calculate gravitational stress from GPE
[grav,Gamma]=grav_stress(gpe,gravs);

if ~isempty(vel)
    % use basis_fns.m to find the stress basis functions
    bases=basis_fns(gpe);
    % use tot_stress.m to find the total stress
    tot_stress(grav,bases,vel,tot,mis,visc,bound,sr,mine,maxm);
    
    % plot using GMT
    system(['./map.txt ' model]);
    % plotting script written in bash for sample data
    % will need amendments for other data or scripting languages
end













