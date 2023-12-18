load NADH.mat;
global NADH_emission;
global NADH_excitation;
global NADH_wave;
NADH_emission=emission;
NADH_excitation=excitation;
NADH_wave=wave;

load FAD.mat;
global FAD_emission;
global FAD_excitation;
global FAD_wave;
FAD_emission=emission;
FAD_excitation=excitation;

% FAD_wave=wave;
load Keratin_WuandQu.mat;
Keratin_emission=emission;
Keratin_excitation=excitation;
Keratin_wave=wave;

%% Write out the text file

% Set the wavelength sampling you would like
NADH_wave = 200:5:800; 
% Here is one fluorophore, read in with that wavelength sampling.
NADH = fluorophoreRead('NADH','wave',NADH_wave);
% Here is the excitation emission matrix
NADH_eem = fluorophoreGet(NADH,'eem');
% FAD_eem_Nonzero = FAD_eem(62:81,2:60);


% Set the wavelength sampling you would like
FAD_wave = 200:5:800; 
% Here is one fluorophore, read in with that wavelength sampling.
FAD = fluorophoreRead('FAD','wave',FAD_wave);
% Here is the excitation emission matrix
FAD_eem = fluorophoreGet(FAD,'eem');
% FAD_eem_Nonzero = FAD_eem(62:81,2:60);

global FADorNADH

global QY;
global Excit_Wave;
global Emission_Wave;

global n;

global mus_air;
mus_air_= 1e-8; % [cm^-1]

global mus_FAD;
mus_FAD_ = 51.9; % [cm^-1] % mus is 204 at 523nm in ref. paper, combined with empirical model in prior work.
global mua_FAD;
global thickness_FAD;

global mus_NADH;
mus_NADH_ = 51.9; % [cm^-1] % mus is 204 at 523nm in ref. paper, combined with empirical model in prior work.
global mua_NADH;
global thickness_NADH;



 for thickness_FAD= (0.008:0.02:0.008)
    for thickness_NADH= (0.012:0.02:0.012)

Record_count=0;

FADorNADH = true;
Record_ExcitWave=zeros(20*3,1);
Record_EmissionWave=zeros(20*3,1);
Record_FluoReEmitPercent=zeros(20*3,1);
 for m = 49 %(25:3:37)%FAD(46:5:56) % FAD m range: (2:1:60)  % NADH: m range: (12:1:39)
     for n= 62:1:81 %(43:1:67) %FAD(62:1:81)  % FAD n range: (62:1:81)  % NADH: n=56;  n range: (43:1:67)  
      Record_count = Record_count +1;
    Record_ExcitWave(Record_count)=FAD_wave(1,m);
    Record_EmissionWave(Record_count)=FAD_wave(1,n);
   k=2;
   if FADorNADH==true    
           Excit_Wave=FAD_wave(1,m) %205 (2:1:60) % FAD: Excit_Wave=260nm
    Emission_Wave=FAD_wave(1,n) %505 (62:1:81) % FAD: Emission_Wave=550nm
    
    mus_air = mus_air_ * ((523 /  Excit_Wave)^0.6);
    mus_FAD = mus_FAD_ * ((523 /  Excit_Wave)^0.6); % [cm^-1]
    mua_FAD = FAD_excitation(m,1)/thickness_FAD;
     
    mus_NADH = mus_NADH_ * ((523 /  Excit_Wave)^0.6); % [cm^-1]
    mua_NADH = 3;

    mua_FAD
    mua_NADH
   QY=FAD_eem(n,m)/(1-10^(-1*FAD_excitation(m,1)))/k;
   %0.3854/(1-10^(-1))/2;
   %0.2124/(1-10^(-0.5949))/2;
   else 
           Excit_Wave=NADH_wave(1,m) %205 (2:1:60) % FAD: Excit_Wave=260nm
    Emission_Wave=NADH_wave(1,n) %505 (62:1:81) % FAD: Emission_Wave=550nm
    
    mus_air = mus_air_ * ((523 /  Excit_Wave)^0.6);
    mus_FAD = mus_FAD_ * ((523 /  Excit_Wave)^0.6); % [cm^-1]
    mua_FAD = 3;
     
    mus_NADH = mus_NADH_ * ((523 /  Excit_Wave)^0.6); % [cm^-1]
    mua_NADH = NADH_excitation(m,1)/thickness_NADH;
   QY=NADH_eem(n,m)/(1-10^(-1*NADH_excitation(m,1)))/k;
   %0.3854/(1-10^(-1))/2;
   %0.2124/(1-10^(-0.5949))/2;
   end
        
% B = G(n) %@func2; % [cm^-1d]

% for n=(1: 1: 121)

%% Geometry definition
MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;

model.G.nx                = 100; % Number of bins in the x direction
model.G.ny                = 100; % Number of bins in the y direction
model.G.nz                = 100; % Number of bins in the z direction
model.G.Lx                = 0.03; % [cm] x size of simulation cuboid
model.G.Ly                = 0.03; % [cm] y size of simulation cuboid
model.G.Lz                = 0.003+thickness_FAD+thickness_NADH; %0.0593; % [cm] z size of simulation cuboid

model.G.mediaPropertiesFunc = @mediaPropertiesFunc; % Media properties defined as a function at the end of this file
model.G.geomFunc          = @geometryDefinition; % Function to use for defining the distribution of media in the cuboid. Defined at the end of this m file.

 model = plot(model,'G'); 

%% Monte Carlo simulation
model.MC.useAllCPUs               = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
model.MC.simulationTimeRequested  = .1; % [min] Time duration of the simulation
model.MC.nExamplePaths            = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes

model.MC.matchedInterfaces        = true; % Assumes all refractive indices are the same
model.MC.boundaryType             = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.MC.wavelength               = Excit_Wave; % [nm] Excitation wavelength, used for determination of optical properties for excitation light

model.MC.lightSource.sourceType   = 2; % 0: Pencil beam, 1: Isotropically emitting line or point source, 2: Infinite plane wave, 3: Laguerre-Gaussian LG01 beam, 4: Radial-factorizable beam (e.g., a Gaussian beam), 5: X/Y factorizable beam (e.g., a rectangular LED emitter)
model.MC.lightSource.theta        = 0; % [rad] Polar angle of beam center axis
model.MC.lightSource.phi          = 0; % [rad] Azimuthal angle of beam center axis

model.MC.useLightCollector        = true;

model.MC.lightCollector.x         = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.MC.lightCollector.y         = 0; % [cm] y position
model.MC.lightCollector.z         = 0.03; % [cm] z position

model.MC.lightCollector.theta     = 0; % [rad] Polar angle of direction the light collector is facing
model.MC.lightCollector.phi       = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

model.MC.lightCollector.f         = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.MC.lightCollector.diam      = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.MC.lightCollector.fieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
model.MC.lightCollector.NA        = 0.22; % [-] Fiber NA. Only used for infinite f.

model.MC.lightCollector.res       = 50; % X and Y resolution of light collector in pixels, only used for finite f

model = runMonteCarlo(model);
 model = plot(model,'MC');

%% Fluorescence Monte Carlo
model.FMC.useAllCPUs              = true; % If false, MCmatlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
model.FMC.simulationTimeRequested = .1; % [min] Time duration of the simulation
model.FMC.nExamplePaths           = 100; % (Default: 0) This number of photons will have their paths stored and shown after completion, for illustrative purposes

model.FMC.matchedInterfaces       = true; % Assumes all refractive indices are the same
model.FMC.boundaryType            = 1; % 0: No escaping boundaries, 1: All cuboid boundaries are escaping, 2: Top cuboid boundary only is escaping, 3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
model.FMC.wavelength              = Emission_Wave; % [nm] Fluorescence wavelength, used for determination of optical properties for fluorescence light

model.FMC.useLightCollector       = true;

model.FMC.lightCollector.x        = 0; % [cm] x position of either the center of the objective lens focal plane or the fiber tip
model.FMC.lightCollector.y        = 0; % [cm] y position
model.FMC.lightCollector.z        = 0.03; % [cm] z position

model.FMC.lightCollector.theta    = 0; % [rad] Polar angle of direction the light collector is facing
model.FMC.lightCollector.phi      = pi/2; % [rad] Azimuthal angle of direction the light collector is facing

model.FMC.lightCollector.f        = .2; % [cm] Focal length of the objective lens (if light collector is a fiber, set this to Inf).
model.FMC.lightCollector.diam     = .1; % [cm] Diameter of the light collector aperture. For an ideal thin lens, this is 2*f*tan(asin(NA)).
model.FMC.lightCollector.fieldSize = .1; % [cm] Field Size of the imaging system (diameter of area in object plane that gets imaged). Only used for finite f.
model.FMC.lightCollector.NA       = 0.22; % [-] Fiber NA. Only used for infinite f.

model.FMC.lightCollector.res      = 50; % X and Y resolution of light collector in pixels, only used for finite f


model = runMonteCarlo(model,'fluorescence');

model = plot(model,'FMC');

P_exc_abs = model.G.dx*model.G.dy*model.G.dz*sum(model.MC.NA(:));
P_flu_emit = model.G.dx*model.G.dy*model.G.dz*sum(model.FMC.sourceDistribution(:));

Record_FluoReEmitPercent(Record_count)=P_flu_emit/P_exc_abs;

if FADorNADH==true
    fprintf('FAD\n');
else
    fprintf('NADH\n');
end
fprintf('%.3g% P_flu_emit.\n',P_flu_emit);
fprintf('%.3g% P_exc_abs.\n',P_exc_abs);
    end
end
    end
end


%% Geometry function(s) (see readme for details)
function M = geometryDefinition(X,Y,Z,parameters)
global thickness_FAD;
   zSurface = 0.003; % 0.003
   M = ones(size(X)); % Air/Keratin
   M(Z > zSurface) = 2; % "FAD" tissue
   M(Z > zSurface+thickness_FAD) = 3; % "NADH" tissue
%   cylinderradius  = 0.0100;
%   M = ones(size(X)); % fill background with fluorescence absorber
%   M(Y.^2 + (Z - 3*cylinderradius).^2 < cylinderradius^2) = 2; % fluorescer
end

%% Media Properties function (see readme for details)
function mediaProperties = mediaPropertiesFunc(parameters)
  mediaProperties = MCmatlab.mediumProperties;

global QY;
global mus_air;
global mus_FAD;
global mua_FAD;
global mus_NADH;
global mua_NADH;
global FADorNADH;

% global Keratin_QY;
   
  j=1;
  mediaProperties(j).name  = 'Air';
  mediaProperties(j).mua = 1e-8;
  mediaProperties(j).mus = mus_air; % [cm^-1]
  mediaProperties(j).g   = 1;
%   mediaProperties(j).QY = Keratin_QY; %@func2; % [cm^-1]
  
  j=2;
  mediaProperties(j).name  = 'Epithelium_FAD';
  mediaProperties(j).mua = mua_FAD; %127.27;
  mediaProperties(j).mus = mus_FAD; % 1; % [cm^-1]
  mediaProperties(j).g   = 0.9;
  
  if FADorNADH==true
  mediaProperties(j).QY = QY; %@func2; % [cm^-1]
  end 
  
  j=3;
  mediaProperties(j).name  = 'Epithelium_NADH';
  mediaProperties(j).mua = mua_NADH; %127.27;
  mediaProperties(j).mus = mus_NADH; % 1; % [cm^-1]
  mediaProperties(j).g   = 0.9;
  
  if FADorNADH==false
  mediaProperties(j).QY = QY; %@func2; % [cm^-1]
  end
%   mediaProperties(j).QY   = 0.4; % Fluorescence quantum yield



%    j=1;
%    mediaProperties(j).name  = 'Air';
%    mediaProperties(j).mua   = 100; % [cm^-1]
%    mediaProperties(j).mus   = 100; % [cm^-1]
%    mediaProperties(j).g     = 1;
% 
%      
%    j=2;
%    mediaProperties(j).name  = 'FAD';
%    mediaProperties(j).mua   = 100; %mua_ep1; % [cm^-1]
%    mediaProperties(j).mus   = 127.47; % [cm^-1]
%   mediaProperties(j).g     = 0.7; % 0.99; % 1; % 0.90;
% 
%     M = FAD_emission./FAD_excitation; % [cm^-1]
%    mediaProperties(j).ES =0; %FAD_emission(n);  
%     mediaProperties(j).QY =0.4;%M(n); %@func2; % [cm^-1]

    
end

