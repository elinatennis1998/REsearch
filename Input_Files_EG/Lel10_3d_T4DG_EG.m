%Elina Geut
%Input file for 3D MRDG case
%Created 8/15/2019


clear

nen = 4;
nel = 4;

xl = [1 0 0 0
      2 2 0 0
      4 0 1 0
      3 2 1 0
      5 0 0 3
      6 2 0 3
      8 0 1 3
      7 2 1 3];
  
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block3d('cart',6,6,6,1,1,1,11,xl,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';
%% Set up material properties 
nummat = 3;
RegionOnElement(1:6) = 3;
RegionOnElement(13:18) = 2;
RegionOnElement(31:36) = 2;
MatTypeTable = [1 2 3
                1 1 1];
MateT = [1 1 1]'*[190e3 0.3 1];

%% Output flags
DHist = 1;
FHist = 1;
SHist = 1;
SEHist = 1;

%% separate the grains by duplicating nodes
numnpCG = numnp;
InterTypes = [0 0 0
              1 0 0
              1 1 0]; % only put CZM between the element edges between materials 1-2
DEIProgram3

%% Insert DG couplers along interfaces and RVE boundary, add a node for each
% grain, make some special coordinate lists
ndm = 3;
ielI = 10;
ielB = 10;
matepropI = [0 0 1 1];
matepropB = [0 1 1 0];
InterDGallG

%% RVE/macroscale BC
uRVE = [0 0 0];
eRVE = [0.2 0 0 0 0 0];
wRVE = [0 0 0];
GrainIntegV
InterIntegV

%% Key Part: set up fluxes and jumps on the grain boundaries and RVE
TFS = 2;5; % FE 5; % Fine Scale 1; %Taylor 3; % Sachs  4; % Combo 
factorTS = 1;0;.5;.25;.75; % factor for Taylor vs Sachs; 1 means Taylor
factorT = factorTS;1/3;
factorS = 1-factorTS;1/3;
FormRVEDirichlet