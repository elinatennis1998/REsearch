% CST solid element with 4x4 grain RVE. Coarse
% scale implemented as a material parameter

clear

% Mesh with nxn tiling
PSPS = 's'; %plane stress condition 
nen = 3;4; %number of nodes per element
nel = 3;4; %max number of nodes per element
numgh = 3;5;4;6;8;7; %number of grains along horiz. edge
numgs = 3;5;4;6;8;7; %number of grains along later. edge
numgrain = 2; %total number of grains in RVE
bCrys = 6; %number of elem. along grain edge
CS_flag = 0; %Flag to collect data for post processing
MRDGF = 2;1;0; %1 - include only interface; 2 - include both
IBP = 0;
if nen == 4
    tfact = 1; %quad 1
    btype = 0; %flag for meshing type (rect, triang, ...)
elseif nen == 3
    tfact = 2; %tri 2
    btype = 6;
end 
numelemg = bCrys^2; %number of elements in a grain
nu = numgh*bCrys; %number of elements along x
nv = numgs*bCrys; %number of elements along y
%% make the mesh of nodes and elements, centered at 0,0
Coordinates = [1 -3 -3
             2 3 -3
             4 -3 3
             3 3 3];
% Coordinates = [1 0 0
%                  2 6 0
%                  4 0 6
%                  3 6 6];
node1 = 1;
elmt1 = 1;
mat = 1;
rskip = 0;
type = 'cart';
[Coordinates,NodesOnElement,RegionOnElement,numnp,numel] = block2d(type,nu,nv,node1,elmt1,mat,rskip,btype,Coordinates,nen);
Coordinates = Coordinates';
NodesOnElement = NodesOnElement';
%% Map element ID onto grain ID
%Assign individual region and material property to each grain 
%Note: alterphase is a built-in function for maerial assignment, id desired
%otherwise, must hard code
m2 = [100 0.25 1];[100 0 1]; %material type 1
m1 = [200 0.25 1];[200 0 1]; %material type 2
MatTypeTable = [1 2; 1 1; 0 0];
MateT = [m2; m1];
nummat = 2;
% grainG = zeros(2,numel);
for j = 1:bCrys
RegionOnElement(2*bCrys^2*numgh+2*bCrys+1+2*bCrys*numgh*(j-1):2*bCrys^2*numgh+2*bCrys+2*bCrys+2*bCrys*numgh*(j-1)) = 2;
end
% nn = 6; %for numg = 4
% n1 = nn+2;
% for j = 1:nn 
% RegionOnElement(2*bCrys^2*numgh+2*bCrys+1+2*bCrys*numgh*(j-1):2*bCrys^2*numgh+2*bCrys+2*bCrys+n1+2*bCrys*numgh*(j-1)) = 2;
% end
% grainG(1,1:numel-bCrys^2) = find(RegionOnElement==1);
% grainG(2,1:bCrys^2) = find(RegionOnElement==2);
% MateT = [m1; m2; m1; m2; m2; m1; m2; m1; m1; m2; m1; m2; m2; m1; m2; m1];

%% provide flags for desired outputs
OptFlag = [0 1 1 0 0 1 1]';
SEHist = 1;
FHist = 1;
%% separate the grains by duplicating nodes
numnpCG = numnp;
InterTypes = tril(ones(nummat),-1); % only put CZM between the element edges between materials 1-2
DEIProgram2

%% Insert DG couplers along interfaces and RVE boundary, add a node for each
% grain, make some special coordinate lists
ndm = 2;
ielI = 10;
ielB = 12;
matepropI = [0 0 1 1];
matepropB = [0 1 1 0];
InterDGallG2

%% RVE/macroscale BC
uRVE = [0 0];
eRVE = [0.2 0 0];
% eRVE = [0 0.2 0];
% eRVE = [0 0 0.2];
wRVE = 0;
GrainIntegV
InterIntegV
%% Key Part: set up fluxes and jumps on the grain boundaries and RVE
TFS = 1; %Taylor 5; % FE 5; % Fine Scale 3; % Sachs  4; % Combo 
factorTS = 1;0;.5;.25;.75; % factor for Taylor vs Sachs; 1 means Taylor
factorT = factorTS;1/3;
factorS = 1-factorTS;1/3;
FormRVEDirichlet
for j = 1:numgrain
    locked_g(1,j) = j;
end
num_locked_g = length(locked_g);
meso_nen = 3; %Number of nodes per mose element
sub = 11; %subroutine number for a meso element
combo = 1; %Flag for plotting FS on CS
[NodeBC,NodesOnElement,RegionOnElement,nummat,MateT,MatTypeTable,numBC,numel,ccell] = Meso_Locking_T(NodesOnElement,num_locked_g,NodeBC,locked_g,....
    meso_nen,nen_bulk,GrainVol,numel,numnpMicro,numnpMeso,MateT,nummat,MatTypeTable,RegionOnElement);

ProbType = [numnp numel nummat 2 2 nen];