
clear

% Mesh with nxn tiling
PSPS = 's'; %plane stress condition 
nen = 3;4; %number of nodes per element
nel = 3;4; %max number of nodes per element
numgh = 3;5;4;8;7;6; %number of grains along horiz. edge
numgs = 3;5;4;8;7;6; %number of grains along later. edge
numgrain = numgh*numgs; %total number of grains in RVE
bCrys = 2;1;3;8; %number of elem. along grain edge
CS_flag = 1; %Flag to collect data for post processing
if nen == 4
    tfact = 1; %quad 1
    btype = 0; %flag for meshing type (rect, triang, ...)
elseif nen == 3
    tfact = 2; %tri 2
    btype = 2;
end 
numelemg = bCrys^2; %number of elements in a grain
nu = numgh*bCrys; %number of elements along x
nv = numgs*bCrys; %number of elements along y
%% make the mesh of nodes and elements, centered at 0,0
Coordinates = [1 0 0
             2 5 0
             4 0 5
             3 5 5];
xlRVE = Coordinates;
GrainA = ((abs(Coordinates(1,2))+abs(Coordinates(4,2)))/numgs)^2; %For Square RVE
% GrainA = ((abs(Coordinates(1,2))+abs(Coordinates(4,2)))/numgs)*((abs(Coordinates(2,2))+abs(Coordinates(3,2)))/numgh);
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
m2 = [100 0.25 1]; %material type 1
m1 = [200 0.25 1]; %material type 2
[RegionOnElement,MatTypeTable,MateT,nummat,grainG] = GeomProp(numgrain,numelemg,tfact,numgs,numgh,bCrys,m2,m1,nel,nu);% MateT = [m1; m2; m1; m2; m2; m1; m2; m1; m1; m2; m1; m2; m2; m1; m2; m1];

%% provide flags for desired outputs
OptFlag = [0 1 1 0 0 1 1]';
SEHist = 1;
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
wRVE = 0;
GrainIntegV
InterIntegV

%% Key Part: set up fluxes and jumps on the grain boundaries and RVE
TFS = 2; % FE 3; % Sachs 1; %Taylor 5; % Fine Scale  4; % Combo 
factorTS = 0;1;.5;.25;.75; % factor for Taylor vs Sachs; 1 means Taylor
factorT = factorTS;1/3;
factorS = 1-factorTS;1/3;
if TFS == 1
% factorTS = 1;
TaylorSet
elseif TFS == 3
% factorTS = 0;
SachsSet
elseif TFS == 4
TaylorSet
SachsSet
elseif TFS == 5
    %     load('CoordinatesIT.mat');
    %     CoordinatesI=CoordinatesIT;
    load('DispFine','DispFine');
    load('GrainEps','GrainEps');
    load('GrainDisp','GrainDisp');
    load('GrainVol','GrainVol');
end

FormRVEDirichlet

%% Zero out meso grains; for this file as it is, I don't want to compute the
% the sequence has to be retained: from low to high order of grain number
    for j = 1:numgrain
        locked_g(1,j) = j;
    end
% locked_g = [];
% locked_g = [5];
% % locked_g = [1 2 3 4 5 6 10 11 15 16 20 21 22 23 24 25];
num_locked_g = length(locked_g); 
meso_nen = 3; %Number of nodes per mose element
sub = 11; %subroutine number for a meso element
combo = 1; %Flag for plotting FS on CS 
[NodeBC,NodesOnElement,RegionOnElement,nummat,MateT,MatTypeTable,numBC,numel,ccell] = Meso_Locking_T(NodesOnElement,num_locked_g,NodeBC,locked_g,....
    meso_nen,nen_bulk,GrainVol,numel,numnpMicro,numnpMeso,MateT,nummat,MatTypeTable,RegionOnElement);
ProbType = [numnp numel nummat 2 2 nen];