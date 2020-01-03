function [NodeBC,NodesOnElement,RegionOnElement,nummat,MateT,MatTypeTable,numBC,numel,ccell] = Meso_Locking_T(NodesOnElement,num_locked_g,NodeBC,locked_g,....
    meso_nen,nen_bulk,GrainVol,numel,numnpMicro,numnpMeso,MateT,nummat,MatTypeTable,RegionOnElement)
%Elina Geut
%Created 3/4/2019
%Last Modified 4/16/2019
%Reassigning BC with locked grains

N_BC_micro = zeros(3,0);
nen = size(NodesOnElement,2);
ccell = zeros(num_locked_g,1);
    
%Data initialization/Node,Material reassignment/Zeroing out micro scale
%and replacing by meso scale
    for i = 1:num_locked_g

        grain = locked_g(i);
        grainG1 = find(RegionOnElement==grain);
        numelemg = length(grainG1);
        micro = unique(reshape(NodesOnElement(grainG1,1:nen_bulk),1,numelemg*nen_bulk)); %Corrected
        lg = length(micro);
        N_BC_micro = [N_BC_micro [micro micro; ones(1,lg) 2*ones(1,lg); zeros(1,2*lg)]];

        %The MateT arrangement is [cR cL fR fL]

        %Value of 1st column 
        index1 = find(MateT(:,1) == grain);
        for k = 1:length(index1)
            flag_meso_top = 1;
            flag_micro_top = 0;

            MateT(index1(k),[5 7]) = [flag_meso_top flag_micro_top];
        end

        %Value of 2nd column
        index2 = find(MateT(:,2) == locked_g(i));
        for j = 1:length(index2)
            flag_meso_l = 1;
            flag_micro_l = 0;

            MateT(index2(j),[6 8]) = [flag_meso_l flag_micro_l];
        end

        % copy material properties
        l_MT = size(MateT,1)+1;
        GrainA = GrainVol(1,grain);
        MateT(l_MT,1:4) = [MateT(grain,1:3) GrainA];
        MatTypeTable(1:3,l_MT) = [l_MT 11 0]';

%             % Deal with grains on the boundary
%             BoundGrain = find(sum(MateT(:,1) == MateT(:,2),2) >= 1);
%             for count = 1:length(BoundGrain)
%                 DBCposition = BoundGrain(count);
%                 MateT(DBCposition,5:8) = [1 1 1 0];
%             end

        % Add coarse elements
        numel = numel + 1;
        NodesOnElement(numel,1:nen) = [numnpMicro+(grain-1)*meso_nen+1 numnpMicro+(grain-1)*meso_nen+2 ...
            numnpMicro+(grain-1)*meso_nen+3 zeros(1,nen-3)];
        RegionOnElement(numel) = l_MT;
        nummat = nummat + 1;
        ccell(i) = numel; %List of meso elements
    end

    NodeBC = [NodeBC; N_BC_micro'];

    %% Constrain the meso nodes of grains that are not locked
    cons_nodes = numnpMicro+1:numnpMicro+numnpMeso;

    for i = 1:num_locked_g
        grain = locked_g(i);
        cons_nodes(3*grain-2:3*grain) = 0;
    end
    cons_nodes = cons_nodes(cons_nodes>0);
    NodeBC2 = [cons_nodes'   ones(length(cons_nodes),1) zeros(length(cons_nodes),1)
               cons_nodes' 2*ones(length(cons_nodes),1) zeros(length(cons_nodes),1)];
    NodeBC = [NodeBC; NodeBC2];
        
   numBC = length(NodeBC); %Reassign number of boudary conditions
   
%    %4/25/2019
%    %Loop added to material property for CS flag 
%    if exist('FSon','var') && FSon == 1
%         for n = 1:num_locked_g
%            MateT(locked_g,9) = 2; %Flag = 2 to turn on both CS and FS
%         end
%    else
%        for n = 1:num_locked_g
%            MateT(locked_g,9) = 1; %Flag = 1 to turn on CS only
%        end
%    end