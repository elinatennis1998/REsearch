%Elina Geut 
%6/27/2019
% Compute Displacement and Traction Jumps for MRDG case

GrainIntegFull

GrainAveD = GrainDisp./(ones(2,1)*GrainVol(1,:)); %computing the FS average to compare to CS
GrainAveE = GrainEps./(ones(3,1)*GrainVol(1,:));
GrainAveS = GrainSig(1:3,:)./(ones(3,1)*GrainVol(1,:));
GrainAveW = GrainRot./(ones(1,1)*GrainVol(1,:));  %Added Rotation
Ug_new = 0*Node_U_V;
% Ug_new = 0*Node_U_V;
for g = 1:numgrain
    % find elements in grain
    grainG1 = find(RegionOnElement==g);
    numelemg = length(grainG1);
    % find nodes in grain
    micro = unique(reshape(NodesOnElement(grainG1,1:nen_bulk),1,numelemg*nen_bulk)); %Corrected
    lg = length(micro);
    microX = Coordinates(micro,:);
    grainX = GrainXmid(:,g)';
    grainU = GrainAveD(:,g)';
    graine = GrainAveE(:,g)';
    grainw = GrainAveW(:,g)';
    grainE = [graine(1) graine(3); graine(3) graine(2)];
    grainW = [0 grainw(1); -grainw(1) 0];
    relat_X = microX - grainX;
    for i = 1:lg
        node = micro(i);
        Ug_new(node,1) = grainU(1) + grainE(1,:)*relat_X(i,:)' + grainW(1,:)*relat_X(i,:)';
        Ug_new(node,2) = grainU(2) + grainE(2,:)*relat_X(i,:)' + grainW(2,:)*relat_X(i,:)';
    end
end

if bCrys == 1 && nel == 4
    DispCoarse = Node_U_V - Ug_new;
    %Save the results to impose onto CS
    save('DispCoarse','DispCoarse');
else
    DispFine = Node_U_V - Ug_new;
    %Save the results to impose onto CS
    save('DispFine','DispFine');
end

