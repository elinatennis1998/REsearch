%Elina Geut 
%Script for imposing the FS stresses onto the CS
%Created 5/10/2019
%Last modified 5/10/2019

count = numelCG;
StreListCF = StreListE;

for j = 1:num_locked_g
    grain_e = find(RegionOnElement (:,1) == locked_g(j));
    for i = 1:count
        StreListCF(1,grain_e) = StreListE(1,ccell(j));
        StreListCF(2,grain_e) = StreListE(2,ccell(j));
        StreListCF(3,grain_e) = StreListE(3,ccell(j));
    end
end
 %plotElemCont2(Coordinates,StreListCF(1,1:numelCG),NodesOnElement(1:numelCG,1:4),1,(1:size(NodesOnElement,1)),[1 0 0])