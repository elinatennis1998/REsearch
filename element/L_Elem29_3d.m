% 02/11/2017
%
% -Linear MPC element using Lagrange multipliers
% -In isw=1 case, the element turns off all dofs that exceed ndm because
%  they are not used for pure-displacement.
%     Purpose: Nodal constraint element using Lagrange multipliers:
%              u_1 - u_2 - a*u_3 - b*u_4 - c*u_5 = 0 for i = 1,...,ndf
% Assumes MPC elements are last in the mesh, and numMPC says how many there
% are.

switch isw %Task Switch
    
    case 1
        
        if ndf > ndm
            
            for i = ndm+1:ndf
                lie(i,1) = 0;
            end
            
        end

    case {3,21}
         
        pbcID = elem - (numel-numMP);
        Coeffs = sign(MPListC(pbcID,1:nel-1)); %only care about direction, not magnitude
        ElemK = zeros(nst);
%         ElemF = zeros(nst,1);
        one = 1.0;

        for i = 1:ndf

            for j = 1:nel-1
              ElemK(nel*ndf-ndf+i,(j-1)*ndf+i) = Coeffs(j);
              ElemK((j-1)*ndf+i,nel*ndf-ndf+i) = Coeffs(j);

            end
            
        end
        
        ElemF = -ElemK*reshape(ul,nst,1);
ElemK;
        
end %Task Switch
