 %The Rayleigh-Ritz Method for finding potential energy
 
%  Fd1(isnan(Fd1))=0;
 ul(isnan(ul))=0;
 elem
 if elem <= numelCG
%      PE_int(elem,1) = PE_int(elem,1) + 0.5*ul(1,1:ndm*nel)*ElemK*ul(1,1:ndm*nel)';
%      PE_ext(elem,1) = PE_ext(elem,1) + ElemF'*ul(1,1:ndm*nel)';
% 
%      PE_int(elem,2) = PE_int(elem,2) + 0.5*ul(2,1:ndm*nel)*ElemK*ul(2,1:ndm*nel)';
%      PE_ext(elem,2) = PE_ext(elem,2) + ElemF'*ul(2,1:ndm*nel)';

 PE_int(elem,1) = PE_int(elem,1) + 0.5*Disp(1,:)*ElemK*Disp(1,:)'; 
 PE_ext(elem,1) = PE_ext(elem,1) + ElemF'*Disp(1,:)';
 
 PE_int(elem,2) = PE_int(elem,2) + 0.5*Disp(2,:)*ElemK*Disp(2,:)'; 
 PE_ext(elem,2) = PE_ext(elem,2) + ElemF'*Disp(2,:)';

    
 end