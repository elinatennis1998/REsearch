%# data
% t = linspace(0,8*pi,200);
% x = 20*t; y = cos(t); z = sin(t);
% 
% %# plot 3D line
% plot3(x,y,z)
% axis tight, grid on, view(35,40)
% c = 1:numel(t);      %# colors
% h = surface([x(:), x(:)], [y(:), y(:)], [z(:), z(:)], ...
%     [c(:), c(:)], 'EdgeColor','flat', 'FaceColor','none');
% colormap( jet(numel(t)))
% 
% x = 0:.05:2*pi;
% y = sin(x);
% z = zeros(size(x));
% col = x;  % This is the color, vary with x in this case.
% 
% surface([x;x],[y;y],[z;z],[col;col],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2)
    
    x1 = Coordinates(NodesOnElement(numelCG+1:numel-2,1),1)';
    x2 = Coordinates(NodesOnElement(numelCG+1:numel-2,2),1)';
    y1 = Coordinates(NodesOnElement(numelCG+1:numel-2,1),2)';
    y2 = Coordinates(NodesOnElement(numelCG+1:numel-2,2),2)';
    z1 = Z;
    z2 = Z;
    figure(1) % Fint_1
    hold on
   surface([x1;x2],[y1;y2],[z1;z2],[Fint_1;Fint_1],'facecol','no','edgecol','interp','linew',2);
    hold off