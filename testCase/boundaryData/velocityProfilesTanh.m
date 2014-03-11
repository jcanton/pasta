clear
close all
clc
%
% Author: Jacopo Canton
% E-mail: jacopo.canton@mail.polimi.it
% Last revision: 28/8/2013
%
%------------------------------------------------------------------------------
%
% Simple script to save to file the values of an analytical velocity profile.
% The file will then be loaded by the continuation program which is capable
% of linearly interpolating these values to the actual nodes of the mesh.
%
%------------------------------------------------------------------------------
% PARAMETERS
%

tube = 2; % 1=inner, 2=outer

if ( tube == 1 )
   side =  1;
   xmin = -5.0;
   xmax = -5.0;
   ymin =  0;
   ymax =  0.5;
else
   side =  5;
   xmin = -5.0;
   xmax = -5.0;
   ymin =  0.6;
   ymax =  1;
end

npoints = 500;

%------------------------------------------------------------------------------
% DEFINITION OF THE VELOCITY PROFILES
%
if ( tube == 1 )

   uMean1 = 1;

   b = 5;

   uX = @(x,y)  tanh( b * ( 1 - y/ymax ) ) * uMean1;
   
   uY = @(x,y)  0;
   
   uZ = @(x,y)  0;

else

   uMean2 = 1;

   b = 5;

   ymed = ( ymin + ymax ) / 2;

   uX = @(x,y)  tanh( b * ( 1 - ( y - ymed )./( (ymin-ymax)/2*(y<ymed) + (ymax-ymin)/2*(y>=ymed) ) ) ) * uMean2;
   
   uY = @(x,y)  0;
   
   uZ = @(x,y)  0;

end

%------------------------------------------------------------------------------
% COMPUTATION OF THE POINT VALUES
%
xP = linspace(xmin,xmax,npoints)';
yP = linspace(ymin,ymax,npoints)';
zeroP = zeros(npoints,1);

uXp = uX(xP,yP);
uYp = uY(xP,yP);
uZp = uZ(xP,yP);

% if the user has specified a constant value, expand it
% in a constant vector of size 'npoints'
%
if ( size(uXp) == [1,1] )
   uXp = uXp * ones(npoints,1);
end
if ( size(uYp) == [1,1] )
   uYp = uYp * ones(npoints,1);
end
if ( size(uZp) == [1,1] )
   uZp = uZp * ones(npoints,1);
end

%------------------------------------------------------------------------------
% PLOT THE PROFILES
%
lin = 20;

% uX
figure(1)
hold on
grid on
plot(1:npoints,zeroP,'-k','lineWidth',2)
plot(1:npoints,uXp,'-b')
for i = 1:npoints/lin:npoints
   plot([i, i], [0, uXp(i)],'-b')
end
%plot([1 npoints],[.98*max(uXp) .98*max(uXp)],'-k')
%plot([1 npoints],[.99*max(uXp) .99*max(uXp)],'-k')
xlim([0 npoints+1])
xlabel('points')
ylabel('uX')
title('uX')

% uY
figure(2)
hold on
grid on
plot(1:npoints,zeroP,'-k','lineWidth',2)
plot(1:npoints,uYp,'-r')
for i = 1:npoints/lin:npoints
   plot([i, i], [0, uYp(i)],'-r')
end
xlim([0 npoints+1])
xlabel('points')
ylabel('uY')
title('uY')

% uZ
figure(3)
hold on
grid on
plot(1:npoints,zeroP,'-k','lineWidth',2)
plot(1:npoints,uZp,'-c')
for i = 1:npoints/lin:npoints
   plot([i, i], [0, uZp(i)],'-c')
end
xlim([0 npoints+1])
xlabel('points')
ylabel('uZ')
title('uZ')

%------------------------------------------------------------------------------
% PRINT TO FILE
%
fid=fopen(['side',num2str(side),'.dat'],'w');

fprintf(fid,'%d\n\n',npoints);
for i = 1:npoints
   fprintf(fid,'%.15f \t%.15f \t%.15f \t%.15f \t%.15f\n', uXp(i), uYp(i), uZp(i), xP(i), yP(i));
end

fclose(fid);
