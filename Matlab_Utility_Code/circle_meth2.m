function H=circle_meth2(center,radius,N,style)
%---------------------------------------------------------------------------------------------
% H=CIRCLE(CENTER,RADIUS,N,STYLE)
% This routine draws a circle with center defined as
% a vector CENTER, radius as a scaler RADIS. NOP is 
% the number of points on the circle. As to STYLE,
% use it the same way as you use the rountine PLOT.
% Since the handle of the object is returned, you
% use routine SET to get the best result.
%
%   Usage Examples,
%
%   circle([1,3],3,1000,':'); 
%   circle([2,4],2,1000,'--');
%%---------------------------------------------------------------------------------------------
if (nargin <3)
 N=256;
 %if (nargin <3),
 % error('Please see help for INPUT DATA.');
elseif (nargin==3)
    style='b-';
end;

t=(0:N)*2*pi/N;
plot( radius*cos(t)+center(1), radius*sin(t)+center(2),style);
axis square

%THETA=linspace(0,2*pi,NOP);
%RHO=ones(1,NOP)*radius;
%[X,Y] = pol2cart(THETA,RHO);
%X=X+center(1);
%Y=Y+center(2);
%H=plot(X,Y,style);
%axis square;
