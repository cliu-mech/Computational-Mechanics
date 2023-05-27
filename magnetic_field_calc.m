%% Magnetic Field Computation
%%
clc;clear;close all;
%%
scalarFactor=0.1;
upperLimit=100;
nSegments=100*upperLimit;
z=linspace(0,upperLimit,nSegments+1);
t=z*pi;
r=1;
x=cos(t);
y=sin(t);

%%
mu0=pi*4e-7;% the permeability of free space [T*m/A]
I=1;% current [A]

wireSegments=[x;y;(z-0.5*max(z))*scalarFactor];% wire segments

upperLimit=9;
d=0.6;
xRange=[-1 1]*upperLimit;
yRange=[-1 1]*upperLimit;
zRange=[-1 1]*upperLimit;

x0=xRange(1):d:xRange(2);
y0=yRange(1):d:yRange(2);
z0=zRange(1):d:zRange(2);

%[X,Y,Z]=meshgrid(x0,y0,z0);
%pts=zeros(3,length(x0)*length(y0)*length(z0));

%%
[X,Y,Z,U,V,W]=b_segments(x0,y0,z0,wireSegments,mu0,I);


%% Display results

figure;
plot3(x,y,(z-0.5*max(z))*scalarFactor,'LineWidth',2);
hold on;
quiver3(X,Y,Z,U,V,W);
%plot3(x,y,z*scalarFactor,'LineWidth',2);
xlim(xRange);
ylim(yRange);
zlim(zRange);
grid minor;


%% local function
%
function [X,Y,Z,U,V,W]=b_segments(x0,y0,z0,segments,mu,I)

[X,Y,Z]=meshgrid(x0,y0,z0);

U=zeros(size(X));
V=zeros(size(Y));
W=zeros(size(Z));

idxPt=1;
for idxZ=1:length(z0)
    for idxY=1:length(y0)
        for idxX=1:length(x0)
            pt=[X(idxX,idxY,idxZ);Y(idxX,idxY,idxZ);Z(idxX,idxY,idxZ)];
            bi=ui_segments(pt,segments,mu*I);
            U(idxX,idxY,idxZ)=bi(1);
            V(idxX,idxY,idxZ)=bi(2);
            W(idxX,idxY,idxZ)=bi(3);
            idxPt=idxPt+1;
        end
    end
end

end


% Function calculating magnetic vector at given points
function ui=ui_segments(pts,segments,gamma)

nPoints=size(pts,2);
nSegments=size(segments,2)-1;
ui=zeros(3,nPoints);

for iPoint=1:nPoints

    for iSegment=1:nSegments
        rs=pts(:,iPoint)-segments(:,iSegment);
        re=pts(:,iPoint)-segments(:,iSegment+1);
        %cross(rs,re)
        if norm(cross(rs,re))~=0
            h=norm(cross(rs,re))/norm(rs-re);
            cosineBs=(rs-re)'*rs/norm(rs-re)/norm(rs);
            cosineBe=(rs-re)'*re/norm(rs-re)/norm(re);
            uiNorm=gamma*(cosineBs-cosineBe)/4/pi/h;
        %pause
            tempCross=cross(rs,re);
            tempCross=tempCross/norm(tempCross);
    
            ui(:,iPoint)=ui(:,iPoint)+uiNorm*tempCross;
        end
    end
end

end