function [p,t,v,l] = SG(level,initialoverlap,overlaprate,numrefinements)
%solve the eigenvalue problem on SG

m = level;
o = initialoverlap;
r = overlaprate/2;

t = sqrt(3);

%calculate side length
l = (1+o*((1-.5^m)/.5))/(2^(m))

%build bottom left triangle
x1 = 0;
x2 = x1+l;
x3 = x1+l/2;
y1 = 0;
y2 = 0;
y3 = y1+t*l/2;

A = [2;3;x1;x2;x3;y1;y2;y3];

%construct the mth level using a recursive type algorigthm where 
%each triangle is a column in a matrix and each iteration builds
%two triangles for each triangle already in the matrix (m iterations)
for i=1:m
    for j=1:3^(i-1)
        x1 = A(3,j)+2^(-i)-o*(r^(i-1))/2;
        y1 = A(6,j);
        x2 = x1+l;
        x3 = x1+l/2;
        y2 = y1;
        y3 = y1+t*l/2;
        
        B = [2;3;x1;x2;x3;y1;y2;y3];
        A = [A B];
        
        x1 = A(3,j)+((2^(-i))-o*(r^(i-1))/2)/2;
        y1 = A(6,j)+t*((2^(-i))-o*(r^(i-1))/2)/2;
        x2 = x1+l;
        x3 = x1+l/2;
        y2 = y1;
        y3 = y1+t*l/2;
        
        B = [2;3;x1;x2;x3;y1;y2;y3];
        A = [A B];
    end  
end

%tell pdetoolbox to build the 'decomposed geometry' and triangulate
g=decsg(A);
[p,e,t]=initmesh(g);

%describe the boundary conditions (Neumann in this case)
%the boundary conditions are described by a matrix with
%one column for each column in the decomposed geometry
B = [1;0;1;1;48;48;48;48;49;48];
C = [1;0;1;1;48;48;48;48;49;48];
for j=1:length(e)
    B = [B C];
end

%perform any refinements to the triangulation
for k=1:numrefinements
    [p,e,t]=refinemesh(g,p,e,t);
end

%pdemesh(p,e,t);                        %display mesh

%solve eigenvalue problem, only looking for eigenvals between 0 and 2000
[v,l]=pdeeig(B,p,e,t,1,0,1,[0,2000]);
l
%pdesurf(p,t,v(:,1))                    %display first eigenmode

