function [p,t,v,l] = sawtooth(level,h,numrefinements)

m = level;
%using preset overlap (o) which shrinks with level
o = (2^(-m))/10;
r = 1;

t = sqrt(3);

%compute the base length
l = (1-o)/(2^m)+o

%build the leftmost triangle
x1 = 0;
x2 = x1+l;
x3 = x1+l/2;
y1 = 0;
y2 = 0;
y3 = y1+h;

A = [2;3;x1;x2;x3;y1;y2;y3];

%build all the other 2^m-1 triangles for level m
for j=1:2^(m)-1
    x1 = j*(l-o);
    y1 = 0;
    x2 = x1+l;
    x3 = x1+l/2;
    y2 = y1;
    y3 = y1+h;

    B = [2;3;x1;x2;x3;y1;y2;y3];
    A = [A B];

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

%pdesurf(p,t,v(:,1))                    %display first eigenmode

