function [p,t,v,l] = carpet1216(level,numrefinements)
%sets up and solves eigenvalue problem on outer approximation of 12/16
%carpet

m = level;
%calculate side length
l = 1/(4^m)

%build lower leftmost square
x1 = 0;
x2 = x1;
x3 = x1+l;
x4 = x1+l;
y1 = 0;
y2 = y1+l;
y3 = y1+l;
y4 = y1;

A = [2;4;x1;x2;x3;x4;y1;y2;y3;y4];

%construct the mth level using a recursive type algorigthm where 
%each square is a column in a matrix and each iteration builds
%twelve squares for each square already in the matrix
for i=1:m
    for j=1:12^(i-1)
        x1 = A(3,j)+1/(4^i);
        y1 = A(7,j);
        x2 = x1;
        x3 = x1+l;
        x4 = x1+l;
        y2 = y1+l;
        y3 = y1+l;
        y4 = y1;
        
        B = [2;4;x1;x2;x3;x4;y1;y2;y3;y4];
        A = [A B];
        
        x1 = A(3,j);
        y1 = A(7,j)+1/(4^i);
        x2 = x1;
        x3 = x1+l;
        x4 = x1+l;
        y2 = y1+l;
        y3 = y1+l;
        y4 = y1;
        
        B = [2;4;x1;x2;x3;x4;y1;y2;y3;y4];
        A = [A B];
        
        x1 = A(3,j)+2/(4^i);
        y1 = A(7,j);
        x2 = x1;
        x3 = x1+l;
        x4 = x1+l;
        y2 = y1+l;
        y3 = y1+l;
        y4 = y1;
        
        B = [2;4;x1;x2;x3;x4;y1;y2;y3;y4];
        A = [A B];
        
        x1 = A(3,j);
        y1 = A(7,j)+2/(4^i);
        x2 = x1;
        x3 = x1+l;
        x4 = x1+l;
        y2 = y1+l;
        y3 = y1+l;
        y4 = y1;
        
        B = [2;4;x1;x2;x3;x4;y1;y2;y3;y4];
        A = [A B];
        
        x1 = A(3,j);
        y1 = A(7,j)+3/(4^i);
        x2 = x1;
        x3 = x1+l;
        x4 = x1+l;
        y2 = y1+l;
        y3 = y1+l;
        y4 = y1;
        
        B = [2;4;x1;x2;x3;x4;y1;y2;y3;y4];
        A = [A B];
        
        x1 = A(3,j)+1/(4^i);
        y1 = A(7,j)+3/(4^i);
        x2 = x1;
        x3 = x1+l;
        x4 = x1+l;
        y2 = y1+l;
        y3 = y1+l;
        y4 = y1;
        
        B = [2;4;x1;x2;x3;x4;y1;y2;y3;y4];
        A = [A B];
        
        x1 = A(3,j)+2/(4^i);
        y1 = A(7,j)+3/(4^i);
        x2 = x1;
        x3 = x1+l;
        x4 = x1+l;
        y2 = y1+l;
        y3 = y1+l;
        y4 = y1;
        
        B = [2;4;x1;x2;x3;x4;y1;y2;y3;y4];
        A = [A B];
        
        x1 = A(3,j)+3/(4^i);
        y1 = A(7,j)+3/(4^i);
        x2 = x1;
        x3 = x1+l;
        x4 = x1+l;
        y2 = y1+l;
        y3 = y1+l;
        y4 = y1;
        
        B = [2;4;x1;x2;x3;x4;y1;y2;y3;y4];
        A = [A B];
        
        x1 = A(3,j)+3/(4^i);
        y1 = A(7,j)+2/(4^i);
        x2 = x1;
        x3 = x1+l;
        x4 = x1+l;
        y2 = y1+l;
        y3 = y1+l;
        y4 = y1;
        
        B = [2;4;x1;x2;x3;x4;y1;y2;y3;y4];
        A = [A B];
        
        x1 = A(3,j)+3/(4^i);
        y1 = A(7,j)+1/(4^i);
        x2 = x1;
        x3 = x1+l;
        x4 = x1+l;
        y2 = y1+l;
        y3 = y1+l;
        y4 = y1;
        
        B = [2;4;x1;x2;x3;x4;y1;y2;y3;y4];
        A = [A B];
        
        x1 = A(3,j)+3/(4^i);
        y1 = A(7,j);
        x2 = x1;
        x3 = x1+l;
        x4 = x1+l;
        y2 = y1+l;
        y3 = y1+l;
        y4 = y1;
        
        B = [2;4;x1;x2;x3;x4;y1;y2;y3;y4];
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

pdemesh(p,e,t);                        %display mesh

%solve eigenvalue problem, only looking for eigenvals between 0 and 1000
%[v,l]=pdeeig(B,p,e,t,1,0,1,[0,1000]);

%pdesurf(p,t,v(:,6))                    %display first eigenmode

