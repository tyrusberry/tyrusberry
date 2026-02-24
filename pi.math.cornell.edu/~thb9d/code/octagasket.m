function [p,t,v,l] = octagasketEigenvalues(level,numrefinements)

m = level;
t = sqrt(2);

%determine the initial position (entirely determined by geometry)
x = 1/(1+2/t);

%compute the side length
l = x/((2+2/t)^(m-1));

%build the left most octagon of the bottom row
x1 = x/t;
y1 = 0;
x2 = x1+l;
y2 = 0;
x3 = x2+l/t;
y3 = l/t;
x4 = x3;
y4 = y3+l;
x5 = x2;
y5 = y4+l/t;
x6 = x1;
y6 = y5;
x7 = x1-l/t;
y7 = y4;
x8 = x7;
y8 = y3;

%note that the last bit of data indicates where this octagon is
%on the smallest 'octagon' containing it (1 being leftmost bottom row 
%and proceeding counter clockwise).  It is not part of the geometry and
%is only used by my code.
A = [2;8;x1;x2;x3;x4;x5;x6;x7;x8;y1;y2;y3;y4;y5;y6;y7;y8;1];

%construct the mth level using a recursive type algorigthm where 
%each octagon is a column in a matrix and each iteration builds
%seven octagons for each octagon already in the matrix.
%The build pattern has to go around the octagons in a confusing manner
%so there is some tricky math here.
for i=1:m-1
    for j=1:8^(i-1)
        n = A(19,j);
        z = x/((2+2/t)^(i-1));
        a = 0;
        b = 0;  
        for k=1:7
            q = n+k;
            if (q > 8)
                q = q-8;
            end

            a = a + cos((q-2)*pi/4);
            b = b + sin((q-2)*pi/4);
            c=0;
            d=0;
            for s=2:q
            c = c + cos((s-2)*pi/4);
            d = d + sin((s-2)*pi/4);
            end
            x1 = A(2+n,j)+a*z;
            y = A(10+n,j)+b*z;
            
            x1 = x1-c*l;
            y = y-d*l;

            x1 = x1;
            y1 = y;
            x2 = x1+l;
            y2 = y1;
            x3 = x2+l/t;
            y3 = y1+l/t;
            x4 = x3;
            y4 = y3+l;
            x5 = x2;
            y5 = y4+l/t;
            x6 = x1;
            y6 = y5;
            x7 = x1-l/t;
            y7 = y4;
            x8 = x7;
            y8 = y3;


            B = [2;8;x1;x2;x3;x4;x5;x6;x7;x8;y1;y2;y3;y4;y5;y6;y7;y8;q];
            
            A = [A B];
        end
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

%pdesurf(p,t,v(:,5))                    %display first eigenmode

