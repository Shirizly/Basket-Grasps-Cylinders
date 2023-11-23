function [N,L,W] = TSNormals(com,p1,p2,p3,nsol)
% Computes force magnitudes to cancel gravity on the object, will not be
% accurate in non-equilibrium stances.
p1t = p1-com;
p2t = p2-com;
n1 = -(p1t-(p1t.'*nsol)*nsol);
n2 = -(p2t-(p2t.'*nsol)*nsol);
n3 = nsol;
N = [n1,n2,n3];
w1 = [n1;p1];
w2 = [n2;p2];
w3 = [n3;p3];
W = [w1,w2,w3];
g = [0;0;1];
L = [n1 n2 n3]\g;
end

