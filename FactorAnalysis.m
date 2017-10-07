X = [77 82 67 67 81; 63 78 80 70 81;
75 73 71 66 81; 55 72 63 70 68;
63 63 65 70 63; 53 61 72 64 73;
51 67 65 65 68; 59 70 68 62 56;
62 60 58 62 70; 64 72 60 62 45;
52 64 60 63 54; 55 67 59 62 44;
50 50 64 55 63; 65 63 58 56 37;
31 55 60 57 73; 60 64 56 54 40;
44 69 53 53 53; 42 69 61 55 45;
62 46 61 57 45; 31 49 62 63 62;
44 61 52 62 46; 49 41 61 49 64;
12 58 61 63 67; 49 53 49 62 47;
54 49 56 47 53; 54 53 46 59 44;
44 56 55 61 36; 18 44 50 57 81;
46 52 65 50 35; 32 45 49 57 64;
30 69 50 52 45; 46 49 53 59 37;
40 27 54 61 61; 31 42 48 54 68;
36 59 51 45 51; 56 40 56 54 35;
46 56 57 49 32; 45 42 55 56 40;
42 60 54 49 33; 40 63 53 54 25;
23 55 59 53 44; 48 48 49 51 37;
41 63 49 46 34; 46 52 53 41 40];
nums = size(X, 1);
dims = size(X, 2);
m = 2;
% normalize
Xbar = mean(X, 1);
X = X - Xbar;
S = X'*X/(nums - 1);
v = 1./sqrt(diag(S));
R = v.*S.*v';
nX = X.*v';
% principle factorization
% initial, 3rd method:
h2 = max(abs(R - diag(ones(dims,1))),[],2);
D = diag(1 - h2);
diff = inf;
epsilon = 1e-2;
while(diff > epsilon)
    R_star = R - D;
    [U, T] = eig(R_star);
    [~,idx] = sort(diag(T),'descend');
    U = U(:,idx(1:m));
    U = U./sqrt(dot(U,U,1));
    T = T(idx(1:m),idx(1:m));
    A = U*sqrt(T);
    h2 = dot(A,A,2);
    sigma2 = 1 - h2;
    diff = norm(sigma2 - diag(D));
    D = diag(sigma2);
end

% orthogonal rotation
A2 = A.*A;
mu = (A2(:,1) - A2(:,2)) ./ h2;
nu = (2 * A(:,1).*A(:,2)) ./ h2;
d = 2 * mu' * nu;
c = sum(mu.^2 - nu.^2);
alpha = sum(mu);
beta = sum(nu);
a = c -(alpha^2 - beta^2)/ dims;
b = d - 2 * alpha * beta / dims;
aphi = atan(b/a)/4;
%rotate = [cos(phi), -sin(phi); sin(phi), cos(phi)];
for i = -4:1:4
    tphi = aphi + i*pi/4;
    if(tphi > -pi/2 && tphi < pi/2 &&  sin(4*tphi)*b >= 0)
        phi = tphi;
        break;
    end
end
O = [cos(phi), -sin(phi); sin(phi), cos(phi)];
B = A*O;
h2 = dot(B,B,2);
q2 = dot(B,B,1);
F = B'*(R\nX');
display(B);
display(h2);
display(q2);
display(F);