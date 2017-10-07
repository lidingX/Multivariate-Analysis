txt1 = ['1', '2', '3', '4', '5', '6'];
txt2 = ['A','B','C','D','E' ,'F','G','H','I'];
txt3 =['1','2','3','4','5','6'];
txt4 = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'];
K = [58 19 1 33 110 163;
5 1 1 1 31 25;
20 3 2 10 78 121;
1 0 0 0 12 14;
1 0 0 2 28 3;
1 0 1 4 16 10;
6 2 1 7 19 17;
3 0 0 10 19 12;
3 0 1 0 30 31];
[n,p] = size(K);
sumK = sum(sum(K));
F = K/ sumK;
Dp = diag(sum(F,1));
Dn = diag(sum(F,2));
S_star = (Dp.^0.5\F')*((Dn\F) / (Dp.^0.5));
[U, T] = eig(S_star);
[~,idx] = sort(diag(T),'descend');
U = U(:,idx);
Phi = (Dp.^0.5)\U;
lambda =diag(T);
lambda = lambda(idx);
sqlambda2 = sqrt(lambda(2));
sqlambda3= sqrt(lambda(3));
phi2 = sqlambda2 * Phi(:,2);
phi3 = sqlambda3 * Phi(:,3);
psi2 = Dn\F*Phi(:,2);
psi3 = Dn\F*Phi(:,3);
x = [phi2; psi2];
y = [phi3; psi3];
c = [repmat([1 0 0],size(phi2,1),1) ; repmat([0 0 0],size(psi2,1),1)];
figure();
dx = phi2(2)*0.03;
dy = phi3(3)*0.03;
for k=1:p
text(x(k)+dx,y(k)+dy,txt1(k));
end
for k=p + 1: n+p
text(x(k)+ dx,y(k) + dy,txt2(k - p));
end
hold on;
scatter(x,y,[],c,'filled');

F = F - sum(F,2)*sum(F,1);
S_star = (Dp.^0.5\F')*((Dn\F) / (Dp.^0.5));
[U, T] = eig(S_star);
[~,idx] = sort(diag(T),'descend');
U = U(:,idx);
Phi = (Dp.^0.5)\U;
lambda =diag(T);
lambda = lambda(idx);
sqlambda1 = sqrt(lambda(1));
sqlambda2= sqrt(lambda(2));
phi1 = sqlambda1 * Phi(:,1);
phi2 = sqlambda2 * Phi(:,2);
psi1 = Dn\F*Phi(:,1);
psi2 = Dn\F*Phi(:,2);
x = [phi1; psi1];
y = [phi2; psi2];
c = [repmat([1 0 0],size(phi1,1),1) ; repmat([0 0 0],size(psi1,1),1)];
figure();
dx = phi2(2)*0.04;
dy = phi3(3)*0.04;
for k=1:p
text(x(k)+dx,y(k)+dy,txt1(k));
end
for k=p + 1: n+p
text(x(k)+ dx,y(k) + dy,txt2(k - p));
end
hold on;
scatter(x,y,[],c,'filled');

K =[190.33 43.77 9.73 60.54 49.01 9.04;
135.20 36.40 10.47 44.16 36.49 3.94;
95.21 22.83 9.30 22.44 22.81 2.80;
104.78 25.11 6.40 9.89 18.17 3.25;
128.41 27.63 8.94 12.58 23.99 3.27;
145.68 32.83 17.79 27.29 39.09 3.47;
159.37 33.38 18.37 11.81 25.29 5.22;
116.22 29.57 13.24 13.76 21.75 6.04;
221.11 38.64 12.53 115.65 50.82 5.89;
144.98 29.12 11.67 42.60 27.30 5.74;
169.92 32.75 12.72 47.12 34.35 5.00;
153.11 23.09 15.62 23.54 18.18 6.39;
144.92 21.26 16.96 19.52 21.75 6.73;
140.54 21.50 17.64 19.19 15.97 4.94;
115.84 30.26 12.20 33.61 33.77 3.85;
101.18 23.26 8.46 20.20 20.50 4.30];
[n,p] = size(K);
sumK = sum(sum(K));
F = K/ sumK;
Dp = diag(sum(F,1));
Dn = diag(sum(F,2));
F = F - sum(F,2)*sum(F,1);
S_star = (Dp.^0.5\F')*((Dn\F) / (Dp.^0.5));
[U, T] = eig(S_star);
[~,idx] = sort(diag(T),'descend');
U = U(:,idx);
Phi = (Dp.^0.5)\U;
lambda =diag(T);
lambda = lambda(idx);
sqlambda1 = sqrt(lambda(1));
sqlambda2= sqrt(lambda(2));
phi1 = sqlambda1 * Phi(:,1);
phi2 = sqlambda2 * Phi(:,2);
psi1 = Dn\F*Phi(:,1);
psi2 = Dn\F*Phi(:,2);
x = [phi1; psi1];
y = [phi2; psi2];
c = [repmat([1 0 0],size(phi1,1),1) ; repmat([0 0 0],size(psi1,1),1)];
figure();
dx = phi2(2)*0.03;
dy = phi3(3)*0.05;
for k=1:p
text(x(k)+dx,y(k)+dy,txt3(k));
end
for k=p + 1: n+p
text(x(k)+ dx,y(k) + dy,txt4(k - p));
end
hold on;
scatter(x,y,[],c,'filled');


