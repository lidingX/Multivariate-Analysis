im1 = imread('im1.JPEG');
im2 = imread('im2.JPEG');
im3 = imread('im3.JPEG');
im1 = rgb2gray(im1);
im2 = rgb2gray(im2);
im3 = rgb2gray(im3);
[n,m] = size(im1);
im1 = double(im1)/255;
im2  = double(im2)/255;
im3 = double(im3)/255;
show(im1,im2,im3);
%A = [0.7,0.2,0.1;0.3,0.4,0.3;0.2,0.15,0.65];
A = [0.19,0.22,0.199;0.57,0.45,0.4;0.3,0.3,0.44];
S = [im1(:)';im2(:)';im3(:)'];
X = A*S;
showmat(X,n,m);
[wX,mX] = whitening(X);
% B1 = llhsolver(wX, 'sub');
% Y1 = B1*(wX + mX);
% showmat(Y1,n,m);
B2 = batchFastsolver(wX);
%B2 = sequentialFastsolver(wX);
Y2 = B2*(wX + mX);
 showmat(Y2,n,m);
% B3 = batchFastsolver(wX);
% Y3 = B3*(wX + mX);
% showmat(Y3,n,m);
% ptn = 16;
% ptm = 16;
% num = 512;
% X1 = randpatch(im1,n,m, ptn, ptm, num);
% X2=  randpatch(im2,n,m, ptn, ptm, num);
% X3=  randpatch(im3,n,m, ptn, ptm, num);
% X4 = vertcat(X1,X2,X3);
% [wX4,mX4] = whitening(X4);
% B4 = sequentialFastsolver(wX4);
% %B4 = batchFastsolver(wX4);
% A4 = inv(B4);
% colnorms = dot(A4,A4);
% [~,idx] = sort(colnorms, 'descend');
% A4  = A4(:, idx);
% A4 = A4(:,1:3);
% Y4 = pinv(A4)*(wX4 + mX4);
% showmat(Y4,ptn,ptm);





function X = randpatch(im,n,m, ptn, ptm,num)
    rowmax = n - ptn + 1;
    colmax = m - ptm + 1;
    rowidx = randi([1,rowmax], num);
    colidx = randi([1,colmax],num);
    X = zeros(num,ptn*ptm);
    for i= 1: num
        sub =  im(rowidx(i): rowidx(i) + ptn - 1,  colidx(i):colidx(i) +ptm - 1);
        X(i,:) = sub(:)';
    end
end
function []=show(im1,im2,im3)
figure();
subplot(1,3,1);
imshow(im1);
subplot(1,3,2);
imshow(im2);
subplot(1,3,3);
imshow(im3);  
end
function [] = showmat(X,n,m)
    minX1 = min(X(1,:));
    maxX1 = max(X(1,:));
     minX2 = min(X(2,:));
    maxX2 = max(X(2,:));
     minX3 = min(X(3,:));
    maxX3 = max(X(3,:));
    X1 = (X(1,:) - minX1)./(maxX1 - minX1);
    X2 = (X(2,:) - minX2)./(maxX2 - minX2);
    X3 = (X(3,:) - minX3)./(maxX3 - minX3);
    show(reshape(X1,[n,m]),reshape(X2,[n,m]),reshape(X3,[n,m]));
end
function [wX, mX] = whitening(X)
mX = mean(X,2);
[~,n] =size(X);
X = X - mX;
SCov = X*X'/(n-1);
[E, D] = eig(SCov);
diagD = diag(D);
idx = diagD > 0;
E = E(:,idx);
Dinvs = diag(diagD(idx).^(-0.5));
T = E*Dinvs*E';
wX = T*X;
mX = T*mX;
end

function B = llhsolver(X, mode, rate)
if strcmp(mode,'super')
    logp = @(y) -2*log(cosh(y));
    g = @(y) -2*tanh(y);
    if (~exist('rate','var') || isempty(rate))
        rate = 0.005;
    end
end
if strcmp(mode,'sub')
     logp = @(y) -(y.^2)*0.5 + log(cosh(y));
    g = @(y) tanh(y) - y;
   if (~exist('rate','var') || isempty(rate))
    rate = 0.005;
   end
end
tol = 0.001;
[d,n] = size(X);
B = diag(ones(d,1));
y = B*X;
llh1 = log(abs(det(B))) + mean(sum(logp(y),1));
Delta = B + (g(y)*y'/n)*B; 
B = B + Delta * rate;
y = B*X;
llh2 = log(abs(det(B))) + mean(sum(logp(y),1));
round =0;
while(abs(llh2 - llh1) > tol * abs(llh1))
    %fprintf('%f %f %f\n',abs(llh2 - llh1),abs(llh1),  tol * abs(llh1));
    llh1 = llh2;
    y = B*X;
    Delta = B + (g(y)*y'/(n-1))*B; 
    B = B + Delta * rate;
    llh2 = log(abs(det(B))) + mean(sum(logp(y),1));
    round = round+1;
end
display(round);
end

function  B = sequentialFastsolver(X)
     g = @(y)  y.*exp(-(y.^2)*0.5);
     dg = @(y)  (1 - (y.^2)).*exp(-(y.^2)*0.5);
     tol = 0.01;
     [d,n] = size(X);  
     w = zeros(d,1);
     w(1) = 1;
     w = w/norm(w);
     Y = w'*X;
     %fprintf('fast\n');
     w_ =X*(g(Y))' / n -  mean(dg(Y))*w;
     w_ = w_/norm(w_);
     round = 0;
     while(abs(1 - abs(w_'*w)) > tol)
         w = w_;
         Y = w'*X;
         w_ =X*(g(Y))' / n -  mean(dg(Y))*w;
         w_ = w_/norm(w_);
         round = round + 1;
     end
    % fprintf('%d %d\n',1, round);
     W = zeros(d,d);
     W(:,1) = w_;
     p = 2;
     while( p <= d)
         w = zeros(d,1);
         w(p) = 1;
         w = w/norm(w);
         Y = w'*X;
         w_ =X*(g(Y))' / n -  mean(dg(Y))*w;
         tmpw_ = w_;
         for i = 1 : p-1
             tmpw_ = tmpw_ - w_'*W(:,i) *W(:,i);
         end
         w_ = tmpw_;
         w_ = w_/norm(w_);
         while(abs(1 - abs(w_'*w)) > tol)
             w = w_;
             Y = w'*X;
             w_ =X*(g(Y))' / n -  mean(dg(Y))*w;
             tmpw_ = w_;
             for i = 1 : p-1
                 tmpw_ = tmpw_ - w_'*W(:,i) *W(:,i);
             end
             w_ = tmpw_;
             w_ = w_/norm(w_);
             round = round + 1;
             %fprintf('%f\n',abs(w_'*w));
         end
          W(:,p) = w_;
         %fprintf('%d %d\n',p, round);
         p = p+1;
     end
     B = W';
end

function  B = batchFastsolver(X)
    tol = 0.01;
    [d,n] =size(X);
    dones  = ones(d,1);
    W = eye(d);
    W = W*(W'*W)^(-0.5);
    tmp = tanh(X'*W);
    c = mean(tmp.^2,1);
    W_ = X*tmp/n - (dones*(1-c)) .*W;
    W_ = W_*(W_'*W_)^(-0.5);
    round = 0;
    while(norm(W-W_) < tol)
        W = W_;
        tmp = tanh(X'*W);
        c = mean(tmp.^2,1);
        W_ = X*tmp/n - (dones*(1-c)) .*W;
        W_ = W_*(W_'*W_)^(-0.5);
        round = round  +  1;
    end
    display(round);
    B = W_'; 
end