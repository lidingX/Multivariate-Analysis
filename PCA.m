data = load('umist_cropped.mat');
facedat = data.facedat;
dirnames = data.dirnames;
face = facedat{1};
face = double(face);
X = reshape(face,112*92,38);
X = X';
n = size(X,1);
mX = mean(X,1); 
Xc = X - mX;
S = Xc*Xc'/(n-1);
[V,D] = eig(S);
[~,idx] = sort(diag(D),'descend');
V = V(:,idx(1:6));
U = Xc'*V;
U = U./sqrt(dot(U,U,1));
mX = reshape(mX,112,92);
imagesc(mX);
colormap('gray');
figure;
for i=1:6
    subplot(2,3,i);
    imagesc(reshape(U(:,i),112,92));
    title(['Eigenfaces ' num2str(i)]); 
    colormap('gray');
end
var = sqrt(sum(Xc.*Xc,1)/(n-1));
Xc = Xc./var;
S = Xc*Xc'/(n-1);
[V,D] = eig(S);
[~,idx] = sort(diag(D),'descend');
V = V(:,idx(1:6));
U = Xc'*V;
U = U./sqrt(dot(U,U,1));
figure;
for i=1:6
    subplot(2,3,i);
    imagesc(reshape(U(:,i),112,92));
    title(['Eigenfaces ' num2str(i)]); 
    colormap('gray');
end