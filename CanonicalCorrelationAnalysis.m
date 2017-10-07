data =[191 155 179 145;
195 149 201 152;
181 148 185 149;
183 153 188 149;
176 144 171 142;
208 157 192 152;
189 150 190 149;
197 159 189 152;
188 152 197 159;
192 150 187 151;
179 158 186 148;
183 147 174 147;
174 150 185 152;
190 159 195 157;
188 151 187 158;
163 137 161 130;
195 155 183 158;
186 153 173 148;
181 145 182 146;
175 140 165 137;
192 154 185 152;
174 143 178 147;
176 139 176 143;
197 167 200 158;
190 163 187 150];
n = size(data,1);
S  = (data -mean(data,1))'* (data -mean(data,1))/(n-1);
Sxx = S(1:2,1:2);
Syy = S(3:4,3:4);
Sxy = S(1:2,3:4);
Syx = Sxy';
W1 = Sxx\Sxy;
W2 = Syy\Syx;
[P, Tp] = eig(W1*W2);
[eigp,idxp] = sort(diag(Tp),'descend');
 P = P(:,idxp(1:2));
 [Q, Tq] = eig(W2*W1);
[eigq,idxq] = sort(diag(Tq),'descend');
 Q = Q(:,idxq(1:2));
fprintf('Sample matrix:\n');
 for i=1:2
     if P(:,i)'*Sxy*Q(:,i) < 0
           P(:,i) = -P(:,i);
     end
     P(:,i) = P(:,i)/sqrt(P(:,i)'*Sxx*P(:,i));
     Q(:,i) = Q(:,i)/sqrt(Q(:,i)'*Syy*Q(:,i));
     fprintf('eig %d:',i);
     display(eigp(i));
     fprintf('a:');
     display(P(:,i));
      fprintf('b:');
     display(Q(:,i));
 end
 
 D =  (data -mean(data,1));
 U = P'*D(:,1:2)';
 V = Q'*D(:,3:4)';
 figure();
 scatter(U(1,:),V(1,:),n);
figure();
 scatter(U(2,:),V(2,:),n);
v = 1./sqrt(diag(S));
S = diag(v)*S*diag(v);
Sxx = S(1:2,1:2);
Syy = S(3:4,3:4);
Sxy = S(1:2,3:4);
Syx = Sxy';
W1 = Sxx\Sxy;
W2 = Syy\Syx;
[P, Tp] = eig(W1*W2);
[eigp,idxp] = sort(diag(Tp),'descend');
 P = P(:,idxp(1:2));
 [Q, Tq] = eig(W2*W1);
[eigq,idxq] = sort(diag(Tq),'descend');
 Q = Q(:,idxq(1:2));
 fprintf('Correlation matrix:\n');
 for i=1:2
     if P(:,i)'*Sxy*Q(:,i) < 0
           P(:,i) = -P(:,i);
     end
     P(:,i) = P(:,i)/sqrt(P(:,i)'*Sxx*P(:,i));
     Q(:,i) = Q(:,i)/sqrt(Q(:,i)'*Syy*Q(:,i));
     fprintf('eig %d:',i);
     display(eigp(i));
     fprintf('a:');
     display(P(:,i));
      fprintf('b:');
     display(Q(:,i));
 end
 
 
 Vx = P'*Sxx;
 fprintf('r1:\n');
 s = Vx(1,:)*Vx(1,:)';
 display(s);
 fprintf('percentage1:\n');
 pv = s/2;
 display(pv);
 Vy = Q'*Syy;
 fprintf('r1:\n');
 s = Vy(1,:)*Vy(1,:)';
 display(s);
 fprintf('percentage1:\n');
 pv = s/2;
 display(pv);
 fprintf('eig percentage:\n');
 display(eigq(1)/(eigq(1) + eigq(2)));
