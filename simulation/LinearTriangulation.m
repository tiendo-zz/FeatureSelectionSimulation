function Y = LinearTriangulation(matches, P1, P2)

Y = zeros(3, size(matches,1));
for i = 1:size(matches,1)
    b1 = [matches(i, 1:2) 1]';
    b1 = b1./norm(b1);
    b2 = [matches(i, 4:5) 1]';
    b2 = b2./norm(b2);
%     A = [skewsymm([matches(i, 1:2) 1]')*P1; skewsymm([matches(i, 4:5) 1]')*P2];
    A = [skewsymm([b1;1])*P1; skewsymm([b2;1])*P2];
%     A = [skewsymm(b1)*P1(:,1:3); skewsymm(b2)*P2(:,1:3)];
%     r = -[skewsymm(b1)*P1(:,4); skewsymm(b2)*P2(:,4)];    
%     V = (A'*A)\(A'*r);  
%     Y(:,i) = V;
    [~,~,V] = svd(A);
    Y(:,i) = V(1:3,end)/V(4,end);
end

end