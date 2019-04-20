function Y = LinearTriangulation_InverseDepth(matches, P1, P2)

Y = zeros(3, size(matches,1));
for i = 1:size(matches,1)
    b1 = [matches(i, 1:2) 1]';
    b1 = b1./norm(b1);
    b2 = [matches(i, 4:5) 1]';
    b2 = b2./norm(b2);
%     A = [skewsymm(b1)*P1; skewsymm(b2)*P2];
%     [~,~,V] = svd(A);
    A = [skewsymm(b1)*P1(:,1:3); skewsymm(b2)*P2(:,1:3)];
    r = -[skewsymm(b1)*P1(:,4); skewsymm(b2)*P2(:,4)];    
    V = (A'*A)\(A'*r);        
    rho = 1/norm(V(1:3));
%     if V(4) < 0
%         V = -V;
%     end
%     rho = (1/norm(V(1:3,end))^2-1)^0.5;
    
    theta = atan2(V(2),V(1));
    phi = atan2(V(3),(V(1)^2+V(2)^2)^0.5);

    Y(:,i) = [theta;phi;rho];
    
end

end