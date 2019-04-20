function VisualizeMatches(I1, I2, matches, f1,f2,linecolor, featurecolor)

figure; clf ;
% imagesc(cat(2, I1, I2)) ;
imshow(cat(2, I1, I2));

for k = 1:size(matches,1)
    xa = f1(1,matches(k,1));
    xb = f2(1,matches(k,2)) + size(I1,2);
    
    ya = f1(2,matches(k,1));
    yb = f2(2,matches(k,2));

    hold on ;
    h = line([xa ; xb], [ya ; yb]) ;
    set(h,'linewidth', 1, 'color', linecolor) ;

    plot(xa, ya, [featurecolor 'o'], 'markerfacecolor',featurecolor, 'markerSize', 3);    hold on;
    plot(xb, yb, [featurecolor 'o'], 'markerfacecolor',featurecolor, 'markerSize', 3);    hold on;
%     pause;
end