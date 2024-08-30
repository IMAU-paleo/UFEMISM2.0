function map = bluewhiteredmap( n)

  ncolors=5;
  c=[0 0 0.5;
   0 0.5 1;
   1 1 1;
   1 0 0;
   0.5 0 0];
  pp=1:(n-1)/(ncolors-1):n;
  r=interp1(pp,c(:,1),1:n);
  g=interp1(pp,c(:,2),1:n);
  b=interp1(pp,c(:,3),1:n);
  map=[r' g' b'];
    
end