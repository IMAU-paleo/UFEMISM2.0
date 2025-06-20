function r95 = calc_95th_percentile( d1, d2)
% Calculate the 95th percentile between d2 and d1

rr = abs( d2 - d1);
rr = sortrows( rr);

i95 = ceil( length( rr) * 0.95);
r95 = rr( i95);

end