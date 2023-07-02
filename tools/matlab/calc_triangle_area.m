function TriA = calc_triangle_area( pq, pr, ps)
% Find the area of the triangle [pq,pr,ps]

TriA = abs( cross2( [pr(1)-pq(1), pr(2)-pq(2)], [ps(1)-pq(1), ps(2)-pq(2)] )) / 2;

end