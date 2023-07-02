function isso = is_in_triangle( pa, pb, pc, p)
% Check if the point p lies inside the triangle abc, or within distance tol_dist of its edges    

tol = 1E-8;

% Check if p lies in the interior of the triangle abc
as_x = p(1)-pa(1);
as_y = p(2)-pa(2);

s1 = ((pb(1)-pa(1))*as_y-(pb(2)-pa(2))*as_x);
s2 = ((pc(1)-pa(1))*as_y-(pc(2)-pa(2))*as_x);
s3 = ((pc(1)-pb(1))*(p(2)-pb(2))-(pc(2)-pb(2))*(p(1)-pb(1)));

isso = false;   

if (s1 > -tol && s2 < tol && s3 > -tol) 
  isso = true;
end

end