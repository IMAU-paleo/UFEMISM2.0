function isso = lies_on_line_segment( pa, pb, pc, tol_dist)
% Test if the point pc lies on the line pb-pc, or within tol_dist of it

d = pb - pa;
e = pc - pa;

d_norm = d / norm(d);

e_par = (e(1)*d_norm(1) + e(2)*d_norm(2)) * d_norm;
e_ort = e - e_par;

isso = true;
if (norm(e_ort) > tol_dist) 
  isso = false;
  return
end

if ((e(1)*d(1) + e(2)*d(2)) > 0) 
  if (norm(e_par) > (norm(d))) 
    isso = false;
    return
  end
else
  if (norm(e_par) > 0) 
    isso = false;
    return
  end
end

end