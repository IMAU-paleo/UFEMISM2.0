fid = fopen('foo/bar.txt','w');
fprintf(fid,'%s\n','Hello file!');
fclose(fid);