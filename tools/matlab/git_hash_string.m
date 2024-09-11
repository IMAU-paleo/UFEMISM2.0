function str = git_hash_string
[~,str] = system('git rev-parse HEAD');
str = str(1:end-1);
end