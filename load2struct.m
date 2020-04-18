function s= load2struct(fn, s)

run(fn)
vars= setdiff(who(), {'fn', 's'});

if ~exist('s', 'var')
    s= [];
end


for i= 1:length(vars)
    if strcmp(vars{i}, 'fn'), continue; end
    
    eval(sprintf('s.%s= %s;', vars{i}, vars{i}));
end