% i=1:100;
% j=1:100;
Y=0;
parfor i=1:10
    for j=1:1:10
        Y=Y+i+j;
    end    
end
