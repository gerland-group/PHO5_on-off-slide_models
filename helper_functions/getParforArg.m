function arg = getParforArg()  
% from https://de.mathworks.com/matlabcentral/answers/196411-run-parfor-as-for-loop#answer_174245
    p = gcp('nocreate'); % get the pool object if it exists, but never open a pool
    if isempty(p)
      arg = 0;
    else
      arg = p.NumWorkers;
    end
end
