function c = myismember(a,b)
    
    input_empty = zeros(2,1);
    if ~isrow(a) && isrow(a')
        a = a';
    elseif isrow(a)
    else
        if isempty(a)
            input_empty(1) = 1;
        else
            error('myismember only works with vectors!')
        end
    end
    if ~isrow(b) && isrow(b')
        b = b';
    elseif isrow(b)
    else
        if isempty(b)
            input_empty(2) = 1;
        else
            error('myismember only works with vectors!')
        end
    end
    if sum(input_empty)==0
        A = ones(length(b),1)*a;
        B = b'*ones(1,length(a));
        c = sum(A-B==0,1)>0;
    else
        c = [];
    end
end
