function massive_to_write_ = massive_to_write(value, type1, type2)
%
s = size(value,1);
%
for i = 1:1:s
%     massive_to_write__ = typecast(swapbytes(typecast(value(i), 'uint8')), 'uint16');
    massive_to_write__ = typecast(swapbytes(typecast(value(i), type1)), type2);
    massive_to_write__ = [massive_to_write__(1) massive_to_write__(2) massive_to_write__(3) massive_to_write__(4)];
    if i ==1
        massive_to_write_ = massive_to_write__;
    else
        massive_to_write_ = [massive_to_write_ massive_to_write__];
    end
end
%    
end