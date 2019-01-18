function massive_to_write_ = massive_to_write_single(value, type1, type2)
%
s = size(value,2);
%
for i = 1:1:s
%     massive_to_write__ = typecast(swapbytes(typecast(value(i), 'uint8')), 'uint16');
    massive_to_write__ = typecast(swapbytes(typecast(value(i), type1)), type2);
    massive_to_write__ = [massive_to_write__(1) massive_to_write__(2)  ];
    if i ==1
        massive_to_write_ = massive_to_write__;
    else
        massive_to_write_ = [massive_to_write_ massive_to_write__];
    end
end
%    
end