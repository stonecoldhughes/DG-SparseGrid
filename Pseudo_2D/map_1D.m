function [ z ] = Map_1D( level, cell, poly, k )

if level==0
    z=poly;
else
    z=k*2^(level-1)+k*cell+poly;
end

end

