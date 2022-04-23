function [count] = SetBits(n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        count = uint8(0); 
        while n~=0
                n = bitand(n,n-1);
                count = count + 1;
        end
end

