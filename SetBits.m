function [count] = SetBits(n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        count = 0; %uint8(0); Changed to be able to do arithmetic operetions with it with other numbers i.e. division with double
        while n~=0
                n = bitand(n,n-1);
                count = count + 1;
        end
end

