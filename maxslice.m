function z = maxslice(v)
%MAXSLICE This function finds a slice with largest tumor volume
sum_max = 0;
% if ~isvector(v)
%     error('Input must be a vector');
% end
for i = 1:5
    sum_temp = sum(sum(v(:,:,i)));
    if sum_max < sum_temp
        z = i;
        sum_max = sum_temp;
    end
end