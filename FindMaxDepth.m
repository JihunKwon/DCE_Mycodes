function MaxDepth = FindMaxDepth(Net,k,i)

MaxDepth=0;
if Net(k,1,i-1) > Net(k,1,i)
    MaxDepth=i+1;
end

end

