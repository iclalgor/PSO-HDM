function y = dogrusalmi(A,B,C)
    u1 = (B-A)/(pdist2(A, B, 'euclidean' ));
    u2 = (C-B)/(pdist2(C, B, 'euclidean' ));
    y = sum(u1==u2)==size(u1,2);
end