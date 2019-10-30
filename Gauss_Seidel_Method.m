function ret = Gauss_Seidel_Method(A,B)
[m,n]=size(A);
for i= 1:m
    diag(i)=A(i,i);
    rest(i,:)=[A(i,1:i-1) A(i,i+1:n)];
    row_sum = sum(abs(rest(i,:)));
    
    if(abs(diag(i))<=row_sum)
        ret='Dominance condition is not satisfied';
        return;
    end
    
end

X=zeros(1,n);
ret(1,:)= X;
for i=2:20
    for j=1:n
        Xrest = [X(1:j-1) X(j+1:n)];
        X(j)=(B(j)-rest(j,:)*Xrest')/diag(j);
    end
    
    if(abs(X(1) - ret(i-1,1)) <= 10^-4)
        ret(i,:) = X;
        return;
    end

    ret(i,:)=X;
    
end

end