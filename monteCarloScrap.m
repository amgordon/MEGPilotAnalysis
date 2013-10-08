clear acc
for i=1:1000
    for j=1:10
        theseGuesses = ceil(rand(size(Y))*3);
        for k=1:3
            acc_h(k) = mean(Y(Y==k)==theseGuesses(Y==k));
        end
        acc(i,j) = mean(acc_h);
    end
end