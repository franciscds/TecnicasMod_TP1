function eta = lambda2eta_sobre(lambda)

tmp_eta=0:0.01:1;
x=log(tmp_eta)./(tmp_eta-1);
tmp_lambda=x.*exp(-x);
[void, pos]=min(abs(lambda-tmp_lambda));
eta=tmp_eta(pos);
end
