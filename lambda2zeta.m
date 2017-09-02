function zeta = lambda2zeta(lambda)
tmp_zeta=0:0.0001:1;
xxx=acos(tmp_zeta)./sqrt(1-tmp_zeta.^2);
tmp_lambda=xxx.*exp(-tmp_zeta.*xxx);
[void, pos]=min(abs(lambda-tmp_lambda));
zeta=tmp_zeta(pos);
