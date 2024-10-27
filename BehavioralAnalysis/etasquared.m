ngroup=9;
% nobs=length(datavec{1})+length(datavec{2})+length(datavec{3})+length(datavec{4})+length(datavec{5})
% nobs=28+52+53
% nobs=19+34+27
nobs=288;
eta_squared=(cell2mat(tbl{3}(2,5))-ngroup+1)./(nobs-ngroup)