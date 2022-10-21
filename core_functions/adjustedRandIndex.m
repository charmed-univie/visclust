function value = adjustedRandIndex(vect1,vect2)
% [VALUE] = EVALUATION(VECT1, VECT2) computes the adjusted Rand index 
% of the vectors VECT1 and VECT2 containing clustering indices.
%
%
% This is part of the clustering algorithm VISCLUST, written by Anna Breger 
% and Clemens Karner.
% University of Vienna, Faculty of Mathematics
% Vienna, Austria
% Copyright (c) 2022
% https://homepage.univie.ac.at/anna.breger/
% https://homepage.univie.ac.at/clemens.karner/
%
% For all questions, bugs and suggestions please email
% clemens.karner@univie.ac.at or anna.breger@univie.ac.at
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create confusion matrix
CM=confusionmat(vect2,vect1);
tpfpm=sum(CM,2);
% calculate true positive and false positive entries
tpfp=0;
for i=1:length(tpfpm)
    if tpfpm(i) > 1
        tpfp=tpfp + nchoosek(tpfpm(i),2);
    end
end
% calulate true positive entries
tp=0;
for i=1:size(CM,1)
    for j=1:size(CM,2)
        if CM(i,j)>1
            tp=tp+nchoosek(CM(i,j),2);
        end
    end
end
% calculate false positive entries
fp=tpfp-tp;
tpfnm=sum(CM,1);
% calculate true positive false negative entries
tpfn=0;
for i=1:length(tpfnm)
    if tpfnm(i) > 1
        tpfn=tpfn + nchoosek(tpfnm(i),2);
    end
end
% calculate false negative entries
fn=tpfn-tp;
% combine to adjusted Rand index
n=length(vect2);
if fn == 0 && fp == 0 % Special cases: empty data or full agreement
    value=1;
else
    value=(tp-(tpfn*tpfp)/(nchoosek(n,2)))/((1/2)*(tpfn+tpfp)-(tpfn*tpfp)/(nchoosek(n,2)));
end
end
