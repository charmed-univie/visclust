function [value] = evaluation(vect1,vect2,mode)
% [VALUE] = EVALUATION(VECT1, VECT2, MODE) computes an error measure
% determined by MODE of the vectors VECT1 and VECT2 containing clustering
% indices. Possible error measure (parameter MODE):
%    'ACC' - Accuracy
%    'ARI' - Adjusted Rand Index
%    'RI' - Rand Index

if mode == "ACC"
    % create confusion matrix
    numclasses=max(length(unique(vect2)),length(unique(vect1)));
    CM=confusionmat(vect1,vect2);
    zerocol=find(all(CM == 0,2));
    zerorow=find(all(CM == 0,1));
    % delete empty rows and columns in confusion matrix
    if ~isempty(zerocol)&&~isempty(zerorow)
        numdel=min(length(zerocol),length(zerorow));
        CM(:,zerorow(1:numdel))=[];
        CM(zerocol(1:numdel),:)=[];
    end
    % iterate through all permutations of the confusion matrix to calculate the optimal accuracy
    classperms=perms(1:numclasses);
    value=0;
    for j=1:length(classperms)
        predsum=0;
        for i=1:numclasses
            predsum=predsum+CM(i,classperms(j,i));
        end
        naccuracy=predsum/length(vect2);
        if naccuracy > value
            value=naccuracy;
        end
    end
else
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
    if mode=="RI"
        tpfpfntn=nchoosek(length(vect2),2);
        tn=tpfpfntn-tp-fp-fn;
        value=(tp+tn)/(tp+tn+fp+fn);
    elseif mode=="ARI"
        n=length(vect2);
        if fn == 0 && fp == 0
            value=1;
        else
            value=(tp-(tpfn*tpfp)/(nchoosek(n,2)))/((1/2)*(tpfn+tpfp)-(tpfn*tpfp)/(nchoosek(n,2)));
        end
    end
end
end