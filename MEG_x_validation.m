function [out] = MEG_x_validation(X,Y,S,lambda,predictionType,classifier)

YSet = unique(Y);
 
%balance training set
if S.balanceTrainingSet
    class1Trials = shuffle(find(Y==YSet(1)));
    class2Trials = shuffle(find(Y==YSet(2)));
    
    classDiff = length(class1Trials) - length(class2Trials);
    
    removedTrials = [];
    %remove unbalanced trials
    if classDiff>0
        removedTrials = class1Trials(1:abs(classDiff));
        class1Trials(1:abs(classDiff)) = [];
    elseif classDiff<0
        removedTrials = class2Trials(1:abs(classDiff));
        class2Trials(1:abs(classDiff)) = [];
    end
    
    %create matched, shuffled cv index arrays
    cv_tr1 = shuffle(ceil(S.nXvals*(1:size(class1Trials,1))/size(class1Trials,1)));
    cv_tr2 = shuffle(cv_tr1);
    
    %combine the cv index arrays into a general cv index array
    cv_tr = zeros(size(Y));
    
    cv_tr(class1Trials) = cv_tr1;
    cv_tr(class2Trials) = cv_tr2;
    cv_tr(removedTrials) = -1*shuffle(ceil(S.nXvals*(1:size(removedTrials,1))/size(removedTrials,1)));
else
    cv_tr = shuffle(ceil(S.nXvals*(1:size(X,1))/size(X,1)));
end

yPred = nan(size(Y));
if length(YSet)==2
    yProb = nan(size(Y'));
else
    yProb = nan(size(Y,2), length(YSet));
end

idxAllTrain = (cv_tr>0);

%penalty range
for l = 1:length(lambda)
    thisL = lambda(l);
    thisLStr = num2str(thisL);
    
    for i=1:S.nXvals
        idxTR = (cv_tr~=i)&(cv_tr>0);
        idxTE = (abs(cv_tr)==i);
        
        XTrain = X(idxTR, :);
        YTrain = Y(idxTR);
        
        XTest = X(idxTE, :);
        YTest = Y(idxTE);
        
        %xval fold id, ensures that each sub-xval is also balanced.
        foldid = cv_tr(idxTR);
        
        %feature selection based on mutual information
        mi = nan(size(XTrain,2),1);
        if S.FS
            if size(XTrain,2)>S.nFeats
                for t = 1:size(XTrain,2)
                    mi(t) = mutualinfo(YTrain, XTrain(:,t));
                end
                
                mi_sorted = sort(mi, 'descend');
                thresh = mi_sorted(S.nFeats);
                XTrain = XTrain(:,mi>=thresh);
                XTest = XTest(:,mi>=thresh);
            end
        end
        
        if strcmp(classifier, 'svr')
            % svr
            trainOpts = [S.trainOptsSVM thisLStr];
            model = svm_train(YTrain, XTrain, trainOpts);
            [yPred(idxTE),~,yProb(idxTE)]  = svmpredict(YTest, XTest, model);
        elseif strcmp(classifier, 'ridge')
            %ridge regression
            trainOpts = [S.trainOptsSVM thisLStr];
            b = ridge(YTrain,XTrain,thisL, 0);
            yPred(idxTE) = [ones(size(XTest,1),1) XTest]*b;
        elseif strcmp(classifier, 'liblinear')
            trainOpts = [S.trainOptsLibLinear thisLStr];
            model = train(YTrain', XTrain, trainOpts);
            [yPred(idxTE),~,yProb(idxTE,:)] = predict(YTest', XTest, model);
        elseif strcmp(classifier, 'glmnet')
            
            % options
            glmnetoptions = glmnetSet();
            %glmnetoptions.alpha_set = S.glmnet.alpha_set;
            glmnetoptions.nlambda = S.glmnet.nlambda;
            
            % xval-within-an-xval on train data to determine parameters to use for actual classification
            
            opts.alpha = S.glmnet.alpha;
            glmnetoptions = glmnetSet(opts);
            
            m = cvglmnet(XTrain,YTrain,S.nXvals,foldid,'class','binomial',glmnetoptions,0);
            
            %m = EF_cvglmnetGrid(XTrain,YTrain,S.nXvals,[],'class','binomial',glmnetoptions,0);
            
            opts.lambda =  m.lambda_1se;
            %opts.alpha = m.alpha_min;
            
            glmnetoptions = glmnetSet(opts);
            
            % actual classification: train on XTrain, test on XTest.
            model = glmnet(XTrain, YTrain, 'binomial', glmnetoptions);
            
            %class guesses
            yPred(idxTE) = glmnetPredict(model, 'class', XTest);
            
            %logit of class guesses
            yProb(idxTE) = glmnetPredict(model, 'link', XTest);
        end
        fprintf('\n xval iteration %g complete', i);
    end
    
    % all data that went into training (should be balanced)
    XallTrain = X(idxAllTrain,:);
    YallTrain = Y(idxAllTrain);
    foldIdAllTrain = cv_tr(idxAllTrain);
    
    % create a model so we can analyze betas
    if strcmp(classifier, 'svr')
        model = svm_train(YTrain, XTrain, trainOpts);
    elseif strcmp(classifier, 'liblinear')
        model = train(Y', X, trainOpts);
    elseif strcmp(classifier, 'glmnet')
        
        glmnetoptions = glmnetSet();
        %glmnetoptions.alpha_set = S.glmnet.alpha_set;
        glmnetoptions.alpha = S.glmnet.alpha;
        glmnetoptions.nlambda = S.glmnet.nlambda;
        
        %m = EF_cvglmnetGrid(X,Y,S.nXvals,[],'class','binomial',glmnetoptions,0);
        m = cvglmnet(XallTrain,YallTrain,S.nXvals,foldIdAllTrain,'class','binomial',glmnetoptions,0);
        
        glmnetoptions = glmnetSet();
        glmnetoptions.lambda =  m.lambda_1se;
        glmnetoptions.alpha = S.glmnet.alpha;
        %glmnetoptions.alpha = m.alpha_min;
        
        model = glmnet(XallTrain, YallTrain, 'binomial', glmnetoptions);
    end
    
    %write this stuff out
    out.mod.model = model;
    out.mod.Y = Y;
    out.mod.YPred = yPred;
    out.mod.YProb = yProb;
    
    if strcmp(predictionType, 'continuous')
        [out.rObsVsPred(l), out.pObsVsPred(l)] = corr(Y,yPred);
    elseif strcmp(predictionType, 'discrete')
        accuracy= Y==yPred;
        for j=1:length(YSet)
            ThisY = YSet(j);
            classAcc(j) = mean(accuracy(Y==ThisY));
        end
        out.meanAccuracyOverClasses(l) = mean(classAcc);
        out.meanAccuracyEachClass(l,:) = classAcc;
    end
end
end