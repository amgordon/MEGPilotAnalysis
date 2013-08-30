function [] = MEG_splitSqdIntoRuns(par)
% Split a single .sqd file into N*3 sqd files, where N is the number of
% miniblocks run for the 3 (Study1, Study2, Test) tasks

tasks = {'Study1' 'Study2' 'Test'};
taskTrialDur = [4000, 4000, 5083];

%% load behavioral data
d = dir(fullfile(par.behavDir, 'MEG*StudyTest_study.mat'));
R = load(fullfile(par.behavDir, d.name));

for i = 1:length(par.dataFiles)
    
    thisDataFile = par.dataFiles{i};
    [~, thisDataName] = fileparts(thisDataFile);
    
    disp('Finding trigger lines and triggers...');
    [~, ~, triggersList, ~] = sqdgettriggers(thisDataFile,par.triggerRange);
    
    
    %% sqd-specifc info
    if strcmp(thisDataName, 'R0487_MEGclassblock2_8.21.13.sqd')
        nSqdsAlreadyWritten = 12;
        nTrialsAlreadyWritten = 144;
    elseif strcmp(thisDataName, 'R0487_MEGclassblock3_8.21.13')
        nSqdsAlreadyWritten = 66;
        nTrialsAlreadyWritten = 792;
    else
        nSqdsAlreadyWritten = 0;
        nTrialsAlreadyWritten = 0;
    end
    
    %% use behavioral data to find where runs begin and end
    
    idx.allTriggers = sort(vertcat(triggersList{:}));
    
    allTrials = (1:(length(idx.allTriggers)/3)) + nTrialsAlreadyWritten/3;
    idx.MB = cell2mat(R.theData.Study1.miniblock(allTrials));
    idx.firstInMB = logical([true diff(idx.MB)]);
    idx.lastInMB = logical([diff(idx.MB) true]);
    
    idx.firstInMBAcrossTasks = repmat(idx.firstInMB,1,3);
    idx.lastInMBAcrossTasks = repmat(idx.lastInMB,1,3);
    
    idx.triggerFirstInMiniBlock = idx.allTriggers(idx.firstInMBAcrossTasks);
    idx.triggerLastInMiniBlock = idx.allTriggers(idx.lastInMBAcrossTasks);
    
    task_code_h = mod(1:length(idx.triggerFirstInMiniBlock),3);
    task_code_h(task_code_h==0) = 3;
    idx.taskCode = task_code_h; % 1 for Study1, 2 for Study2, 3 for Test.
    
    % for each miniblock of each task
    for i=1:length(idx.triggerFirstInMiniBlock)
        clear trigList
        trigList{1} = idx.triggerFirstInMiniBlock(i);
        nPretriggerSamples = par.nPretriggerSamplesForRun;
        nPosttriggerSamples = idx.triggerLastInMiniBlock(i)  - idx.triggerFirstInMiniBlock(i) + taskTrialDur(idx.taskCode(i));
        outputFilenameRoot = strrep(thisDataFile,'.sqd','.epoch.mat');
        epochAvgFlag = 'epoch';
        
        % read out the data from the old .sqd file
        [epochFiles, epochInfo] = sqdmakeepochs(thisDataFile, trigList, nPretriggerSamples, nPosttriggerSamples, outputFilenameRoot, epochAvgFlag);
        disp(epochInfo)
        
        y = load(epochFiles{1});
        
        % write out the new .sqd file
        thisTask = tasks{idx.taskCode(i)};
        thisNum = ceil((i+nSqdsAlreadyWritten)/3);
        thisSqd = fullfile(par.runsDir, ['MEGDat' '_' thisTask '_' prepend(thisNum) '.sqd']);
        
        % delete previous versions of the file, otherwise sqdwrite will want to
        % append to them.
        if exist(thisSqd, 'file')
            delete(thisSqd)
        end
        
        sqdwrite(thisDataFile, thisSqd, y.data)
    end
end
end

