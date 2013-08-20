function [trl, event] = detectTrigger(cfg)

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);
%event = ft_read_event(cfg.dataset,'trigindx',cfg.trialdef.trigchannel,'threshold',3,'detectflank','down');
event = ft_read_event(cfg.dataset,'trigindx',str2num(cfg.trialdef.trigchannel),'threshold',3,'detectflank','down');

% search for "trigger" events according to 'trigchannel' defined outside the function
value  = [event(find(strcmp(cfg.trialdef.trigchannel, {event.type}))).value]';
sample = [event(find(strcmp(cfg.trialdef.trigchannel, {event.type}))).sample]';


pretrig = cfg.trialdef.prestim*-1000;
posttrig = cfg.trialdef.poststim*1000;

% pretrig = -250;
% posttrig = 500;
trl = [];

%load in behavioral data for sub mem analysis
load(['behav/encodingDataFileSub' cfg.sub 'Block' num2str(cfg.block)])
load(['behav/testDataFileSub' cfg.sub 'Block' num2str(cfg.block)])

% make a vector representing memory accuracy
memVec = zeros(23,1);
memVec(strmatch('1!',[testDataFile(1:end,8)]))=1;
memVec(strmatch('2@',[testDataFile(1:end,8)]))=2;
memVec(strmatch('3#',[testDataFile(1:end,8)]))=3;
memVec(strmatch('4$',[testDataFile(1:end,8)]))=4;
memVec(1)=[];

memVec(:,2)=cell2mat(testDataFile(2:end,12)); %#ok<COLND>

% calculate subsequent memory for each test item
for imem=1:length(memVec)
    
    if memVec(imem,1)==1&&memVec(imem,2)==0||memVec(imem,1)==4&&memVec(imem,2)==1
        
        memVec(imem,3)=1;
        
    elseif memVec(imem,1)==2&&memVec(imem,2)==0||memVec(imem,1)==3&&memVec(imem,2)==1
        
        memVec(imem,3)=2;
        
    elseif memVec(imem,1)==1&&memVec(imem,2)==1||memVec(imem,1)==4&&memVec(imem,2)==0
        
        memVec(imem,3)=3;
        
    elseif memVec(imem,1)==2&&memVec(imem,2)==1||memVec(imem,1)==3&&memVec(imem,2)==0
        
        memVec(imem,3)=4;
        
    end
    
end


% make the order condition labels unique
memVec(12:end,3)=memVec(12:end,3)+4;

memVec=[0 0 0; memVec];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Subsequent Memory Condition Labels%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1 = Source HC correct
% 2 = Source LC correct
% 3 = Source HC incorrect
% 4 = Source LC incorrect
% 5 = Order HC correct 1st item
% 6 = Order LC correct 1st item
% 7 = Order HC incorrect 1st item
% 8 = Order LC incorrect 1st item
% 9 = Order HC correct 2nd item
% 10 = Order LC correct 2nd item
% 11 = Order HC incorrect 2nd item
% 12 = Order LC incorrect 2nd item

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for imem = 1:length(testDataFile)
    
    if imem==1
        
        
    elseif imem<=12
        
        tmpMemVec(encodingDataFile{strncmp(testDataFile{imem,2},[encodingDataFile{2:end,2}],3),1}+1)=memVec(imem,3);
        
        
    elseif imem>12
        
        tmpMemVec(encodingDataFile{strncmp(testDataFile{imem,3},[encodingDataFile{2:end,2}],3),1}+1)=memVec(imem,3);
        tmpMemVec(encodingDataFile{strncmp(testDataFile{imem,4},[encodingDataFile{2:end,2}],3),1}+1)=memVec(imem,3)+4;
        
    end
    
end

tmpMemVec=tmpMemVec';

bIdx = 1:6:36;
bIdx(1)=[];
nbIdx = 4:6:36;
tmpMemVec(bIdx,2)=tmpMemVec(bIdx-2,1);
tmpMemVec(nbIdx,2)=tmpMemVec(nbIdx-2,1);

nbOrderIdx1 = [2:6:36];
nbOrderIdx2 = [6:6:36];

bOrderIdx1 = [5:6:36];
bOrderIdx1(6)=[];
bOrderIdx2 = [3:6:36];
bOrderIdx2(1)=[];

tmpMemVec(nbOrderIdx1,2)=tmpMemVec(nbOrderIdx1+2,1);
tmpMemVec(nbOrderIdx2,2)=tmpMemVec(nbOrderIdx2-2,1);

tmpMemVec(bOrderIdx1,2)=tmpMemVec(bOrderIdx1+2,1);
tmpMemVec(bOrderIdx2,2)=tmpMemVec(bOrderIdx2-2,1);



% creating your own trialdefinition based upon the events
for j = 1:length(value);
    trlbegin = sample(j) + pretrig;
    trlend   = sample(j) + posttrig;
    offset   = pretrig;
    newtrl   = [ trlbegin trlend offset j];
    trl      = [ trl; newtrl];
end

% % subject specific changes
% if cfg.block==1
% 
%     trl = [trl tmpMemVec(2:36,:)];
% 
% else
% 
 trl = [trl tmpMemVec];
% 
% end
% 

tmp = trl(find(mod(trl(:,4),6)==2),5);


eventSubmem = [repmat(tmp(1),6,1); repmat(tmp(2),6,1); repmat(tmp(3),6,1); ...
    repmat(tmp(4),6,1); repmat(tmp(5),6,1); repmat(tmp(6),6,1)];

% if cfg.block==1
% 
% trl = [trl eventSubmem(2:36,:)];
% 
% else

trl = [trl eventSubmem];

if any(trl(:,1)<0)
    
    trl(find(trl(:,1)<0),:)=[];
    
end

% end

