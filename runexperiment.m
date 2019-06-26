%% Novel objects in fMRI 
% Tijl Grootswagers

%% set up psychtoolbox if needed (should only be needed for testing)
if isempty(which('PsychtoolboxVersion.m'))
    addpath(genpath('~/PHD/MatlabToolboxes/Psychtoolbox/'))
end
close all; %close all figure windows
clearvars; %clear all variables
p=struct();

%% START define parameters

%debugging/testing parameters: CHECK THESE!
p.isrealexperiment = 0; %should only be 0 for testing; skips overwrite checks and doesn't ask for subject nr
p.fullscreen = 0; %should only be 0 for testing; if 0 it does not create a fullscreen window
p.isFMRIexperiment = 0; %should only be 0 for testing; does not send or wait for triggers
p.synctest = 0; %should only be 0 for testing; skips synchronisation tests
p.this_is_not_a_drill = 0; %should only be 0 for practice; gives feedback and exits after 4 blocks

%timing parameters
p.stimulusduration = 0.200; %stimulus duration (secs);
p.responseduration = 2; %response duration from stimulus onset (secs);
p.halfrefresh = .5/60; %to subtract from fliptimes
p.blockstartfixationduration = 2; %fixation duration at start of block (before first stim)
p.block_num_message_time = 20; % duration of break between blocks (secs)
p.feedbackduration = 16; %minimum time to leave feedback on screen to capture full HRF (secs)

%display parameters
p.fixationsize = 20; %diameter of fixation cross (pixels)
p.backgroundcolor = [100 100 100]; %backround (same grey as stimulus background)
p.fontsize = 20; %font size for all the things
p.windowsize = [0 0 800 600]; %windowsize if not running fullscreen
p.stimulussize = [0 0 201 201]; %size of stimulus in pixels

%responsekey parameters (eg 1! and 2@ are the response buttons in MEG)
p.responsekeysleft = [KbName('1') KbName('1!')];
p.responsekeysright = [KbName('2') KbName('2@')];

% fMRI parameters
p.tr = 1.1135; % TR in s
p.num_baseline_triggers = floor(p.block_num_message_time/p.tr)-1; % Automaticall calculate the number of triggers we record as a baseline at the start of each run

%% END define parameters

% repetition numbers
p.nruns = 8; %number of runs (this cannot be changed easily)
p.trialsperrun = 144; %how many trials per block (this cannot be changed easily)

%% numbers
% for x = [24 26 28 30]
% p.trialsperrun = (x^2)*2 / 8;

%estimates
p.estimatedrunduration_minutes = p.trialsperrun * p.responseduration/60;
p.estimatedexperimentduration_minutes = p.nruns * p.estimatedrunduration_minutes;
% 
% fprintf('%i %.2fm %.2fm (%i stimuli)\n',x, p.estimatedrunduration_minutes, p.estimatedexperimentduration_minutes, x^2)
% end

% set stimulus parameters
p.stimuli = struct();i=0;
for len = [-12:-1 1:12]
    for ori = [-12:-1 1:12]
        if p.this_is_not_a_drill
            fn = dir(sprintf('exp_stimuli%s*len_%02i_ori_%02i.png',filesep,len,ori));
        else
            fn = dir(sprintf('exp_stimuli_practice%s*len_%02i_ori_%02i.png',filesep,len,ori));
        end
        assert(~isempty(fn),'stimuli not found')
        i=i+1;
        p.stimuli(i).stimnum = i;
        p.stimuli(i).len = len;
        p.stimuli(i).ori = ori;
        p.stimuli(i).fn = fn.name;
        p.stimuli(i).fol = fn.folder;
        p.stimuli(i).run = 1+mod(abs(len)-abs(ori),4);
    end
end
p.nstimuli = length(p.stimuli);
assert(p.nstimuli==576)

%% get subject info
if ~p.isrealexperiment
    p.subjectnr = 0;
    disp('    +-----------------------------------------+')
    disp('    | WARNING: p.isrealexperiment is set to 0 |')
    disp('    | If this is a real fMRI run, set to 1    |')
    disp('    +-----------------------------------------+')
else
    p.subjectnr = input('\n    Subject number: ','s');
    p.subjectnr = str2double(p.subjectnr);
    if isnan(p.subjectnr)
        error('ERROR: invalid subject number');
    end
end
if p.this_is_not_a_drill
    p.datafilename = sprintf('sub-%02i_task-categorisation_events.mat',p.subjectnr);
    p.csvdatafilename = sprintf('sub-%02i_task-categorisation_events.csv',p.subjectnr);
else
    disp('    +--------------------------------------------+')
    disp('    | WARNING: p.this_is_not_a_drill is set to 0 |')
    disp('    | If this is a real fMRI run, set to 1       |')
    disp('    +--------------------------------------------+')
    p.datafilename = sprintf('sub-%02i_task-categorisation_practice.mat',p.subjectnr);
    p.csvdatafilename = sprintf('sub-%02i_task-categorisation_practice.csv',p.subjectnr);
end

%check if we are possibly overwriting data
if p.isrealexperiment
    if exist(p.datafilename,'file')
        error(['ERROR: ' p.datafilename ' exists. Move (or delete) this file or use a different subject number']);
    end
    if exist(p.csvdatafilename,'file')
        error(['ERROR: ' p.csvdatafilename ' exists. Move (or delete) this file or use a different subject number']);
    end
end

%This should only be used to test/debug. CHECK THIS!!
Screen('Preference', 'SkipSyncTests', double(~p.synctest));
if ~p.synctest
    disp('    +-----------------------------------------+')
    disp('    | WARNING: p.synctest is set to 0         |')
    disp('    | If this is an fMRI experiment, set to 1 |')
    disp('    +-----------------------------------------+')
end

% set up seed
if p.this_is_not_a_drill
    p.randomseed = rng(p.subjectnr);
else
    p.randomseed = rng(100+p.subjectnr); %use a different seed for practice
end

%% open i/o port
if ~p.isFMRIexperiment
    disp('    +-----------------------------------------+')
    disp('    | WARNING: p.isFMRIexperiment is set to 0 |')
    disp('    | If this is an fMRI experiment, set to 1 |')
    disp('    +-----------------------------------------+')
end

%make a copy of the raw experiment code and store it so we can recreate the exact script
p.experimentcode = fileread('runexperiment.m');

%% create trialstructure for all runs and blocks

% create yoked stimulus orders (for the two dimensions)
eventlist = table();
p.rundims = repmat({'len','ori'},1,4);
p.runqs = {};
for i = 1:2
    p.runqs = [p.runqs reshape([randsample({'shorter','longer'},2); randsample({'pointed up','pointed down'},2)],1,[])];
end
if mod(p.subjectnr,2) %balance the dimension order across pp
    p.rundims = fliplr(p.rundims);
    p.runqs = fliplr(p.runqs);
end
p.runstims = ones(1,8);
p.runstims(1:2:8) = randsample(1:4,4);
while ~all(diff(p.runstims)) || ~all(all(diff(sort(reshape(p.runstims,4,[])))==1))
    p.runstims(2:2:8) = randsample(1:4,4);
end

for runnr = 1:p.nruns
    % find stimuli to present this run
    stimulithisrun = [p.stimuli([p.stimuli(:).run]==p.runstims(runnr)).stimnum];
    sequence = randsample(stimulithisrun',length(stimulithisrun));
    % add sequence to trialstruct
    r=struct();
    r.eventnumber = size(eventlist,1) + (1:length(sequence))';
    r.runnumber = runnr*ones(size(sequence));
    r.presentationnumber = (1:length(sequence))';
    r.stimset = p.runstims(runnr)*ones(size(sequence));
    r.stimnumber = sequence;
    r.stimname = {p.stimuli(sequence).fn}';
    r.len = [p.stimuli(sequence).len]';
    r.ori = [p.stimuli(sequence).ori]';
    cr = {'left','right'};
    if ismember(p.runqs{runnr},{'shorter','pointed up'})
        cr = fliplr(cr);
    end
    if strcmp(p.rundims(runnr),'len')
        r.taskdim = r.len;
    else
        r.taskdim = r.ori;
    end
    r.dimension = repmat(p.rundims(runnr),length(sequence),1);
    r.question = repmat(p.runqs(runnr),length(sequence),1);
    r.correctresp = cr(1 + (r.taskdim < 0))';
    %concat
    eventlist = [eventlist;struct2table(r)]; %#ok<AGROW>
end
% coding of reponse questions is a little confusing so here we check that 
% the responses fit (e.g., -12,-12 should be right longer & right pointed down)
idx = strcmp(eventlist.dimension,'len') & strcmp(eventlist.question,'shorter') & eventlist.len<0;
assert(all(strcmp(eventlist.correctresp(idx),'left')))
idx = strcmp(eventlist.dimension,'len') & strcmp(eventlist.question,'longer') & eventlist.len<0;
assert(all(strcmp(eventlist.correctresp(idx),'right')))
idx = strcmp(eventlist.dimension,'ori') & strcmp(eventlist.question,'pointed up') & eventlist.ori<0;
assert(all(strcmp(eventlist.correctresp(idx),'left')))
idx = strcmp(eventlist.dimension,'ori') & strcmp(eventlist.question,'pointed down') & eventlist.ori<0;
assert(all(strcmp(eventlist.correctresp(idx),'right')))

%if practice, only do the first 4 blocks, and exit nicely (without errors)
if ~p.this_is_not_a_drill 
    eventlist = eventlist(eventlist.runnumber<=4,:);
end

%% setup abort function
abort = ['sca;Priority(0);ListenChar(0);ShowCursor();fprintf(''\n\nABORTING...\n'');writetable(eventlist,p.csvdatafilename);save(''aborted.mat'');sca;'...
    'error(struct(''stack'',struct(''file'',''runexperiment'',''name'',''runexperiment'',''line'',0),'...
    '''message'',''############################## EXPERIMENT ABORTED BY USER ##############################''));'];

%% disable keyboard input
ListenChar(2);
KbName('UnifyKeyNames') %cause old version of psychtoolbox
KbCheck(); % make sure kbcheck is compiled so it is always accurate on first call

%% summarize parameters
while KbCheck()
    WaitSecs(0.01);
end
disp('    +------------------------+')
disp('    | Experiment parameters: |')
disp('    +------------------------+')
disp(p);
disp('If this is correct, press <SPACE> to start.')
[~, keycode, ~] = KbWait(); %wait for a keypress
if ~keycode(KbName('space'));eval(abort);end % check for return key
% Wait for key release

%% open window, and wait for the photodiode setup
p.screennumber=min (Screen('Screens'));
if p.screennumber>1 && p.isrealexperiment
    ListenChar(0);error('Found too many screens! Exit matlab, setup mirrored-display mode, and try again')
end
p.black=BlackIndex(p.screennumber);
p.gray=GrayIndex(p.screennumber);
p.white=WhiteIndex(p.screennumber);
if p.fullscreen
    p.windowsize=[];
elseif ismac && max(Screen('Screens',1))==2
    p.windowsize=[];
end
[p.window, p.windowrect]=Screen('OpenWindow',p.screennumber, p.backgroundcolor,p.windowsize,32,2);
p.stimrect = CenterRect(p.stimulussize,p.windowrect);
Screen('TextFont',p.window, 'Arial');
Screen('TextSize',p.window,p.fontsize);
Screen('BlendFunction',p.window,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

%fixation fuction
p.fixationlocation = .5*p.windowrect([3 3 3 3;4 4 4 4])+.5*[-p.fixationsize p.fixationsize 0 0;0 0 -p.fixationsize p.fixationsize];
drawfixation = @() Screen('DrawLines', p.window, p.fixationlocation,2,p.white);
drawfixationred = @() Screen('DrawLines', p.window, p.fixationlocation,2,[200 0 0]);
drawfixationgreen = @() Screen('DrawLines', p.window, p.fixationlocation,2,[0 200 0]);

%% disable mouse cursor
if p.isrealexperiment
    HideCursor();
    Priority(MaxPriority(p.window));
end

%% load stimuli into textures
tstart = GetSecs();
for i = 1:p.nstimuli
    p.stimuli(i).tex = Screen('MakeTexture',p.window,imread([p.stimuli(i).fol filesep p.stimuli(i).fn]));
    if ~mod(i,4)
        progress = ceil(100*i/p.nstimuli);
        eta = (p.nstimuli-i)*(GetSecs()-tstart)/i;
        DrawFormattedText(p.window,sprintf('Loading stimuli\n %i%% %.2fs',progress,eta),'center','center',p.white);
        Screen('Flip',p.window);
    end
end

%% Connect fmri
if p.isFMRIexperiment
    %DrawFormattedText(p.window,'Waiting for scanner','center','center',p.white);
    Screen('Flip',p.window);
    scansync('reset',p.tr);
end

%% ready to go
% Wait for key release
while KbCheck()
    WaitSecs(0.01);
end
keycode=[];
DrawFormattedText(p.window, '<SPACE> to start the experiment or <ESCAPE> to abort', 'center', 'center', p.white);
Screen('Flip',p.window);
while isempty(keycode) || ~keycode(KbName('space'))
    [~, keycode, ~] = KbWait();
    if keycode(KbName('escape'));eval(abort);end
end
    
%% we're go - record start time 
p.time_experiment_start = Screen('Flip', p.window); %record experiment start time

%% start experiment loop
nevents = size(eventlist,1);
nruns = eventlist.runnumber(end);
p.triggerlist = {};
currentrun = 0;
for eventnr=1:nevents
    isrunstart = eventlist.runnumber(eventnr)>currentrun;
    %things to do at the start of each run
    if isrunstart
        currentrun = eventlist.runnumber(eventnr);
        %save data
        writetable(eventlist,p.csvdatafilename)
        save(p.datafilename,'p');

        % Wait for key release
        while KbCheck()
            WaitSecs(0.01);
        end
        %instructions
        questiontext = sprintf('On which side are the spikes %s?', eventlist.question{eventnr});
        instructiontext = sprintf('Block %i/%i\n\n%s\n\nPlease keep still and wait for the experimenter to continue the experiment', currentrun, nruns, questiontext);
        if currentrun > 1 %show feedback if not first block
            idx = eventlist.runnumber==(currentrun-1);
            acc = mean(eventlist.correct(idx))*100;
            instructiontext = sprintf('%.2f%% correct\n\n\n%s',acc, instructiontext);
        end
        DrawFormattedText(p.window, instructiontext, 'center', 'center', p.white);
        Screen('Flip', p.window);
        if p.isFMRIexperiment && currentrun > 1 % wait to get full HRF of last stim
            WaitSecs(p.feedbackduration);
        end
        %wait for any keypress to start the run
        keycode=[];
        while isempty(keycode) || ~keycode(KbName('space'))
            [~, keycode, ~] = KbWait();
            if keycode(KbName('escape'));eval(abort);end
        end

        % Get n triggers as a baseline
        triggerlist = nan(p.num_baseline_triggers,2);
        for trigger_num = 1:p.num_baseline_triggers
            instructiontext = sprintf('Block %i/%i\n\n%s\n\nWaiting for scanner', currentrun, nruns, questiontext);
            DrawFormattedText(p.window, instructiontext, 'center', 'center', p.white);
            Screen('Flip', p.window);
            if p.isFMRIexperiment
                % Get trigger
                [pulse_time,~,daqstate] = scansync(1,Inf);
                triggerlist(trigger_num,:) = [pulse_time daqstate.nrecorded(1,1)];
            else
                if p.this_is_not_a_drill %wait for scanner
                    WaitSecs(p.tr);
                end
                triggerlist(trigger_num,1) = GetSecs();
            end
        end
        time_runstart = triggerlist(1,1);
        time_trigger_mean = mean(diff(triggerlist(:,1)));
        time_trigger_std = std(diff(triggerlist(:,1)));
        p.triggerlist{currentrun} = triggerlist;
        if p.this_is_not_a_drill %wait for scanner
            drawfixation();
            time_runfixstart = Screen('Flip', p.window, time_runstart + p.block_num_message_time - p.halfrefresh);
        else
            drawfixation();
            time_runfixstart = Screen('Flip', p.window);
        end
        drawfixation();
        Screen('Flip', p.window, time_runfixstart + p.blockstartfixationduration - p.halfrefresh);
    end
    
    %draw stim    
    Screen('DrawTexture',p.window,p.stimuli(eventlist.stimnumber(eventnr)).tex,[],p.stimrect);
    DrawFormattedText(p.window, questiontext, 'center', .75*p.windowrect(4), p.white);
    time_stimon = Screen('Flip',p.window);
    
    drawfixation();
    DrawFormattedText(p.window, questiontext, 'center', .75*p.windowrect(4), p.white);
    time_stimoff = Screen('Flip', p.window, time_stimon + p.stimulusduration - p.halfrefresh);
    
    %check for reponse
    [response,RT,correct]=deal(0);
    if p.isFMRIexperiment
        % Wait the allowed time for a button to be pressed
        rt = scansync([2,3],time_stimon+p.responseduration);
        % Get the response button and rt
        if any(isfinite(rt))
            [RT,resp] = min(rt);
            RT = RT-time_stimon; % Subtract stimon time to get the RT
            if resp==1
                response = 'left';
            elseif resp==2
                response = 'right';
            end
            correct = strcmp(eventlist.correctresp(eventnr),response);
        end
        % for practice, show feedback
        if ~p.this_is_not_a_drill
            if correct
                drawfixationgreen();
            else
                drawfixationred();
            end
            DrawFormattedText(p.window, questiontext, 'center', .75*p.windowrect(4), p.white);
            Screen('Flip',p.window);
        end            
        % Use remaining time to check for an escape key to allow exit
        while GetSecs()-time_stimon < p.responseduration
            [keydown,~,keycode] = KbCheck();
            if keydown && keycode(KbName('escape')) % emergency break
                eval(abort);
            end
        end
    else %not fmri - for debugging or training
        while GetSecs()-time_stimon < p.responseduration
            [keydown, secs, keycode] = KbCheck();
            if keydown
                if ~response
                    if any(keycode(p.responsekeysleft))
                        response = 'left';
                    elseif any(keycode(p.responsekeysright))
                        response = 'right';
                    end
                    if response
                        RT = secs-time_stimon;
                        correct = strcmp(eventlist.correctresp(eventnr),response);
                        if ~p.this_is_not_a_drill
                            if correct
                                drawfixationgreen();
                            else
                                drawfixationred();
                            end
                            DrawFormattedText(p.window, questiontext, 'center', .75*p.windowrect(4), p.white);
                            Screen('Flip',p.window);
                        end          
                    end
                end
                if keycode(KbName('escape')) %emergency exit
                    eval(abort);
                end
            end
        end
    end
    %update the eventlist with results
    eventlist.response{eventnr} = response;
    eventlist.RT(eventnr) = RT;
    eventlist.correct(eventnr) = correct;
    eventlist.time_stimon(eventnr) = time_stimon-time_runstart;
    eventlist.time_stimoff(eventnr) = time_stimoff-time_runstart;
    eventlist.stimdur(eventnr) = time_stimoff-time_stimon;
    eventlist.time_stimon_abs(eventnr) = time_stimon-p.time_experiment_start;
    eventlist.time_stimoff_abs(eventnr) = time_stimoff-p.time_experiment_start;
    eventlist.time_runstart(eventnr) = time_runstart-p.time_experiment_start;
    eventlist.time_trigger_mean{eventnr} = time_trigger_mean;
    eventlist.time_trigger_std{eventnr} = time_trigger_std;
end
finish up 
writetable(eventlist,p.csvdatafilename)
save(p.datafilename,'p','eventlist');

acc = mean(eventlist.correct);
DrawFormattedText(p.window, sprintf('%.2f%% correct\n\n\nExperiment complete!\n\nRelax and wait for experimenter...\n\nExperimenter: press <space> to exit',acc), 'center', 'center', p.white);

Screen('Flip', p.window);
keycode=[];
while isempty(keycode) || ~keycode(KbName('space'))
    [~, keycode, ~] = KbWait();
end
Priority(0);ListenChar(0);ShowCursor();
Screen('CloseAll');



