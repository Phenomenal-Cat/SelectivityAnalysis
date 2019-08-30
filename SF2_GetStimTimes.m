function Stim = SF2_GetStimTimes(TDTdata, NoStim)

%========================= MF3D2_GetSTimTimes.m ===========================
% This function retreives the timing and stimulus identity information for
% every stimulus presentation that occured within the TDT block specified
% by the input directory.

persistent PD

if ischar(TDTdata) & exist(TDTdata, 'file')
    load(TDTdata);
end

PlotAnalogData = 1;
EventCodes  = TDTdata.epocs.PtAB.data;
EventTimes  = TDTdata.epocs.PtAB.onset;
PDdata      = TDTdata.streams.x_Pd.data;
PD_srate    = TDTdata.streams.x_Pd.fs;
PDTimes     = linspace(0, (1/PD_srate)*numel(PDdata), numel(PDdata));

PDthresh    = 2.5;
if ~exist('PD','var')
    PD          = CheckPhotodiode_v2(PDdata, PD_srate, [], PDTimes(end));      	% Extract photodiode state changes from signal
else
    PD          = CheckPhotodiode_v2(PDdata, PD_srate, PD, PDTimes(end));
end

%================= Optional data plotting
if PlotAnalogData == 1
    AllFields   = fieldnames(TDTdata.streams);
    AllFields   = AllFields(~cellfun(@isempty, strfind(AllFields, 'x_')));
    fh      = figure('units','normalized');
    axh     = tight_subplot(numel(AllFields), 2, 0.02, 0.02, 0.02);
    axh     = reshape(axh, [numel(AllFields), 2]);
    Ypos    = linspace(0.88, 0.05, numel(AllFields));
    for n = 1:numel(AllFields)
        axes(axh(n,1));
        data = eval(sprintf('TDTdata.streams.%s.data', AllFields{n}));
        plot(PDTimes, data);
        axis tight
        box off
        grid on
        title(AllFields{n}(3:end));
        axes(axh(n,2));
        hist(data, 1000);
        axis tight
        box off
        grid on;
        set(axh(n,1),'position', [0.05, Ypos(n), 0.65, 0.1]);
        set(axh(n,2),'position', [0.75, Ypos(n), 0.24, 0.1],'yticklabel',[]);
    end
    set(axh(1:n-1,:),'xticklabel',[], 'tickdir','out');
    linkaxes(axh(:,1),'x')
end

%================= Get eye calibration offset & gain values
Events          = SCNI_LoadEventCodes();
E_OffsetIndx    = find(~cellfun(@isempty, strfind({Events.String}, 'Sending_EyeOffsets')));
E_GainIndx      = find(~cellfun(@isempty, strfind({Events.String}, 'Sending_EyeGains')));
OffsetIndx      = find(EventCodes==Events(E_OffsetIndx).TDTnumber);
GainIndx        = find(EventCodes==Events(E_GainIndx).TDTnumber);
for n = 1:numel(OffsetIndx)
    OffsetValues(n,:)   = (EventCodes(OffsetIndx(n)+[2,4])-4000)/1000;
    Gainvalues(n,:)     = (EventCodes(GainIndx(n)+[2,4])-4000)/1000;
end

%================= Find stimulus onset time and ID# for each presentation
% StimIndices         = unique(EventCodes(EventCodes<2000 & EventCodes >0));
% NoStim              = numel(StimIndices);
StimOnVal           = Events(find(~cellfun(@isempty, strfind({Events.String}, 'Stim_On')))).TDTnumber;
StimOffVal          = Events(find(~cellfun(@isempty, strfind({Events.String}, 'Stim_Off')))).TDTnumber;
StimOnIndx          = find(EventCodes==StimOnVal);
StimOnTime          = EventTimes(StimOnIndx);

if ~exist('Stim','var') || isempty(Stim)
    Stim.Count  	= zeros(1,NoStim);
end
if StimOnIndx(end) == numel(EventCodes)
    StimOnIndx(end) = [];
end

Stim.Times = cell(1,NoStim);
for s = 1:numel(StimOnIndx)                                     % For each 'StimOn' event code...
    StimOffTemp  = StimOnIndx(s) + find(EventCodes((StimOnIndx(s)+1):end)==StimOffVal, 1, 'first');
    if ~isempty(StimOffTemp)
        StimOffIndx(s)  = StimOffTemp;
        StimOffTime(s)  = EventTimes(StimOffIndx(s));
        DistTemp       	= find(EventCodes((StimOnIndx(s)+1):end)>0 & EventCodes((StimOnIndx(s)+1):end)<=NoStim, 1, 'first');
        if ~isempty(DistTemp)
            Dist(s)     = DistTemp;
            StimID(s)   = EventCodes(StimOnIndx(s) + Dist(s));
            PDonIndx(s) = find(PD.OnTimes > StimOnTime(s),1,'first');
            PDoffIndx(s)= find(PD.OffTimes > StimOffTime(s),1,'first');
            PDdur(s)    = PD.OffTimes(PDoffIndx(s))-PD.OnTimes(PDonIndx(s));
            TimeDiff(s) = PD.OnTimes(PDonIndx(s))-StimOnTime(s);
            if PDdur(s) > 0.2
                Stim.Count(StimID(s))                           = Stim.Count(StimID(s))+1;
                Stim.Times{StimID(s)}(Stim.Count(StimID(s)))    = PD.OnTimes(PDonIndx(s));
            end
        end
    end
end
%StimOnIndx(Dist > 2)    = [];          % Remove 'stim on' if stim ID# was not found


% StimOnTimes         = EventTimes(StimOnIndx);
% StimIndices         = EventCodes(StimOnIndx+1);
% StimOffsetEvntTimes = EventTimes(EventCodes==StimOffVal);

%T = readtable('/projects/murphya/Stimuli/AvatarRenders_2018/StereoShape/StereoShape_StimSummary.csv');