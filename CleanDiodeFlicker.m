function NewSignal = CleanDiodeFlicker(Signal)

%======================== CleanDiodeFlicker.m =============================
% Takes a phodiode signal that has been poluted by picking up screen
% refresh (typically 60Hz) and removes the pollution, leaving only genuine
% state transitions.
%

Thresh          = 1;
AboveThresh     = Signal>Thresh;                              
Onsets          = find(diff(AboveThresh)==1)+1;
Offsets         = find(diff(AboveThresh)==-1)+1;
Onsets(1)       = [];
Offsets(end)    = [];
OffDurs         = Onsets-Offsets;
ValidOffDurs    = find(OffDurs>20);
NewSignal       = zeros(size(Signal));
for n = 1:numel(ValidOffDurs)
    NewSignal(Offsets(ValidOffDurs(n)):Onsets(ValidOffDurs(n))) = 1;
end
