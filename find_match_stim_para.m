% rm2 = runmarker 2 from digi
% rm3 = runmarker 3 from digi
% cframe = frame from digi (-rm2) you want the values from
% s are the input parameters used for the stimulus creation
% out = structure with 6 values (for ref and matching stimulus)

function out = find_match_stim_para(s, rm2, rm3, cframe)


% s.distance=5000; %initial distance of the reference stimulus in mm
% s.speed=200;%the speed of the reference stimulus in mm/s
% s.equal_distance=445.04;%the distance that both stimuli would be there at the exact0 same time
% s.frameRate=60;%monitors frame rate
% s.dr=22;%diameter of the reference stimulus in mm
% s.d=45;%diameter of the matching stimulus in mm
% s.elevation=15; %elevation of stimulus at initial position in deg
% s.XYangle=45;%Angle of end position in XY plane relative to X axis, you only need this if you have a near miss stimulus
% s.finalDistance=1.5;%final distance of stimulus from the observer
% % s.passingDistance=0; 
% s.traveling_distance=5000-s.finalDistance;
% s.travel_time = s.traveling_distance/s.speed;


noplot = 1;
s = matchingPureLoomingExpansionSpeed(s, noplot);
ldat =length(s.AsizeMatching_deg); % length of calculated entry (1501)
cfframes = cframe/(rm3 - rm2)*(ldat-1)+1;
 
out.asize_match = linterp(1:ldat, s.AsizeMatching_deg, cfframes);
out.asize_ref = linterp(1:ldat, s.AsizeReference_deg, cfframes);
out.dist_match = linterp(1:ldat, s.MatchingDistances, cfframes);
out.dist_ref = linterp(1:ldat, s.Distances_reference, cfframes);
out.aspeed_match = linterp(1:ldat, s.AspeedMatching, cfframes);
out.aspeed_ref = linterp(1:ldat, s.AspeedReference, cfframes);


