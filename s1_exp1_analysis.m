%% Read in data

dataset = 'E:\Helasha\EX1';
% dataset = 'F:\juan\data\behavioural_data\different_parameters';

delay = 0.14; % neural delay acc. to Callums paper
framerate = -1; 

% Which column of the crab markers in digi was used to mark behaviours
% (crab), predator approach direction (pred) and the centre of the
% treadmill ball (ball)
pred = [1 2];
ori = []; % ORIENTATIONB OF the CRAB
crab = 4;
ball = 3;
imcorrect = 0;
claw_pos = [];
make_defence = 0;


fictrac_dir = 'E:\Helasha\FicTrac_analysis\fictrac_output';
% fictrac_dir = 'C:\data\fictrac_run\Chapter3\dat';
% fictrac_dir2 = [];

redo_dl = 0;
redo_dt = 1;

if redo_dl == 1
  cdata = dataset_loop(dataset, fictrac_dir, delay, pred, crab, ball, ori, imcorrect, claw_pos);
  save('./dres.mat', 'cdata')
else
  load('./dres.mat', 'cdata');
end  

if redo_dt == 1
  res = data_trans(cdata, delay, framerate, make_defence);
  save('./res.mat', 'res')
else  
  load('./res.mat', 'res');
end  


% converts angles to 0 - 2*pi
running = mod(res.dir_track, 2*pi);





%% Circular data

% Place larger predator at the bottom (270 deg)
% run_rot = Inf(numel(running),1);
% ori_rot = Inf(numel(running),1);
% all_dir_rot = Inf(size(all_dir2,1), size(all_dir2,2));
% 
% run_rot2 = NaN(numel(running),1);
% run_rot2(run_rot2 > pi) = run_rot2(run_rot2 > pi) - (2*pi);

% ori_rot2 = NaN(numel(running),1);
% % ori_rot2(ori_rot2 > pi) = ori_rot2(ori_rot2 > pi) - (2*pi);
% 
% all_dir_rot2 = all_dir2;
% all_dir_rot2(all_dir_rot2 > pi) = all_dir_rot2(all_dir_rot2 > pi) - (2*pi);


%% 180 degree seperation

sel9 = res.pred_pos(:,1) == 1;
run_rot(sel9) = running(sel9) - 0;
run_rot2(sel9) = running(sel9) - 0;
% ori_rot(sel9) = crab_ori(sel9) - 0;
% all_dir_rot(:,sel9) = all_dir2(:,sel9) -0;
% ori_rot2(sel9)  = crab_ori(sel9) - 0;

% pred_pos(sel9) = 0;

sel10 = res.pred_pos(:,1) == 3;
run_rot(sel10) = running(sel10) + pi;
run_rot2(sel10) = running(sel10) + pi;
% ori_rot(sel10) = crab_ori(sel10) + pi;
% all_dir_rot(:,sel10) = all_dir2(:,sel10) + pi;
% ori_rot2(sel10) = crab_ori(sel10) + pi;

% pred_pos(sel10) = 0;

sel11 = res.pred_pos(:,1) == 2;
run_rot(sel11) = running(sel11) - (pi/2);
run_rot2(sel11) = running(sel11) - (pi/2);
% ori_rot(sel11) = crab_ori(sel11) + (pi/2);
% all_dir_rot(:,sel11) = all_dir2(:,sel11) + (pi/2);
% ori_rot2(sel11) = crab_ori(sel11) + (pi/2);

% pred_pos(sel11) = 0;

sel12 = res.pred_pos(:,1) == 4 ;
run_rot(sel12) = running(sel12) + (pi/2);
run_rot2(sel12) = running(sel12) + (pi/2);
% ori_rot(sel12) = crab_ori(sel12) - (pi/2);
% all_dir_rot(:,sel12) = all_dir2(:,sel12) - (pi/2);
% ori_rot2(sel12)  = crab_ori(sel12) - (pi/2);

% pred_pos(sel12) = 0;


run_rot = mod(run_rot, (2*pi))*180/pi;
run_rot2 = mod(run_rot2, (2*pi))*180/pi;
% ori_rot = mod(ori_rot, (2*pi));
% ori_rot2 = mod(ori_rot2, (2*pi));
% all_dir_rot =  mod(all_dir_rot, (2*pi));
run_rot_diff = run_rot;

sel = res.mm == 2;
subplot(2,2,1)
ed = -7.5:15:367.5;
a=histc(running(sel)*180/pi, ed)
bar(ed, a, 'histc')
set(gca, 'xlim', [0 360])

subplot(2,2,2)
plot(res.pred_pos(sel,1), 'ko')
subplot(2,2,3)
ed = -7.5:15:367.5;
a=histc(run_rot(sel), ed)
bar(ed, a, 'histc')
set(gca, 'xlim', [0 360])
subplot(2,2,4)
hist(res.pred_pos)





















