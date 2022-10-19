function make_excel_behaviour(path,file)
% this function makes an excel spreadsheet from the res structure of a behavioural experiment so
% that you can later analyse (and plot) in R.
% the spreadsheet is saved as 'res.xlsx' in the same folder as the res.mat file
% so far it only works for a single predator, will update later.
% so far it only works for freeze, run and limb tuck.

% load data
% path= 'C:\data\juan\MEGAsync\data\behavioural_data\different_parameters\';
% file= 'res.mat';

load(strcat(path,file));

% make the excel spreadsheet with the desired variables
% I only selected the few I was interested in
% first the column names

names={'protocol','crab','sex','trial','csize','group','crab_in_group','rundir','pred_angle','pred_pos',...
  'psize','pspeed','sasize','ftime','rtime','ltime','ft2c','rt2c','lt2c','fdist','rdist',...
  'ldist','fasize','rasize','lasize','fespeed','respeed','lespeed','fanginc','ranginc','langinc'};

data= cell(length(res.crab)+1,length(names));
data(1,:)= names;
data(2:end,:)= num2cell([res.mm,res.crab,res.sex,res.trial,res.csize,res.group,res.crab_in_group,...
  res.rundir(:,1),res.pred_angle,res.pred_pos,res.size,res.speed,res.start_asize,...
  res.response_time(:,1),res.response_time(:,2),res.response_time(:,3),...
  res.t2c(:,1),res.t2c(:,2),res.t2c(:,3),...
  res.dist(:,1),res.dist(:,2),res.dist(:,3),...
  res.asize(:,1),res.asize(:,2),res.asize(:,3),...
  res.espeed(:,1),res.espeed(:,2),res.espeed(:,3),...
  res.anginc(:,1),res.anginc(:,2),res.anginc(:,3)]);

xlswrite(strcat(path,'res.xlsx'), data)