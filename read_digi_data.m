function out = read_digi_data(dataset, sname, fictrac_dir, ~, pred, crab, ball, ori, imcorrect, claw_pos)

%  if nargin < 2; sname = 'G1T1_2'; end
%  if nargin < 1 || isempty(dataset); dataset = 'F:\Callum\Digi'; end

pred = [1 2] ;
ori = [];
crab = 4;
ball = 3;
imcorrect = 0;
claw_pos = [];
make_defence = 0;

% Reads files with .data extension from the resfile folder in the directory
% provided to the function
% dataset= 'F:\juan\data\behavioural_data\different_parameters';
% sname= 'crab10_protocol1_06242022123720_0000';
fname = fullfile(dataset, 'resfiles', sname, [sname '.data']);
% Prints the current dataset of the loop - check that all have been read
fprintf('digi %s %s loadfull\n', dataset, sname)

% will load save_data
load(fname, '-mat');


if ~isempty(fictrac_dir)
  if contains(fname, '_combined.data') % find the fictrac names for combined digi files
    [fd fn fe] = fileparts(fname);
    fname
    a = dir(fullfile(fictrac_dir, [fn(1:4) '*.dat']));
    for fi = 1:numel(a)
      fic_name = fullfile(a(fi).folder, a(fi).name);
      protocol = str2double(a(fi).name(6:7));
      fic_data(protocol).dat = dlmread(fic_name);
    end      
  else
    a = dir(fullfile(fictrac_dir, [sname '*.dat']));
    if numel(a)==1
      fic_name = fullfile(a(1).folder, a(1).name);
    else
      fullfile(fictrac_dir, [sname '*.dat'])
      error('ambiguous or empty fic_track file');
    end
    fic_data(1).dat = dlmread(fic_name);
  end
end


% Defines the markers used
frametext = save_data.markers.frametext;
mmarker = save_data.markers.marker;
runmarker = save_data.markers.runmarker;
pmarker = save_data.markers.pointmarker;


% Need to define which runmarker corresponds to the different states
% present during a run (i.e. run start and end, fade start and end)

% Finds the runmarker identifiers
if isempty(runmarker)
  error('no runmarker found, fix input file');
end
temp = unique(runmarker(:,2));

% If statement checks if a fade was used or not. If so, runmarker 2 and 3
% correspond to the start and end of fade. If not, runmarker 2 and 3
% correspond the the start and end of the predator stimulus
if numel(temp) == 5
  runs = runmarker(runmarker(:,2) == temp(1),:); % start of run
  runas = runmarker(runmarker(:,2)== temp(2),:); % fade appear
  runae = runmarker(runmarker(:,2)== temp(3),:); % fade end
  runss = runmarker(runmarker(:,2)== temp(4),:); % stimulus start
  runse = runmarker(runmarker(:,2)== temp(5),:); % stimulus ends
  rune =runse;    %runend
elseif numel(temp) == 3
  runs = runmarker(runmarker(:,2) == temp(1),:); % start of run
  runss = runmarker(runmarker(:,2)== temp(2),:); % stimulus start
  runse = runmarker(runmarker(:,2)== temp(3),:); % stimulus end
  rune =runse;    %runend
  videos = 0;
elseif numel(temp) == 4
  runs = runmarker(runmarker(:,2) == temp(1),:); % start of run
  runss = runmarker(runmarker(:,2)== temp(2),:); % stimulus start
  runse = runmarker(runmarker(:,2)== temp(3),:); % stimulus end
  videos = runmarker(runmarker(:,2)== temp(4),:); % video start
  rune =runse;    %runend
end
% prints the digi file where the number of run starts and run ends are not
% balanced. The printed directory needs to be balanced
if size(runs,1) ~=size(rune,1)
  fprintf('digi %s %s  not balanced\n', dataset, sname)
  pause
end




% This loop builds the output structure for all relevant information in
% each dataset file

for i = size(runs,1):-1:1
  %     [i runs(i,1)   rune(i,1)]
  
  temp_text = frametext.text{i};
  sel = strfind(temp_text, ',');
  if isempty(sel)
    sel = strfind(temp_text , ' ');
  end
  
  
  
  % Finds the number before the delimiter and places them in seperate
  % vectors. Using the ':' to from the logical position until the sel
  % position allows the number to be greater than one digit
  out(i).crab = str2double(temp_text(1:sel(1)-1));
  
  if out(i).crab <= inf % here you can make two groups
    out(i).group= 1;
    out(i).crab_in_group= out(i).crab;
  else
    out(i).group= 2;
    out(i).crab_in_group= out(i).crab-10;
  end
  
  % out(i).group = str2double(temp_text(sel(1)+1: sel(2)-1));
  out(i).csize = str2double(temp_text(sel(2)+1: sel(3)-1));
  out(i).trial = str2double(temp_text(sel(3)+1:sel(4)-1));
  out(i).sex = str2double(temp_text(sel(4)+1:end));
  
  
  
  if ~isempty(crab)
    % the frame number and behaviour associated with each pointmarker
    temp = pmarker{crab};
    
    % Saves point data for the ball centre and predator direction. This
    % assumes that the ball and predator was marked on the first runmarker.
    % If not, the obeject runs need to change to reflect where it was
    % marked. The +1 is to shift the frame by one as digi is zero based.
    out(i).ballcentre = save_data.data(ball).pos(runs(i)+1,1:2);
    out(i).rundir1=nan;
    
    
    
    if imcorrect > 0
      if out(i).group == 1 || out(i).group == 2
        imdistort = 0.75;
      else
        imdistort = 1;
      end
    elseif imcorrect == 0
      imdistort = 1;
    end
    
    out(i).rundir_digi = [nan nan nan nan nan]; % allow for 5 runs
    rno=1;
    
    if ~isempty(temp)
      sel_run = temp(:,1) >= (runss(i,1)) & temp(:,1) <=(runse(i,1));
      behav = temp(sel_run,:);
      
      for ii=1:size(behav,1)
        % Only calculates if the behaviour is a run and occurs between the
        % stimulus start and end. delay accounts for the neural lag
        % if rno is equal to 1, then it checks for a run behav and does the
        % calc, else if it is greater than 1, it checks for run behav or
        % direction change behav
        if rno == 1
          if behav(ii,2)=='r'
            % Position of the run behaviour (x,y co-ord)
            rpos=save_data.data(crab).pos(behav(ii,1)+1,1:2);
            % Difference between the pos of the run behaviour of the crab and the
            % centre of the treadmill ball
            rv = rpos-out(i).ballcentre;
            out(i).rundir_digi(1, rno) = atan2(-rv(2)*(imdistort), rv(1))*180/pi; % running direction (deg)
            rno=rno+1;
            
          end
        elseif rno > 1
          if behav(ii,2)=='r' || behav(ii,2)=='d'
            rpos=save_data.data(crab).pos(behav(ii,1)+1,1:2);
            % Difference between the pos of the run behaviour of the crab and the
            % centre of the treadmill ball
            rv = rpos-out(i).ballcentre;
            out(i).rundir_digi(1, rno) = atan2(-rv(2)*(imdistort), rv(1))*180/pi; % running direction (deg)
            rno=rno+1;
          end
        end
      end
    end
  end
  
  
  if ~isempty(claw_pos)
    
    claw_position = [NaN NaN];
    
    for ii=1:size(behav,1)
      if behav(ii,2)=='r' && out(i).sex == 1
        claw_position =  save_data.data(claw_pos).pos(behav(ii,1)+1,1:2);
      end
    end
    out(i).claw_position = claw_position;
    
  elseif isempty(claw_pos)
    
    out(i).claw_position = NaN;
    
  end
  
  % line specific

  if ~isempty(fictrac_dir)&& ~isempty(fic_data(i).dat)
    run_track = [NaN NaN];
    DirVectorY = NaN;
    DirVectorX = NaN;
    sel2 = [];
    sel3 = [];
    
    sel = find(behav(:,2) == 'r', 1,'first');
    if (sel + 1) > size(behav,1)
      rune_ind = rune(i);
    else
      rune_ind = behav(sel+1,1);
    end
    run_ind = behav(sel,1);
    if ~isempty(run_ind)
      sel2 = find(run_ind-videos(i) == fic_data(i).dat(:,1), 1);
      sel3 = find(rune_ind-videos(i) == fic_data(i).dat(:,1), 1);
      %         [sel2 sel3 run_ind rune_ind size(run_track) size(fic_data(i).dat)]
      
      run_track(1:numel(sel2:sel3),1:2) = fic_data(i).dat(sel2:sel3,[21 20]); % x y integrated
      %             run_track(:,1) = cumsum(fic_data(i).dat(sel2:sel3,6))*-1; % x integrated
      %             run_track(:,2) = cumsum(fic_data(i).dat(sel2:sel3,7)); % y integrated
      try
        DirVector = run_track(end,:) - run_track(1,:);
      catch ME
        ME
        keyboard
      end
      
      run_dir = atan2(DirVector(2), DirVector(1));
      
      for jj = 1:size(run_track(:,1),1)-1
        
        DirVectorX = run_track(jj+1,1) - run_track(jj,1);
        DirVectorY = run_track(jj+1,2) - run_track(jj,2);
        
        run_all(jj) = atan2(DirVectorY, DirVectorX);
      end
      %             run_dir = mod(run_dir-pi/2, pi*2); % rotated by 90 degrees to fit digi
      
      
    else
      run_track = NaN(1,2);
      run_dir = NaN;
      run_all = NaN;
      
    end
    
    % x and y co-ordinates of run track (not in deg)
    out(i).run_track = run_track;
    
    % Vector between first and last run behaviour (mean direction)
    out(i).dir_track = run_dir;
    
    % Each small vector (i.e., all vectors between all points)
    out(i).run_all = run_all;
    
    out(i).run_speed = fic_data(i).dat(sel2:sel3,19);
  else
        % x and y co-ordinates of run track (not in deg)
    out(i).run_track = NaN;
    
    % Vector between first and last run behaviour (mean direction)
    out(i).dir_track = NaN;
    
    % Each small vector (i.e., all vectors between all points)
    out(i).run_all = NaN;
    
    out(i).run_speed = NaN;
  end
  
  
  % Calculates the position of the predator
  if    ~isempty(ball)
    for ii = 1:numel(pred)
      out(i).pred_dat(ii,1:2) = save_data.data(pred(ii)).pos(runs(i)+1, 1:2);
      
      % This loop uses trig between the ball centre and the digitized predator approach
      % position
      if sum(out(i).pred_dat(ii,1:2)) ~= 0
        ttemp = out(i).pred_dat(ii,1:2)-out(i).ballcentre;
        %         out(i).pred_a(ii,:) = atan2(abs(ttemp(2))*(imdistort), abs(ttemp(1)))*180/pi;
        %             out(i).pred_a(ii,:) = atan2((ttemp(2)*imdistort), ttemp(1))*180/pi;
        out(i).pred_a(ii,:) = atan2((ttemp(2)*imdistort), ttemp(1));
      else
        out(i).pred_a(ii,:) = NaN;
      end
      % Finds the direction of the approaching predator. Currently in x,y
      % co-ordindates
%       out(i).dir(ii,:) = save_data.data(pred(ii)).pos(runs(i,1)+1, 1:2);
      
      
      % find which predator the large claw is closest to
      %             claw_x = out(i).dir(ii,1) - out(i).claw_position(1);
      %             claw_y = out(i).dir(ii,2) - out(i).claw_position(2);
      %             out(i).claw_angle(ii) = atan2(claw_y, claw_x);
      
      if ~isempty(claw_pos)
        claw_x = out(i).claw_position(1) - out(i).ballcentre(1);
        claw_y = out(i).claw_position(2) - out(i).ballcentre(2);
        out(i).claw_angle(ii) = atan2(claw_y, claw_x);
      elseif isempty(claw_pos)
        out(i).claw_angle(ii) = NaN;
      end
      
    end
    
  end
  
  
  
  
  
  %behaviours % adjust selection for reaction time
  %    out(i).pm = temp(temp(:,1)>=(runss(i,1)+delay*60) & temp(:,1)<=(runse(i,1)+delay*60), :);
  
  % all the behaviours between the start and end of the stimulus
  if ~isempty(temp)
    out(i).pm = temp(temp(:,1)>=(runss(i,1)) & temp(:,1)<=(runse(i,1)), :);
    
    % All behaviours from the whole run
    out(i).pm_all = temp(temp(:,1) <= (runse(i,1)),:);
    
    % Number of frames past the start of stimulus when that behaviour
    % was observed
    out(i).pm(:,1) = out(i).pm(:,1) - runss(i, 1);
    % Number of frames before (-) and after the strat of stimulus in
    % which that behaviour was observed
    out(i).pm_all(:,1) = out(i).pm_all(:,1) - runss(i, 1);
  else
    out(i).pm =[];
    out(i).pm_all =[];
  end
  
  
  % Frametext - should indicate excess information e.g. sex and crad i.d.
  temp = frametext;
  out(i).ft = temp.text{i};
  
  % This is the marker marker which represent condition (Important!)
  % subset for the marker present at the frame of run start
  temp = mmarker;
  
  if sum(temp(:,1)==runs(i,1))==0
    fprintf('marker marker missing\n');
    keyboard
  end
  out(i).mm = temp(temp(:,1)==runs(i,1), 2);
  
  % Finds the row number where a behaviour was marked and concatinates the
  % row number with the x,y co-ordindates of the behaviour, i.e. the
  % running/walking direction (only after the stimulus starts)
  temp = save_data.data(crab).pos(runss(i,1)+1:runse(i,1)+1, 1:3);
  sel = find(sum(temp>0, 2)~=0);
  temp = [sel temp(sel,1:2)];
  
  % The same thing but for the entire run i.e. before the stimulus starts
  temp2 = save_data.data(crab).pos(runs(i,1)+1:rune(i,1)+1, 1:3);
  sel = find(sum(temp2>0, 2)~=0);
  temp2 = [sel temp2(sel,1:2)];
  
  % Loads the directions into the structure
  out(i).rdir = temp;
  out(i).rdir_full = temp2;
  
  
  % Put here for tidyness. Each ori point contains x and y data and
  % therefore has a total of 8 data points. This defines if you want both
  % x and y data (1:2) or just x data (1) or just y data (2)
  
  % Selects data between run marker 2 and 3 (start and end of stim)
  sel = [];
  temp = [];
  temp2 = [];
  nruns = 3;
  
  temp2 = cell(nruns, numel(ori));
  for ii = 1:numel(temp2)
    temp2{ii} = [NaN NaN];
  end
  
  sel = runss(i,1):runse(i,1)+1;
  
  % check to see if user has oritation argument
  if ~isempty(ori)
    
    for ii = 1:numel(ori)
      
      temp = save_data.data(ori(ii)).pos(sel, 1:2);
      sel_ori = find(temp(:,1) > 0);
      
      if ~isempty(sel_ori)
        for j = 1:numel(sel_ori)
          
          temp2{j, ii} = temp(sel_ori(j),:);
          
        end
      end
    end
    
    
    %         if out(i).trial == 3 || out(i).trial == 4
    %
    %             out(i).crab
    %         end
    temp2 = cell2mat(temp2)';
    
    crab_ori = reshape(temp2, 1, size(temp2,1), size(temp2,2));
    ori_xy = [(crab_ori(:,2,1) - crab_ori(:,4,1)) (crab_ori(:,1,1) - crab_ori(:,3,1))];
    out(i).crab_ori = atan2(ori_xy(1), -ori_xy(2));
    
    if sum(isnan(crab_ori(:,:,nruns))) ~= 0
      out(i).crab_ori_middle = NaN;
      
      ori_xy_end = [(crab_ori(:,2,2) - crab_ori(:,4,2)) (crab_ori(:,1,2) - crab_ori(:,3,2))];
      out(i).crab_ori_end = atan2(ori_xy_end(1), -ori_xy_end(2));
      
    elseif sum(isnan(crab_ori(:,:,nruns))) == 0
      ori_xy_middle = [(crab_ori(:,2,2) - crab_ori(:,4,2)) (crab_ori(:,1,2) - crab_ori(:,3,2))];
      out(i).crab_ori_middle = atan2(ori_xy_middle(1), -ori_xy_middle(2));
      
      ori_xy_end = [(crab_ori(:,2,nruns) - crab_ori(:,4,nruns)) (crab_ori(:,1,nruns) - crab_ori(:,3,nruns))];
      out(i).crab_ori_end = atan2(ori_xy_end(1), -ori_xy_end(2));
    end
  end
  
  
  
  
  % Loads all runmarkers into the structure
  out(i).rm = runmarker(runmarker(:,1)>=runs(i,1) & runmarker(:,1)<=rune(i,1),:);
  out(i).rm(:,1) = out(i).rm(:,1) - runss(i, 1);
  %   out(i).runas = runas(i,1);
  out(i).runs = runs(i,1);
  out(i).rune = rune(i,1);
  out(i).runse = runse(i,1);
  out(i).runss = runss(i,1);
  %   out(i).runae = runae(i,1);
  out(i).sname = sname;
  out(i).dataset = dataset;
  out(i).imagesize = save_data.imagesize;
  
end

end

