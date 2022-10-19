function res = data_trans(cdata, delay, framerate, make_defence)

% Input camera framerate
res.framerate = framerate;
res.frametime= 1000/res.framerate;

%turn runmarker data into a matrix
for i = 1:numel(cdata)
  res.rm(i,:) = cdata(i).rm(:,1);
end

%take data out of stucture array and turn into a matrix
res.mm = [cdata.mm]';
% res.dir = [cdata.dir]';
res.crab = [cdata.crab]';
res.sex = [cdata.sex]';
res.trial = [cdata.trial]';
res.csize = [cdata.csize]';
res.group = [cdata.group]';
res.runss = [cdata.runss]';
res.runse = [cdata.runse]';
res.dir_track = [cdata.dir_track]';
res.crab_in_group = [cdata.crab_in_group]';

temp = reshape([cdata.rundir_digi],5,[])';
res.rundir = temp;


% res.runddir_mean =
% if size(temp,2)>2
%     warning('you have more than 2 runs for a crab');
% end

res.pred_angle = (rad2deg([cdata.pred_a]'));

% Create rule to determine the position of predator (i.e., is it apporaching from mon 1, 2 or 3)
% res.pred_pos= res.pred_angle >= 45;
% res.pred_dir = vertcat(cdata.pred);
% res.pred_dir = res.pred_dir(:,1);


res.pred_pos= zeros(length(res.pred_angle), size(res.pred_angle, 2));
for j= 1:size(res.pred_angle,2)
  for i= 1:size(res.pred_angle,1)
    if res.pred_angle(i,j) <= 40 && res.pred_angle(i,j) >= -40
      res.pred_pos(i,j)= 1;
      
    elseif res.pred_angle(i,j) >= 50 && res.pred_angle(i,j) <= 130
      res.pred_pos(i,j)= 2;
      
    elseif res.pred_angle(i,j) >= 140 || res.pred_angle(i,j) <= -140
      res.pred_pos(i,j)= 3;
      
    elseif res.pred_angle(i,j) <= -50 && res.pred_angle(i,j) >= -130
      res.pred_pos(i,j)= 4;
      
    else
      res.pred_pos(i,j)= NaN;
      
    end
  end
  
  res.bcentre = reshape([cdata.ballcentre], 2, [])';
  
  
  % for i = 1:numel(res.crab)
  %
  %     if res.group(i) == 1
  %         res.ind(i) = res.crab(i);
  %     elseif res.group(i) == 2
  %         res.ind(i) = res.crab(i) + 10;
  %     elseif res.group(i) == 3
  %         res.ind(i) = res.crab(i) + 20;
  %     elseif res.group(i) == 4
  %         res.ind(i) = res.crab(i) + 30;
  %     end
  %
  % end
  
  % res.ind = res.ind';
  
  %create day variable for random effects model
  res.day = zeros(size(res.trial));
  res.day(ismember(res.trial, [1 2]) & ismember(res.group, 1)) = 1;
  res.day(ismember(res.trial, [3 4]) & ismember(res.group, 1)) = 3;
  res.day(ismember(res.trial, [5 6]) & ismember(res.group, 1)) = 5;
  res.day(ismember(res.trial, [7 8]) & ismember(res.group, 1)) = 7;
  res.day(ismember(res.trial, 9) & ismember(res.group, 1)) = 9;
  
  res.day(ismember(res.trial, [1 2]) & ismember(res.group, 2)) = 2;
  res.day(ismember(res.trial, [3 4]) & ismember(res.group, 2)) = 4;
  res.day(ismember(res.trial, [5 6]) & ismember(res.group, 2)) = 6;
  res.day(ismember(res.trial, [7 8]) & ismember(res.group, 2)) = 8;
  res.day(ismember(res.trial, 9) & ismember(res.group, 2)) = 10;
  
  res.day(ismember(res.trial, [1 2]) & ismember(res.group, 3)) = 1;
  res.day(ismember(res.trial, [3 4]) & ismember(res.group, 3)) = 2;
  res.day(ismember(res.trial, [5 6]) & ismember(res.group, 3)) = 3;
  res.day(ismember(res.trial, [7 8]) & ismember(res.group, 3)) = 4;
  res.day(ismember(res.trial, [9 10]) & ismember(res.group, 2)) = 5;
  
  
  for i = 1:numel(res.mm)
    % parameters that never change can be here for all of them
    % for matching stimulus calculations
    res.stim_para(i,1).distance=5000; %initial distance of the reference stimulus in mm
    res.stim_para(i,1).speed=200;%the speed of the reference stimulus in mm/s
    res.stim_para(i,1).equal_distance=445.04;%the distance that both stimuli would be there at the exact0 same time
    res.stim_para(i,1).frameRate=60;%monitors frame rate
    res.stim_para(i,1).elevation=15; %elevation of stimulus at initial position in deg
    res.stim_para(i,1).finalDistance=15; %final distance of stimulus from the observer
    res.stim_para(i,1).traveling_distance=5000-res.stim_para(i).finalDistance; % keep at end
    res.stim_para(i,1).travel_time = res.stim_para(i).traveling_distance/res.stim_para(i).speed;
    res.stim_para(i,1).XYangle=45;%Angle of end position in XY plane relative to X axis, you only need this if you have a near miss stimulus
    res.stim_para(i,2) = res.stim_para(i,1);
    % for normal and reference stimuli
    res.time(i,:) = 24.925;
    res.speed(i,:) = 200;  % keep at reference speed even for matching stimuli
    res.size(i,:) = [30 nan];
    res.elevation(i,:) = 15;
    
    switch(res.mm(i))
      case 1 % 22mm constant speed
        res.ismatching_stim(i,:) = [0 nan];
        res.size(i,:) = [22 nan];
        
      case 2 % 22mm mimicking 30mm constant speed
        res.ismatching_stim(i,:) = [1 nan];
        res.size(i,:) = [22 nan];
        res.stim_para(i).dr = 30;%diameter of the reference stimulus in mm
        res.stim_para(i).d = res.size(i,1);%diameter of the matching stimulus in mm
        zdistance = 704630;
        
      case 3 % 30mm constant speed
        res.ismatching_stim(i,:) = [0 nan];
        
      case 4 % 30mm mimicking 22mm constant speed
        res.ismatching_stim(i,:) = [1 nan];
        res.stim_para(i,1).dr = 22;%diameter of the reference stimulus in mm
        res.stim_para(i,1).d = res.size(i,1);%diameter of the matching stimulus in mm
        zdistance = 1341.6;
        
      case 5 % 45mm mimicking 22mm constant speed
        res.ismatching_stim(i,:) = [1 nan];
        res.size(i,:) = [45 nan];
        res.stim_para(i,1).dr = 22;%diameter of the reference stimulus in mm
        res.stim_para(i,1).d = res.size(i,1);%diameter of the matching stimulus in mm
        zdistance = 803.3;
        
       case 6 % Dual, 22mm constant speed + 30mm constant speed
        res.ismatching_stim(i,:) = [0, 0];
        res.size(i,:) = [22 30];
        
       case 7 % Dual, 22mm constant speed + 30mm mimicking 22mm constant speed
        res.ismatching_stim(i,:) = [0, 1];
        res.size(i,:) = [22, 30];
        res.stim_para(i,2).dr = 22;%diameter of the reference stimulus in mm
        res.stim_para(i,2).d = res.size(i,2);%diameter of the matching stimulus in mm
        zdistance = [5000, 1341.6];
        
       case 8 % Dual, 22mm constant speed + 45mm mimicking 22mm constant speed
        res.ismatching_stim(i,:) = [0, 1];
        res.size(i,:) = [22, 45];
        res.stim_para(i,2).dr = 22;%diameter of the reference stimulus in mm
        res.stim_para(i,2).d = res.size(i,2);%diameter of the matching stimulus in mm
        zdistance = [5000, 803.3];
        
       case 9 % Dual, 30mm constant speed + 22mm mimicking 30mm constant speed
        res.ismatching_stim(i,:) = [0, 1];
        res.size(i,:) = [30, 22];
        res.stim_para(i,2).dr = 30;%diameter of the reference stimulus in mm
        res.stim_para(i,2).d = res.size(i,2);%diameter of the matching stimulus in mm
        zdistance = [5000, 704630];
        
       case 10 % Dual, 30mm constant speed + 45mm mimicking 30mm constant speed
        res.ismatching_stim(i,:) = [0, 1];
        res.size(i,:) = [30, 45];
        res.stim_para(i,2).dr = 30;%diameter of the reference stimulus in mm
        res.stim_para(i,2).d = res.size(i,2);%diameter of the matching stimulus in mm
        zdistance = [5000, 1134.9];
        
       case 11 % Dual, 30mm mimicking 22mm constant speed + 45mm mimicking 22mm constant speed
        res.ismatching_stim(i,:) = [1, 1];
        res.size(i,:) = [30, 45];
        res.stim_para(i,1).dr = 22;%diameter of the reference stimulus in mm
        res.stim_para(i,1).d = res.size(i,1);%diameter of the matching stimulus in mm
        res.stim_para(i,2).dr = 22;%diameter of the reference stimulus in mm
        res.stim_para(i,2).d = res.size(i,2);%diameter of the matching stimulus in mm
        zdistance = [1341.6, 803.3];
        
    end
    % do common calculations
    for j = 1:numel(zdistance)
      res.startpos(i,(1:3)+(j-1)*3) = [960 zdistance(j)*sind(res.elevation(i,:)) zdistance(j)];
    end
  end
  res.duration = (res.startpos(:,3)./res.speed);
  res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;

  if framerate == -1
    est_FR = (res.runse - res.runss)./res.duration;
    res.framerate = est_FR;
    
  end
  
  
  % Turning pm data into a data matrix that is meaningful for analysis
  % Extract behaviours, return to a matrix and warn if more than one type of behaviour occurs in one trial
  
  res.pmr=[];
  res.pml=[];
  res.pmo=[];
  res.pmw=[];
  res.pme=[];
  res.pms=[];
  res.pmc=[];
  res.pmf=[];
  nr = inf;
  
  for i = 1:numel(cdata)
    %     if cdata(i).crab == 2 && cdata(i).trial==3 && cdata(i).group==1
    %         keyboard
    %     end
    
    if ~isempty(cdata(i).pm)
      cdata(i).pm_all(cdata(i).pm(:,1)<=-200,:)=[]; % remove all points before -200
      cdata(i).pm(cdata(i).pm(:,1)<=-200,:)=[]; % remove all points before -200
      
      sel = find(cdata(i).pm(:,2) == 'r');
      if isempty(sel)
        res.pmr(i,1) = nr;
      elseif numel(sel)>1
        %         fprintf('warning, crab %d trial %d has more than one run\n', res.crab(i), res.trial(i));
        res.pmr(i,1) = cdata(i).pm(sel(1), 1);
      else
        res.pmr(i,1) = cdata(i).pm(sel(1), 1);
      end
      
      sel = find(cdata(i).pm(:,2) == 'l');
      if isempty(sel)
        res.pml(i,1) = nr;
      elseif numel(sel)>1
        %         fprintf('warning, crab %d trial %d has more than one l\n', res.crab(i), res.trial(i));
        res.pml(i,1) = cdata(i).pm(sel(1), 1);
      else
        res.pml(i,1) = cdata(i).pm(sel, 1);
      end
      
      sel = find(cdata(i).pm(:,2) == 'o');
      if isempty(sel)
        res.pmo(i,1) = nr;
      elseif numel(sel)>1
        %         fprintf('warning, crab %d trial %d has more than one o\n', res.crab(i), res.trial(i));
        res.pmo(i,1) = cdata(i).pm(sel(1), 1);
        %        res.pmo(i,1) = cdata(i).pm(sel(1), 1);
      else
        res.pmo(i,1) = cdata(i).pm(sel, 1);
      end
      
      
      sel = find(cdata(i).pm(:,2) == 'w');
      if isempty(sel)
        res.pmw(i,1) = nr;
      elseif numel(sel)>1
        %         fprintf('warning, crab %d trial %d has more than one w\n', res.crab(i), res.trial(i));
        res.pmw(i,1) = cdata(i).pm(sel(1), 1);
      else
        res.pmw(i,1) = cdata(i).pm(sel, 1);
      end
      
      
      sel = find(cdata(i).pm(:,2) == 'e');
      if isempty(sel)
        res.pme(i,1) = nr;
      elseif numel(sel)>1
        %         fprintf('warning, crab %d trial %d has more than one e\n', res.crab(i), res.trial(i));
        res.pme(i,1) = cdata(i).pm(sel(1), 1);
      else
        res.pme(i,1) = cdata(i).pm(sel, 1);
      end
      
      sel = find(cdata(i).pm(:,2) == 's');
      if isempty(sel)
        res.pms(i,1) = nr;
      elseif numel(sel)>1
        %         fprintf('warning, crab %d trial %d has more than one s\n', res.crab(i), res.trial(i));
        res.pms(i,1) = cdata(i).pm(sel(1), 1);
      else
        res.pms(i,1) = cdata(i).pm(sel, 1);
      end
      
      
      sel = find(cdata(i).pm(:,2) == 'c');
      if isempty(sel)
        res.pmc(i,1) = nr;
      elseif numel(sel)>1
        %         fprintf('warning, crab %d trial %d has more than one c\n', res.crab(i), res.trial(i));
        res.pmc(i,1) = cdata(i).pm(sel(1), 1);
      else
        res.pmc(i,1) = cdata(i).pm(sel, 1);
      end
      
      sel = find(cdata(i).pm(:,2) == 'f');
      if isempty(sel)
        res.pmf(i,1) = nr;
      elseif numel(sel)>1
        %         fprintf('warning, crab %d trial %d has more than one f\n', res.crab(i), res.trial(i));
        res.pmf(i,1) = cdata(i).pm(sel(1), 1);
      else
        res.pmf(i,1) = cdata(i).pm(sel, 1);
      end
    end
  end
  %find crab/trial number for early double behaviour responses
  
  
  
  %%%%%%%%% Now we can get time of earliest response for each trial
  res.behav_order = [{'Freeze'}, {'Run'}, {'Limb'}];
  allresp = [res.pmf res.pmr res.pml];
  
  if numel(res.behav_order) ~= size(allresp,2)
    
    warning('Number of elements of behaviour order string does not match the size of the behaviour matrix')
    
  end
  
  %%%%%%%%%% exclude early (assume a crab response time of 160ms)
  sel = allresp <= (delay*res.framerate) & isfinite(allresp);
  fprintf('removed %d early response(s)\n', sum(sel(:)));
  allresp(sel)=inf; %% now for responsetime, only consider times when the stimulus is on.
  % check=[]; %% creates a dataframe with the crab, trial, and protocol with
  % the early responses. Only used for double checking yourself.
  % for uuu=1:size(sel,1)
  %  if sel(uuu,:)~=zeros(1,size(sel,2))
  %    check=[check; cdata(uuu).crab cdata(uuu).trial cdata(uuu).mm];
  %  end
  %end
  
  % %%%%%%%%%%%% exclude late ones (with a crab response time of delay)
  sel = (res.rm(:,size(res.rm,2)*ones(1,size(allresp,2))) <= (allresp-(delay*res.framerate))) & isfinite(allresp);
  fprintf('removed %d late response(s)\n', sum(sel(:)));
  allresp(sel)=inf; %%%%%%% now for responsetime, only consider times when the stimulus is on.
  
  if make_defence == 1
    % Locates where in the behaviour order string, there is not run and
    % freeze
    locate = contains(res.behav_order, 'run', 'IgnoreCase', true) | contains(res.behav_order, 'freeze', 'IgnoreCase', true);
    locate_def = locate ~= 1;
    
    temp_resp = allresp(:,locate_def);
    
    for i = 1:size(allresp,1)
      
      sel = ~isinf(temp_resp(i,:));
      
      if sum(sel)> 0
        
        temp_defence(i) = min(temp_resp(i,sel), [],2);
      else
        temp_defence(i) = Inf;
      end
    end
    allresp = [allresp(:,locate), temp_defence'];
    res.behav_order = [{'Freeze'}, {'Run'}, {'Defence'}];
    
    if numel(res.behav_order) ~= size(allresp,2)
      warning('Number of elements of behaviour order string does not match the size of the behaviour matrix')
    end
  end
  
  
  
  % Adjust all bhaviours for delay
  % loop because of mimicking stimuli
  res.behav_frame = allresp - (delay*res.framerate);
  res.t2c=nan(size(res.behav_frame,1),size(res.behav_frame,2));
  res.dist=res.t2c;
  res.asize=res.t2c;
  res.edgespeed=res.t2c;
  res.espeed=res.t2c;
  res.area_expansion=res.t2c;
  res.anginc=res.t2c;
  res.anginc2=res.t2c;
  for stim = 0:numel(res.size(i,:))-1
    for run = 1:size(res.behav_frame, 1)
      for i = 1:size(res.behav_frame, 2)
        response_time = res.behav_frame(run,i)./res.framerate;
        res.response_time(run,i) = response_time(run);
        if res.ismatching_stim(run, stim+1) == 1
          if isfinite(res.behav_frame(run, i))
            out = find_match_stim_para(res.stim_para(run, stim+1), res.rm(run, 2), res.rm(run,3), res.behav_frame(run,i));
            res.asize(run,i) = out.asize_match;
            res.espeed(run,i) = out.aspeed_match;
          else
            res.asize(run,i) = inf;
            res.espeed(run,i) = inf;
          end
          res.t2c(run,i) = -((res.rm(run,3)-res.behav_frame(run,i))./res.framerate(run)); % ref stim only
          res.dist(run,i) = res.startpos(run,3+stim*3)-(res.speed(run).*res.response_time(run,i));
          res.edgespeed(run,i) = nan;
          res.area_expansion(run,i) = nan;
          res.anginc(run,i) = nan;
          res.anginc2(run,i) = nan;
        else
          res.t2c(run,i) = -((res.rm(run,3)-res.behav_frame(run,i))./res.framerate(run));
          res.dist(run,i) = res.startpos(run,3)-(res.speed(run).*res.response_time(run,i));
          res.asize(run,i) = (2.*(atan(res.size(run, stim+1)/2./res.dist(run,i))))/pi*180;
          res.area(run,i) = (2.*(atan((pi.*(res.size(run, stim+1)/2).^2)./res.dist(run,i))))/pi*180;
          dtemp = (res.startpos(run,3+stim*3)-(res.speed(run).*((res.behav_frame(run,i)-1)./res.framerate(run)))); % one timestep before response
          temp = (2.*(atan(res.size(run, stim+1)/2./dtemp)))/pi*180;
          res.edgespeed(run,i) = abs(((res.asize(run,i) - temp)./2).*res.framerate(run)); % deg/sec for one edge %gives absolute values
          res.espeed(run,i) = abs((res.asize(run,i) - temp).*res.framerate(run)); % deg/sec for one edge %gives absolute values
          dtemp2 = (res.startpos(run,3+stim*3)-(res.speed(run).*((res.behav_frame(run,i)-1)./res.framerate(run)))); % one timestep before response
          areatemp = (2.*(atan((pi.*(res.size(run, stim+1)/2).^2)./dtemp2)))/pi*180;
          res.area_expansion(run,i) = abs((res.area(run,i) - areatemp).*res.framerate(run)./2); % deg/sec for one edge %gives absolute values
          res.anginc(run,i) = res.asize(run,i) - (2.*atand(res.size(run, stim+1)/2./res.startpos(run,3+stim*3))); %angsize at response minus initial angsize
          res.anginc2(run,i) = res.asize(run,i) - res.start_asize(run); %angsize at response minus initial angsize
        end
      end
    end
  end
  sel = isinf(res.response_time);
  res.asize(sel) = Inf;
  res.espeed(sel) = Inf;
  res.area(sel) = Inf;
  res.area_expansion(sel) = Inf;
  res.anginc(sel) = Inf;
  
end
