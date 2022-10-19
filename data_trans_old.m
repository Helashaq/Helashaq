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
    switch(res.mm(i))
        case 1 % full black square, start size different (10 mm more than calc. prob.)
            res.size(i,:) = 30;
            res.speed(i,:) = 200;
%            res.threedspeed(i,:) = [0 -51.7638 -193.1852];
            res.time(i,:) = 24.925;
            res.elevation(i,:) = 15;
            zdistance = 5000;
            res.startpos(i,:) = [960 zdistance*sind(res.elevation(i,:)) zdistance];
            res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;
            
            
        case 2 % empty square, edges only, start size different 
            res.size(i,:) = 30;
            res.speed(i,:) = 200;
%            res.threedspeed(i,:) = [0 -51.7638 -193.1852];
            res.time(i,:) = 24.925;
            res.elevation(i,:) = 15;
            zdistance = 5000;
            res.startpos(i,:) = [960 zdistance*sind(res.elevation(i,:)) zdistance];
            res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;
            
        case 3 % empty square, 2x edges only, start size different
            res.size(i,:) = 30;
            res.speed(i,:) = 200;
%            res.threedspeed(i,:) = [0 -51.7638 -193.1852];
            res.time(i,:) = 24.925;
            res.elevation(i,:) = 15;
            zdistance = 5000;
            res.startpos(i,:) = [960 zdistance*sind(res.elevation(i,:)) zdistance];
            res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;
            
        case 4 % square, expanding horizontal only 
            res.size(i,:) = 30;
            res.speed(i,:) = 200;
%            res.threedspeed(i,:) = [0 -51.7638 -193.1852];
            res.time(i,:) = 24.925;
            res.elevation(i,:) = 15;
            zdistance = 30000;
            res.startpos(i,:) = [960 zdistance*sind(res.elevation(i,:)) zdistance];
            res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;
            
        case 5 % square, expanding vertical only 
            res.size(i,:) = 30;
            res.speed(i,:) = 200;
%            res.threedspeed(i,:) = [0 -51.7638 -193.1852];
            res.time(i,:) = 24.925;
            res.elevation(i,:) = 15;
            zdistance = 30000;
            res.startpos(i,:) = [960 zdistance*sind(res.elevation(i,:)) zdistance];
            res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;
            
        case 6 % standard circle 
            res.size(i,:) = 30;
            res.speed(i,:) = 200;
%            res.threedspeed(i,:) = [0 -51.7638 -193.1852];
            res.time(i,:) = 24.925;
            res.elevation(i,:) = 15;
            zdistance = 5000;
            res.startpos(i,:) = [960 zdistance*sind(res.elevation(i,:)) zdistance];
            res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;

        case 7 % circle, end point 10 cm
            res.size(i,:) = 30;
            res.speed(i,:) = 200;
%            res.threedspeed(i,:) = [0 -51.7638 -193.1852];
            res.time(i,:) = 24.925;
            res.elevation(i,:) = 15;
            zdistance = 5000;
            res.startpos(i,:) = [960 zdistance*sind(res.elevation(i,:)) zdistance];
            res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;
            

        case 8 % circle, end point 10 cm, mimicking speed 
            res.size(i,:) = 30;
            res.speed(i,:) = 200;
%            res.threedspeed(i,:) = [0 -51.7638 -193.1852];
            res.time(i,:) = 24.925;
            res.elevation(i,:) = 15;
            zdistance = 5000;
            res.startpos(i,:) = [960 zdistance*sind(res.elevation(i,:)) zdistance];
            res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;
            

       case 9  % circle, end point 10 cm
            res.size(i,:) = 30;
            res.speed(i,:) = 200;
%            res.threedspeed(i,:) = [0 -51.7638 -193.1852];
            res.time(i,:) = 24.925;
            res.elevation(i,:) = 15;
            zdistance = 5000;
            res.startpos(i,:) = [960 zdistance*sind(res.elevation(i,:)) zdistance];
            res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;
            

        case 10 % standard circle + end point 10 cm
            res.size(i,:) = 30;
            res.speed(i,:) = 200;
%            res.threedspeed(i,:) = [0 -51.7638 -193.1852];
            res.time(i,:) = 24.925;
            res.elevation(i,:) = 15;
            zdistance = 5000;
            res.startpos(i,:) = [960 zdistance*sind(res.elevation(i,:)) zdistance];
            res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;
            

        case 11 % standard circle + end point 10 cm, mimicking speed 
            res.size(i,:) = 30;
            res.speed(i,:) = 200;
%            res.threedspeed(i,:) = [0 -51.7638 -193.1852];
            res.time(i,:) = 24.925;
            res.elevation(i,:) = 15;
            zdistance = 5000;
            res.startpos(i,:) = [960 zdistance*sind(res.elevation(i,:)) zdistance];
            res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;
            

%          case 12 % direct circle mimicking speed + end point 10 cm
%             res.size(i,:) = 30;
%             res.speed(i,:) = 200;
% %            res.threedspeed(i,:) = [0 -51.7638 -193.1852];
%             res.time(i,:) = 24.925;
%             res.elevation(i,:) = 15;
%             zdistance = 5000;
%             res.startpos(i,:) = [960 zdistance*sind(res.elevation(i,:)) zdistance];
%             res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;
%             
% 
%          case 13 % standard circle + end point 1.5 cm
%             res.size(i,:) = 30;
%             res.speed(i,:) = 200;
% %            res.threedspeed(i,:) = [0 -51.7638 -193.1852];
%             res.time(i,:) = 24.925;
%             res.elevation(i,:) = 15;
%             zdistance = 5000;
%             res.startpos(i,:) = [960 zdistance*sind(res.elevation(i,:)) zdistance];
%             res.start_asize(i,:) = atand(res.size(i,:)/zdistance)*2;
            


            
    end
end

res.duration = (res.startpos(:,3)./res.speed);

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
res.behav_frame = allresp - (delay*res.framerate);


for i = 1:size(res.behav_frame,2)
    
    response_time = res.behav_frame(:,i)./res.framerate;
    res.response_time(:,i) = response_time;
    
    res.t2c(:,i) = -((res.rm(:,3)-res.behav_frame(:,i))./res.framerate);
    res.dist(:,i) = res.startpos(:,3)-(res.speed.*res.response_time(:,i));
    res.asize(:,i) = (2.*(atan(res.size/2./res.dist(:,i))))/pi*180;
    res.area(:,i) = (2.*(atan((pi.*(res.size/2).^2)./res.dist(:,i))))/pi*180;
    
    
    dtemp = (res.startpos(:,3)-(res.speed.*((res.behav_frame(:,i)-1)./res.framerate))); % one timestep before response
    temp = (2.*(atan(res.size/2./dtemp)))/pi*180;
    
    res.edgespeed(:,i) = abs(((res.asize(:,i) - temp)./2).*res.framerate); % deg/sec for one edge %gives absolute values
    
    res.espeed(:,i) = abs((res.asize(:,i) - temp).*res.framerate); % deg/sec for one edge %gives absolute values
    
    % res.espeed(:,i) = rad2deg(((res.size./2)./res.speed)./(res.t2c(:,i).^2 + ((res.size./2)./res.speed).^2));
    
    dtemp2 = (res.startpos(:,3)-(res.speed.*((res.behav_frame(:,i)-1)./res.framerate))); % one timestep before response
    areatemp = (2.*(atan((pi.*(res.size/2).^2)./dtemp2)))/pi*180;
    res.area_expansion(:,i) = abs((res.area(:,i) - areatemp).*res.framerate./2); % deg/sec for one edge %gives absolute values
    
    res.anginc(:,i) = res.asize(:,i) - (2.*atand(res.size/2./res.startpos(:,3))); %angsize at response minus initial angsize
    res.anginc2(:,i) = res.asize(:,i) - res.start_asize; %angsize at response minus initial angsize

    
end

sel = isinf(res.response_time);
res.asize(sel) = Inf;
res.espeed(sel) = Inf;
res.area(sel) = Inf;
res.area_expansion(sel) = Inf;
res.anginc(sel) = Inf;

end
