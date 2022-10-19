function allout = dataset_loop(dataset_dir, fictrac_dir, delay, pred, crab, ball, ori, imcorrect, claw_pos)


dataset_dir = {dataset_dir};
allout= struct([]);

for ds = 1:numel(dataset_dir)
    nn = dir(fullfile(dataset_dir{ds}, 'resfiles'));
    nn(1:2)=[];

    for wname = 1:numel(nn)
        if isfolder(fullfile(dataset_dir{ds}, 'resfiles', nn(wname).name))
            %          digi(dataset, nn(wname).name, 'loadfull')
            %          s = savedata(s);
            out = read_digi_data(dataset_dir{ds}, nn(wname).name, fictrac_dir, delay, pred, crab, ball, ori, imcorrect, claw_pos);
            allout = [allout out];
        end
    end
    
    
    
    fprintf('%10s\n', allout.sname)
    fprintf('collected %d runs from directory\n', numel(allout))
end

 