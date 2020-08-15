function block__plot(blocklist, mcluster, children)
%
%   drawing the blocks given in blocklist
%   
    N = mcluster(end,2) - mcluster(end,1) + 1;
    
    %   Frame
    fig = figure;
    rectangle('Position', [0,0,N,N]);
    axis([0 N 0 N]);
    hold on
    %   Off-diagonal blocks
    for i = 1 : size(blocklist, 1);
        c1 = blocklist(i, 1);
        c2 = blocklist(i, 2);
        x0 = mcluster(c2, 1) - 1;
        y0 = (N+1) - mcluster(c1, 2) - 1;        
        h = (N+1-mcluster(c1, 1)) - (N+1-mcluster(c1, 2)) + 1;
        w = mcluster(c2, 2) - mcluster(c2, 1) + 1;
        rectangle('Position', [x0,y0,w,h], 'FaceColor',[220,220,220]/255)
    end
    
    %   Diagonal blocks
    if (nargin == 3)
        for i = 1 : size(mcluster, 1)
            if ~isnan(children(i,1)); continue; end;
            x0 = mcluster(i, 1) - 1;
            y0 = (N+1) - mcluster(i, 2) - 1;        
            h = (N+1-mcluster(i, 1)) - (N+1-mcluster(i, 2)) + 1;
            w = mcluster(i, 2) - mcluster(i, 1) + 1;
            rectangle('Position', [x0,y0,w,h], 'FaceColor',[0 0 0])
        end
    end
end


