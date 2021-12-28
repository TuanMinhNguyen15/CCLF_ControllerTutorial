function dVdz = dVdzCal(x,y,F,d,n)
    %d = 0.001;
    %n = 8;

    %x = -1;
    xmax = x+d;
    xmin = x-d;
    x_grid = xmin:2*d/n:xmax;

    %y = 1;
    ymax = y+d;
    ymin = y-d;
    y_grid = ymin:2*d/n:ymax;

    Pq = [];

    for i=1:length(x_grid)
        for j = 1:length(y_grid)
            Pq = [Pq;x_grid(i) , y_grid(j)];
        end
    end


    Vq = F(Pq);
    mdl = fitlm(Pq,Vq);
    result_table = mdl.Coefficients;
    result_array = table2array(result_table);

    dVdz = [result_array(2,1) , result_array(3,1)];
end