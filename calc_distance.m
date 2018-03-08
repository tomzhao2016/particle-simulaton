%% 
fName = strcat('gpu.txt');
fid = fopen(fName);
C = textscan(fid, '%s', 'delimiter', '\n');

%%
num_time = 100;  % saved every 10 time steps
num_particle = 2000;
myx = zeros(num_particle, num_time);
myy = zeros(num_particle, num_time);
for k = 1:size(C{1,1}, 1) - 1
    tmp = C{1, 1}{k + 1, 1};
    tmp = strsplit(tmp);
    current_time = ceil(k/num_particle);
    if (mod(k, num_particle) == 0)
        current_particle = num_particle;
    else
        current_particle = mod(k, num_particle);
    end
    myx(current_particle, current_time) = str2double(tmp{1});
    myy(current_particle, current_time) = str2double(tmp{2});
end
fclose(fid);

%%
% calculate the average/min distance in each time step
cutoff = 0.01;
distance_avg = zeros(num_time, 1);
distance_min = zeros(num_time, 1);
for i = 1:num_time
    waitbar(i/num_time);
    temp_avg = 0;
    temp_min = Inf;
    count = 0;
    for j = 1:num_particle
        for k = j+1:num_particle
            temp_dis = sqrt((myx(j, i) - myx(k, i))^2 + (myy(j, i) - myy(k, i))^2);
            if(temp_dis <= cutoff)
                temp_avg = temp_avg + temp_dis;
                count = count + 1;
            end
            if(temp_dis < temp_min)
                temp_min = temp_dis;
            end
        end
    end
%     if(count ~= (num_particle^2 - num_particle)/2)
%         warning('possible errors');
%     end
    if(count > 0)
        temp_avg = temp_avg/count;
    end
    distance_avg(i) = temp_avg;
    distance_min(i) = temp_min;
end



%% divide by the cutoff
distance_avg = distance_avg/cutoff;
distance_min = distance_min/cutoff;
distance_avg_time = mean(distance_avg)
distance_min_time = min(distance_min)


%%
scatter(myx(:, 1), myy(:, 1))

%%
scatter(myx(:, 100), myy(:, 100))


