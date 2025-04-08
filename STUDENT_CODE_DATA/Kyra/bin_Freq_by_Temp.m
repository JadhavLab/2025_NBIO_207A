
load kyra_proc_data1.mat

%Temp_bins = [10:1:32];
Temp_bins = [10:2:32];

Freq1_hist_bytemp=zeros(1,length(Freq_bins)-1);
Freq1_bytemp = [];
% Temp1, Freq1
for t=1:length(Temp_bins)-1
    curr_trange = Temp_bins(t):Temp_bins(t)+1;
    x_bins = find(Temp1>curr_trange(1) & Temp1<curr_trange(2));
    if ~isempty(x_bins)
        curr_freq = Freq1(x_bins);
        Freq1_bytemp(t) = mean(curr_freq);
        %[n_freq_curr] = histcounts(curr_freq, Freq_bins);
        %Freq1_hist_bytemp = Freq1_hist_bytemp + n_freq_curr;
    else
        Freq1_bytemp(t) = 0;
    end

end

figure(100); hold on; 
plot(Temp_bins(1:end-1),Freq1_bytemp,'ro-', 'MarkerSize', 8)


% Temp2, Freq2
Freq2_bytemp = [];
for t=1:length(Temp_bins)-1
    curr_trange = Temp_bins(t):Temp_bins(t)+1;
    x_bins = find(Temp2>curr_trange(1) & Temp2<curr_trange(2));
    if ~isempty(x_bins)
        curr_freq = Freq2(x_bins);
        Freq2_bytemp(t) = mean(curr_freq);
    else
        Freq2_bytemp(t) = 0;
    end
end

figure(100); hold on; 
plot(Temp_bins(1:end-1),Freq2_bytemp,'bo-', 'MarkerSize', 8)


% Temp3, Freq3
Freq3_bytemp = [];
for t=1:length(Temp_bins)-1
    curr_trange = Temp_bins(t):Temp_bins(t)+1;
    x_bins = find(Temp3>curr_trange(1) & Temp3<curr_trange(2));
    if ~isempty(x_bins)
        curr_freq = Freq3(x_bins);
        Freq3_bytemp(t) = mean(curr_freq);
    else
        Freq3_bytemp(t) = 0;
    end
end

figure(100); hold on; 
plot(Temp_bins(1:end-1),Freq3_bytemp,'gsq-', 'MarkerSize', 12)


% Temp4, Freq4
Freq4_bytemp = [];
for t=1:length(Temp_bins)-1
    curr_trange = Temp_bins(t):Temp_bins(t)+1;
    x_bins = find(Temp4>curr_trange(1) & Temp4<curr_trange(2));
    if ~isempty(x_bins)
        curr_freq = Freq4(x_bins);
        Freq4_bytemp(t) = mean(curr_freq);
    else
        Freq4_bytemp(t) = 0;
    end
end

figure(100); hold on; 
plot(Temp_bins(1:end-1),Freq4_bytemp,'csq-', 'MarkerSize', 12)



% Temp5, Freq5
Freq5_bytemp = [];
for t=1:length(Temp_bins)-1
    curr_trange = Temp_bins(t):Temp_bins(t)+1;
    x_bins = find(Temp5>curr_trange(1) & Temp5<curr_trange(2));
    if ~isempty(x_bins)
        curr_freq = Freq5(x_bins);
        Freq5_bytemp(t) = mean(curr_freq);
    else
        Freq5_bytemp(t) = 0;
    end
end
figure(100); hold on; 
plot(Temp_bins(1:end-1),Freq5_bytemp,'mx-', 'MarkerSize', 12)


% Temp6, Freq6
Freq6_bytemp = [];
for t=1:length(Temp_bins)-1
    curr_trange = Temp_bins(t):Temp_bins(t)+1;
    x_bins = find(Temp6>curr_trange(1) & Temp6<curr_trange(2));
    if ~isempty(x_bins)
        curr_freq = Freq6(x_bins);
        Freq6_bytemp(t) = mean(curr_freq);
    else
        Freq6_bytemp(t) = 0;
    end
end
figure(100); hold on; 
%plot(Temp_bins(1:end-1),Freq6_bytemp,'kx', 'MarkerSize', 12)




