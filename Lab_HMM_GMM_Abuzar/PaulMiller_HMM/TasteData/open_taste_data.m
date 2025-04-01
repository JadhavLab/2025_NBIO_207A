% open_taste_data.m

%ps_array's is a spike array with dimensions (4,30,30,7000) corresponding 
% to (Taste,Trial,Neuron,Millisecond) with taste deliveries occurring at 
% 2000ms, and the order of taste dimension is NaCl(1), Sucrose(2), 
% Citric Acid(3) and Quinine(4). 
%
% These data corresponds to one session of tastes delivered in random order 
% 20s apart. (Data from the intertrial interval, or the "true" linear order 
% of the taste deliveries is destroyed in the dataset). 
%
% Also, you will have to be careful, since a handful of trials have no data 
% (just composed of zeros) due to noise that occurred during the trial, 
% so you will have to filter those trials out, or deal with them somehow. 


peristim_mat = load('BS26_4Tastes_180204_102132_Peristim_all_tastes.mat');
ps_init = fieldnames(peristim_mat.all_tastes); 
ps_array = peristim_mat.all_tastes.(ps_init{1}); 
clear ps_init

