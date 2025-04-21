% fourier_transform.m
%
% Code to do the discrete Fourier transform in a simple-minded (slow) way 
% that allows the code to be altered easily.
%
% data_array should contain one or more columns of data, where each column
% is a time series
load('data_array.mat');

[Nt Nvar] = size(data_array);   % Nt = no. of rows, Nvar = No. of columns
tvec = 1:Nt;                    % Time vector in units of time-step (mins)

f_vec = 1/Nt:1/Nt:1;              % frequencies to calculate Fourier Transform
Nf = length(f_vec);             % no. of frequencies
sin_array = zeros(Nt,Nf);          % sine array initialized as NtxNf
cos_array = zeros(Nt,Nf);          % cosine array initialized as NtxNf

i = 0;                          % to count column in sine/cosine array
for f = f_vec                   % update frequency to calculate
    i = i+1;                    % column number
    sin_array(:,i) = sin(2*pi*f*tvec);  % vector of sine at given freq
    cos_array(:,i) = cos(2*pi*f*tvec);  % vector of cosine at given freq
end

first_col = 1;                      % which columns of array to look at
last_col = 1;
Ncols = last_col-first_col+1;       % Number of columns to look at

sf = zeros(Nf,Ncols);               % will contain the sine coefficients
cf = zeros(Nf,Ncols);               % will contain the cosine coefficients
power = zeros(Nf,Ncols);            % will contain the power

j = 0;                              % column to store results
for col_num = first_col:last_col;   % loop through data_array columns
    j = j + 1;                      % store results in next column
    for i = 1:Nf;                   
        sf(i,j) = mean(sin_array(:,i).*data_array(:,col_num));  % overlap with sine
        cf(i,j) = mean(cos_array(:,i).*data_array(:,col_num));  % overlap with cosine
    end
end

power = sf.*sf + cf.*cf;        % power is sum of squared coefficients

j = 0;
for col_num = first_col:last_col;   % for each column of data used 
    j = j+1
    figure(j)                       % plot a new figure
    plot(f_vec,power(:,j))          % of power versus frequency
end
