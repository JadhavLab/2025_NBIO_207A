% GLM fitting

%load(ripplemodfile);
load('HParipplemod01.mat');
%load(cellinfofile);
load('HPAcellinfo.mat')
day=1; epoch=4;

filterString = '(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && ($meanrate < 7)';
cellindices = evaluatefilter(cellinfo{day}{epoch}, filterString);
cellsi = [repmat([day epoch], size(cellindices,1),1 ), cellindices]; % day-epoch-tet-cell for CA1 cells
usecellsi = 1:size(cellsi,1); nCA1cells = size(cellsi,1);

PFCtet=[];
% PFC cells
if ~isempty(PFCtet),
    cellsp = [day epoch PFCtet PFCcell];
    usecellsp = 1;
else
    %filterString = 'strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')';
    filterString = 'strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($ripmodtag2, ''y'') ';
    pcellindices = evaluatefilter(cellinfo{day}{epoch}, filterString);
    cellsp = [repmat([day epoch], size(pcellindices,1),1 ), pcellindices]; % day-epoch-tet-cell for CA1 cells
    usecellsp = 1:size(cellsp,1); nPFCcells = size(cellsp,1);   
end



% Get ripplemoddata
for i=1:size(cellsi,1)
    i;
    eval(['ripplemodi{',num2str(i),'}= ripplemod{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
        '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.trialResps;']);
end

for i=1:size(cellsp,1)
    i;
    eval(['ripplemodp{',num2str(i),'}= ripplemod{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
        '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.trialResps;']);
%     eval(['spiketimep{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
%         '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,1);']);
end




% Get glmfits for the model using Poisson distr, which uses log link function by default:
% PFC trialresps(y), X = constant and CA1 trialResps, (log(u)) = Xb

% X will be the same for all PFC cells
Xmat = [];
for i=1:size(cellsi,1)
    eval(['currResps = ripplemodi{',num2str(i),'};']);
    Xmat = [Xmat, currResps]; % Rows are observations, Columns are triaslResps for current CA1 cell
end

allmodelfits = []; allmodelb = []; allmodelp=[]; allidxs = []; 
nsig = []; nsigpos = []; nsigneg = []; fracsig =[]; fracsigpos=[]; fracsigneg=[];
% For each PFC cell
for i=1:size(cellsp,1)
    cellsp(i,:);
    eval(['y = ripplemodp{',num2str(i),'};']);
    [b, ~, stats] = glmfit(Xmat,y,'poisson');
    allmodelfits(i).stats = stats;
    allmodelb = [allmodelb; b(2:end)]; % Coefficients. Same as stats.beta
    allmodelp = [allmodelp; stats.p(2:end)];  
    
    % Save index as [day epoch CA1tet CA1cell PFCtet PFCcell]
    PFCdet = repmat([cellsp(i,3), cellsp(i,4)], nCA1cells,1 );
    curridxs = [cellsi,PFCdet];
    allidxs = [allidxs; curridxs];   
    
    % For each PFC neuron, what fraction of CA1 cells were significantly predictive - positive or negative?
    currsig = find(stats.p(2:end) < 0.05);
    nsig(i) = length(currsig);
    fracsig(i) = nsig(i)./nCA1cells;
    bsig = b(currsig+1);
    nsigpos(i) = length(find(bsig>0)); fracsigpos(i) = nsigpos(i)./nCA1cells;
    nsigneg(i) = length(find(bsig<0)); fracsigneg(i) = nsigneg(i)./nCA1cells;
    
end

% Get pairwise correlations between CA1 and PFC - to compare to model coeffs
% Skip getting significance from shuffling for now.
corridxs = []; rcorr = []; pcorr = []; nsimul = [];
for pp=1:size(cellsp,1)
     eval(['y = ripplemodp{',num2str(pp),'};']); % PFC cell
     for ii=1:size(cellsi,1)
         eval(['x = ripplemodi{',num2str(ii),'};']); % CA1 cell
         [r, p] = corrcoef(x,y);
         rcorr = [rcorr; r(1,2)];
         pcorr = [pcorr; p(1,2)];
         corridxs = [corridxs; day epoch cellsi(ii,3) cellsi(ii,4) cellsp(pp,3) cellsp(pp,4)];
         
         % Get number of "trials/ripples" with co-occurences as well
         nsimul = [nsimul; length(find((x~=0) & (y~=0)))];
     end
end


% Skip bad fits, and corrlns. where no. of co-occurences are <10
%rem = find( ((allmodelb>1) || (allmodelb<-1)) && allmodelp > 0.99); % Corresponding p will be 1
rem2 = find(nsimul<10);
allrem = rem2;
%allrem = union(rem, rem2);
allmodelb(allrem)=[]; allmodelp(allrem)=[]; rcorr(allrem)=[]; pcorr(allrem)=[]; allidxs(allrem,:)=[]; corridxs(allrem,:)=[];
sigglm = find(allmodelp < 0.05);
sigcorr = find(pcorr < 0.05);

% Sig GLM vs Sig Corr
glmvec = zeros(size(allmodelb)); glmvec(sigglm)=1;
corrvec = zeros(size(allmodelb)); corrvec(sigcorr)=1;
[rvec,pvec] = corrcoef(glmvec,corrvec)


figure; hold on;
plot(rcorr, allmodelb, 'k.','MarkerSize',24);
plot(rcorr(sigglm), allmodelb(sigglm), 'r.','MarkerSize',24);
plot(rcorr(sigcorr), allmodelb(sigcorr), 'bo','MarkerSize',20);
title(sprintf('GLM fits vs Corr Coeff'),'FontSize',24,'Fontweight','normal');
xlabel(['Corr Coeff'],'FontSize',24,'Fontweight','normal');
ylabel(['GLM coeffs'],'FontSize',24,'Fontweight','normal');
legend('All Pairs','Sig GLM','Sig Corr');

[rall,pall] = corrcoef(rcorr,allmodelb) 
[rglm,pglm] = corrcoef(rcorr(sigglm),allmodelb(sigglm)) 
[rc,pc] = corrcoef(rcorr(sigcorr),allmodelb(sigcorr))

figure; hold on;
plot(abs(rcorr), abs(allmodelb), 'k.','MarkerSize',24);
plot(abs(rcorr(sigglm)), abs(allmodelb(sigglm)), 'r.','MarkerSize',24);
plot(abs(rcorr(sigcorr)), abs(allmodelb(sigcorr)), 'bo','MarkerSize',20);
title(sprintf('ABS GLM fits vs Corr Coeff'),'FontSize',24,'Fontweight','normal');
xlabel(['Corr Coeff'],'FontSize',24,'Fontweight','normal');
ylabel(['GLM coeffs'],'FontSize',24,'Fontweight','normal');
legend('All Pairs','Sig GLM','Sig Corr');


