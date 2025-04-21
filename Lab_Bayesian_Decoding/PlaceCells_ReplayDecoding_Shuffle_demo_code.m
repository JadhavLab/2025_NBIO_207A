% let's first calculate how sequential the real sequence is
% calculate weighted correlation
mloc = sum(sum(pMat_all.*x_position))/sum(sum(pMat_all));
mt = sum(sum(pMat_all.*timevec))/sum(sum(pMat_all));
dloc = x_position-mloc;
dt = timevec - mt;
cov_loc_t = sum(sum(pMat_all.*dloc.*dt))/sum(sum(pMat_all));
cov_loc = sum(sum(pMat_all.*(dloc.^2)))/sum(sum(pMat_all));
cov_t = sum(sum(pMat_all.*(dt.^2)))/sum(sum(pMat_all));
rvalue = cov_loc_t/sqrt(cov_loc*cov_t);
%%
nonzerobins = find(sum(pMat_all) > 0);% non-zeros bins
nperm = 1000;%1000 times permutation

%-------Shuffling time bins to get the pvalue------%
rvalue_tshuf = zeros(nperm,1);
permbins = nonzerobins;
for iteration = 1:nperm 
     % circularly shift position columns
     pMatShuf = zeros(size(pMat_all));
     shiftVal = randperm(length(permbins)); %permute the temporal bins
     pMatShuf(:,permbins(shiftVal)) = pMat_all(:,permbins);
   
     % calculate weighted correlation of the shuffled data
     mloc = sum(sum(pMatShuf.*x_position))/sum(sum(pMatShuf));
     mt = sum(sum(pMatShuf.*timevec))/sum(sum(pMatShuf));
     dloc = x_position-mloc;
     dt = timevec - mt;
     cov_loc_t = sum(sum(pMatShuf.*dloc.*dt))/sum(sum(pMatShuf));
     cov_loc = sum(sum(pMatShuf.*(dloc.^2)))/sum(sum(pMatShuf));
     cov_t = sum(sum(pMatShuf.*(dt.^2)))/sum(sum(pMatShuf));
     rvalue_tshuf(iteration) = cov_loc_t/sqrt(cov_loc*cov_t);
end
pvalue_tr =sum(rvalue < rvalue_tshuf)/length(rvalue_tshuf); %p-value
figure('Position',[100,100,600,200]),
subplot(121)
imagesc(pMatShuf)
ax = gca;
ax.YDir = 'normal';
title('Temporal shuffle')
colormap(hot)
caxis([0,0.1])%probability 0-0.1
%%
%-------Shuffling positions to get the pvalue------%
rvalue_xshuf = zeros(nperm,1);
permbins = nonzerobins;
shiftVal = randi([1,length(x_position)],nperm,length(permbins)); %shift value for the 1000 permutations
for iteration = 1:nperm 
     % circularly shift position columns
     pMatShuf = zeros(size(pMat_all));
     for t = 1:length(permbins)
         pMatShuf(:,permbins(t)) = circshift(pMat_all(:,permbins(t)),shiftVal(iteration,t));
     end
                
     % calculate weighted correlation of the shuffled data
     mloc = sum(sum(pMatShuf.*x_position))/sum(sum(pMatShuf));
     mt = sum(sum(pMatShuf.*timevec))/sum(sum(pMatShuf));
     dloc = x_position-mloc;
     dt = timevec - mt;
     cov_loc_t = sum(sum(pMatShuf.*dloc.*dt))/sum(sum(pMatShuf));
     cov_loc = sum(sum(pMatShuf.*(dloc.^2)))/sum(sum(pMatShuf));
     cov_t = sum(sum(pMatShuf.*(dt.^2)))/sum(sum(pMatShuf));
     rvalue_xshuf(iteration) = cov_loc_t/sqrt(cov_loc*cov_t);
end

pvalue_xr =sum(rvalue < rvalue_xshuf)/length(rvalue_xshuf);%p-value
subplot(122)
imagesc(pMatShuf)
ax = gca;
ax.YDir = 'normal';
title('Position shuffle')
colormap(hot)
caxis([0,0.1])%probability 0-0.1
%%
% visualization of the shuffling result
figure('Position',[100,100,600,200]),
subplot(121)
[counts,ind]=hist(rvalue_tshuf); %reverse the axis
counts = counts./sum(counts);
bar(ind,counts)
hold on
plot([rvalue,rvalue],[0,max(counts)],'r')
hold off
title(['Temporal shuffle, pval = ',num2str(round(pvalue_tr*100)/100)])

subplot(122)
[counts,ind]=hist(rvalue_xshuf); %reverse the axis
counts = counts./sum(counts);
bar(ind,counts)
hold on
plot([rvalue,rvalue],[0,max(counts)],'r')
hold off
title(['Position shuffle, pval = ',num2str(round(pvalue_xr*100)/100)])
