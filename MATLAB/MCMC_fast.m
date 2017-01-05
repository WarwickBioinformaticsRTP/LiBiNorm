clear all;

close all;
% minus is the 3' end and plus is the 5' end


%% load data as a vector of reads
 

%fid = fopen('WoldSmartSeq'); %all6
%fid = fopen('Hebenstreit_A.plus.minus'); %all6 A and B
%fid = fopen('Tang31March2015'); %all6
%fid = fopen('Deng31March2015'); %all6
%fid = fopen('Th2-poly-A-tagging.NoA.plus.minus'); %all6
%fid = fopen('BhargavaSmartseq31March2015'); %all6
%fid = fopen('CEL-seq.NoA.plus.minus');
%fid = fopen('Smartseq2'); %all6
%fid = fopen('Bhargava_SMF_25pg.NoA.plus.minus'); %all6
%fid = fopen('Bhargava_AA_25pg.NoA.plus.minus'); %all6
%fid = fopen('Bhargava_AA_1ng.NoA.plus.minus'); %all6
%fid = fopen('SasagawaQuartzSeq'); %all6
%fid = fopen('SRR557798.NoA.plus.minus'); %all6
%fid = fopen('r42ss72.NoA.plus.minus'); %all6
%fid = fopen('r25ss42.NoA.plus.minus'); %all6
%fid = fopen('r42ss42.NoA.plus.minus'); %all6
%fid = fopen('r25ss72.NoA.plus.minus'); %all6
%fid = fopen('ERR489030.NoA.plus.minus'); %all6
%fid = fopen('ERR488955.NoA.plus.minus'); %all6
%fid = fopen('ERR488972.NoA.plus.minus'); %all6
%fid = fopen('ERR489013.NoA.plus.minus'); %all6
%fid = fopen('ERR488982.NoA.plus.minus'); %all6
%fid = fopen('R25p72.NoA.plus.minus'); %all6
%fid = fopen('R42p42.NoA.plus.minus'); %all6
%fid = fopen('R42p72.NoA.plus.minus'); %all6
%fid = fopen('R25p42.NoA.plus.minus'); %all6
%fid = fopen('R25ss72E.NoA.plus.minus'); %all6
%fid = fopen('R42ss42E.NoA.plus.minus'); %all6
%fid = fopen('R42ss72E.NoA.plus.minus'); %all6
%fid = fopen('R25ss42E.NoA.plus.minus'); %all6
%fid = fopen('R25ss72EII.NoA.plus.minus'); %all6
%fid = fopen('R42ss42EII.NoA.plus.minus'); %all6
%fid = fopen('R42ss72EII.NoA.plus.minus'); %all6
%fid = fopen('R25ss42EII.NoA.plus.minus'); %all6
%fid = fopen('R25ss42.17Dec2015'); %all6
%fid = fopen('R42ss72.17Dec2015'); %all6
%fid = fopen('R25ss72.17Dec2015'); %all6
%fid = fopen('R42ss42.17Dec2015'); %all6
%fid = fopen('R10SS72.ExtremeTemp.17Dec2015'); %all6
%fid = fopen('R25SS25.ExtremeTemp.17Dec2015'); %all6
%fid = fopen('R25SS72.ExtremeTemp.17Dec2015'); %all6
%fid = fopen('R42SS25.ExtremeTemp.17Dec2015'); %all6
%fid = fopen('R42SS42.ExtremeTemp.17Dec2015'); %all6
%fid = fopen('R10SS72R2.ExtremeTemp2.17Dec2015'); %all6
fid = fopen('R25SS25R2.ExtremeTemp2.17Dec2015'); %all6




tline = fgetl(fid);


 i = 1;
 
 while ischar(tline)
   % disp(tline)
    [token, remain] = strtok(tline);
    [token, remain] = strtok(remain);
    gene_l(i) = str2num(token);
    [token, remain] = strtok(remain);
    plus{i} = str2num(remain);
    
    tline = fgetl(fid);
    [token, remain] = strtok(tline);
    [token, remain] = strtok(remain);
    [token, remain] = strtok(remain);
    minus{i} = str2num(remain);
    
   tline = fgetl(fid); 
    
   i = i+1;
 end
 
 
% global gene_l plus minus; 
 
 for i = 1:length(gene_l)
    plus{i} = plus{i}(find(plus{i}>0));
    plus{i} = plus{i}(find(plus{i}<gene_l(i)));
    minus{i} = minus{i}(find(minus{i}>0)); 
    minus{i} = minus{i}(find(minus{i}<gene_l(i)));
 end

% putting all the reads in one vector for fast LL evaluations
 
 MaxRead = 100;
 Reads = [];
 Reads_l = [];
 freq_l = [];
 
 bins = linspace(-500, 10000, 22);
 bins(1) = 0;
 bins(2) = 300;
 bins(23) = 11000;
 bins(24) = 12000;
 bins(25) = 15000;
 bins(26) = 30000;
 
% bins = logspace(2, 4.5, 20);
 
 [freq,ind] = histc(gene_l, bins);
 
 freq(1) = freq(1)*2;
 freq(2) = freq(2)*2;
 freq(22) = freq(22)/2;
 freq(23) = freq(23)/2;
 freq(24) = freq(24)/6;
 freq(25) = freq(25)/6;
 
 for i = 1:length(gene_l)
    if (length(plus{i}) > 0 && length(plus{i}) <= MaxRead)
      Reads = [Reads, plus{i}];
      Reads_l = [Reads_l, ones(size(plus{i}))*gene_l(i)]; 
      freq_l = [freq_l, ones(size(plus{i}))*freq(ind(i))];
  elseif (length(plus{i}) > MaxRead)
      perm = randperm(length(plus{i}),MaxRead);
      Reads = [Reads, plus{i}(perm)];
      Reads_l = [Reads_l, ones(1, MaxRead)*gene_l(i)]; 
      freq_l = [freq_l, ones(1, MaxRead)*freq(ind(i))];
  end
  if (length(minus{i}) > 0 && length(minus{i}) <= MaxRead)
      Reads = [Reads, minus{i}];
      Reads_l = [Reads_l, ones(size(minus{i}))*gene_l(i)];
      freq_l = [freq_l, ones(size(minus{i}))*freq(ind(i))];
  elseif (length(minus{i}) > MaxRead)
      perm = randperm(length(minus{i}),MaxRead);
      Reads = [Reads, minus{i}(perm)];
      Reads_l = [Reads_l, ones(1, MaxRead)*gene_l(i)]; 
      freq_l = [freq_l, ones(1, MaxRead)*freq(ind(i))];
   end
 end
 
 data(1, :) = Reads;
 data(2, :) = Reads_l;
 data(3, :) = freq_l;
 
%% MCMC
 
 % select the method used

method = 'mh'; % 'mh','am','dr', or 'dram', see below

JumpSize = 0.01;
options.nsimu = 2000;
Nruns = 100;

mse = 1;
model.sigma2 = mse;
%model.S20 = 10;    %  prior for sigma2
%model.N0 = 10;

switch method
 case 'mh'

   drscale  = 0;
   adaptint = 0;
 case 'dr'
  drscale  = 2; 
  adaptint = 0;
 case 'am'
  drscale  = 0; 
  adaptint = 100;
 case 'dram'
  options.drscale  = 2; 
  options.adaptint = 100;
end

options.updatesigma = 0;

options.method = method;


for ii = 1:6
    
Model = ii; % 1:ModelA; 2:ModelB, 3:Smart(Deng), 4:Tang, 5:Dan, 6:ModelB+Tang
RejectionRate(ii) = 0;

switch Model
    case 4
    model.ssfun = @FLL_Tang;
    case 3
    model.ssfun = @FLL_Deng;
    case 5
    model.ssfun = @FLL_Dan;
    case 1
    model.ssfun = @FLL_ModelA;    
    case 2
    model.ssfun = @FLL_ModelB;
    case 6
    model.ssfun = @FFL_ModelE;
end

clear Chain;

for kk = 1:Nruns
ii
kk
% create input arguments for the dramrun function
p0 = [rand(1)*3, rand(1)*3, rand(1)*4-5, rand(1)*4-5, rand(1)];

switch Model
    case {2, 4, 5}
    options.qcov = eye(4)*JumpSize;
    
    params = {
    {'d', p0(1), -1 , 2};    % average length of fragments
    {'h',  p0(2), 0 , 3};   % the minimum length of fragmenation
    {'t1', p0(3), -5 , -1};    % theta1
    {'t2', p0(4), -5, -1};  % theta2
   % {'sig', p0(5), 0, 3}; % sigma
    };
    case 3
    options.qcov = eye(3)*JumpSize;
    
    params = {
    {'d', p0(1), -1 , 2};    % average length of fragments
    {'h',  p0(2), 0 , 3};   % the minimum length of fragmenation
  %  {'t1', p0(3), -5 , -1};    % theta1
    {'t2', p0(4), -5, -1}; 
    };
    case 1
    options.qcov = eye(2)*JumpSize;
    
    params = {
    {'d', p0(1), -1 , 2};    % average length of fragments
    {'h',  p0(2), 0 , 3};   % the minimum length of fragmenation
  %  {'t1', p0(3), p0(3)/50 , p0(3)*50};    % theta1
  %  {'t2', p0(4), p0(4)/50, p0(4)*50};  % theta2
    };
    case 6
    options.qcov = eye(6)*JumpSize; 
    
    params = {
    {'d', p0(1), -1 , 2};    % average length of fragments
    {'h',  p0(2), 0 , 3};   % the minimum length of fragmenation
    {'t1', p0(3), -5 , -1};    % theta1
    {'t2', p0(4), -5, -1};  % theta2
    {'a', p0(5), 0, 1}; % alpha strength of model B
    };
end
                
                  
[results,chain, s2chain, sschain] = mcmcrun(model,data,params,options);


Chain(kk, :) = chain(end, :);
SSchain(ii, kk) =  sschain(end);
RejectionRate(ii) = RejectionRate(ii)+ results.rejected;

end

%save Tang_onlyFrag_Dram_TangFit results chain s2chain sschain;

RejectionRate(ii) = RejectionRate(ii)/Nruns;

%load Bhargava_onlyFrag_DRAM_TangFit.mat
%options.nsimu = 1800;
burn = options.nsimu/2;

%param = mean(10.^chain(burn:options.nsimu, :))
%s_param = std(10.^chain(burn:options.nsimu, :))
switch ii
    case 1
    param(ii, 1) = median(10.^Chain(:, 1));
    param(ii, 2) = median(10.^Chain(:, 2));
    param(ii, 3) = 0;
    param(ii, 4) = 0;
    param(ii, 5) = 0;

    s_param(ii,1) = mad(10.^Chain(:,1));
    s_param(ii,2) = mad(10.^Chain(:,2));
    s_param(ii,3) = 0;
    s_param(ii,4) = 0;
    s_param(ii, 4) = 0;
    case {2, 4, 5}    
    
    param(ii, 1) = median(10.^Chain(:, 1));
    param(ii, 2) = median(10.^Chain(:, 2));
    param(ii, 3) = median(10.^Chain(:, 3));
    param(ii, 4) = median(10.^Chain(:, 4));
    param(ii, 5) = 0;


    s_param(ii,1) = mad(10.^Chain(:,1));
    s_param(ii,2) = mad(10.^Chain(:,2));
    s_param(ii,3) = mad(10.^Chain(:,3));
    s_param(ii,4) = mad(10.^Chain(:,4));
    s_param(ii, 5) = 0;

    
    case 3
        
    param(ii, 1) = median(10.^Chain(:, 1));
    param(ii, 2) = median(10.^Chain(:, 2));
    param(ii, 3) = 0;
    param(ii, 4) = median(10.^Chain(:, 3));
    param(ii, 5) = 0;


    s_param(ii,1) = mad(10.^Chain(:,1));
    s_param(ii,2) = mad(10.^Chain(:,2));
    s_param(ii,3) = 0;
    s_param(ii,4) = mad(10.^Chain(:,3));  
    s_param(ii, 5) = 0;
    
    case 6
    
    param(ii, 1) = median(10.^Chain(:, 1));
    param(ii, 2) = median(10.^Chain(:, 2));
    param(ii, 3) = median(10.^Chain(:, 3));
    param(ii, 4) = median(10.^Chain(:, 4));
    param(ii, 5) = median(Chain(:, 5));


    s_param(ii,1) = mad(10.^Chain(:,1));
    s_param(ii,2) = mad(10.^Chain(:,2));
    s_param(ii,3) = mad(10.^Chain(:,3));
    s_param(ii,4) = mad(10.^Chain(:,4));
    s_param(ii, 5) = mad(Chain(:, 5));
        

end

%param = [5   50    0.03    0.01]
d = param(ii,1);
h = param(ii,2);
t1 = param(ii,3);
t2 = param(ii,4);
a = param(ii, 5);

%LL_median(ii) = median(SSchain)





%BayesF = (options.nsimu - burn) / sum(exp(sschain(burn:options.nsimu)))

figure(5*ii-4);

mcmcplot(Chain, [], results, 'hist');

figure(5*ii-3); 
mcmcplot(chain,[],results,'chainpanel');

figure(5*ii-2);

hist(SSchain(ii,:));

figure(5*ii-1);
%mcmcplot(chain(burn:options.nsimu, :),[],results,'pairs');
mcmcplot(Chain,[],results,'pairs');


geneL = [500 1000 2000 4000 8000];
color = ['b', 'r', 'y', 'k', 'g'];

figure(5*ii); 
subplot(2, 1, 1);


for j = 1:length(geneL)
    Total_reads = [];
    for i = 1:length(gene_l)
        if (abs(gene_l(i) - geneL(j))<0.1*geneL(j))
            Total_reads = cat(2, Total_reads, minus{i}/gene_l(i), plus{i}/gene_l(i));
        end
    end
    
    [f_L,xi_L] = ksdensity(Total_reads); 
    xi = linspace(0, geneL(j), 100);
    
    switch Model
        case 4
            LL_L =  exp(Tang_logL_minus(xi, geneL(j), d, h, t1, t2))*geneL(j);
        case 3
            LL_L = exp(Deng_logL_minus(xi, geneL(j), d, h, t1, t2))*geneL(j);
        case 5
            LL_L = exp(Dan_logL_plus(xi, geneL(j), d, h, t1, t2))*geneL(j);
        case 1
            x=xi;
            l = geneL(j);
            LL_L =  ((x> h).*(x < l-h) + 1/d)./((2*h<l).*(l-2*h) + l/d)*geneL(j);
        case 2
            x=xi;
            l = geneL(j);
            LL_L = ((x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l-h+x)) + t2*exp(-l*(t1+t2))) + ...
                    1./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l+x)) + t2*exp(-l*(t1+t2)))/d)./ ...
               ((2*h<l).*(t1.*(exp(-2*h*(t1+t2))-exp(-l*(t1+t2)))+t2*(t1+t2).*(l-2*h).*exp(-l*(t1+t2)))/(t1 + t2)^2 + ...
                exp(-l*(t1+t2)).*(l.*t2^2+t1*(exp(l*(t1 + t2))+l.*t2-1))/(t1 + t2)^2/d)*geneL(j);
        case 6
            x=xi;
            l = geneL(j);
            LL_L = (a*( (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l-h+x)) + t2*exp(-l*(t1+t2))) + ...
         1./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l+x)) + t2*exp(-l*(t1+t2)))/d) + ...
         (1-a)*( (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2.*exp(-t1*l-t2*(x+h))) + ...
         1./(t1+ t2) .* (t1.*exp(-t1*(l-x)) + t2.*exp(-t1*l-t2*(x)))/d))./...
         (a*( (2*h<l).*(t1.*(exp(-2*h*(t1+t2))-exp(-l.*(t1+t2)))+t2*(t1+t2).*(l-2*h).*exp(-l.*(t1+t2)))/(t1 + t2)^2 + ...
        (exp(-l.*(t1+t2)).*(l.*t2^2+l.*t2*t1-t1)+t1)/(t1 + t2)^2/d) + ...
        (1-a)*( (2*h<l).*(exp(-2*h*(t1 + t2)) - exp(-l.*(t1 + t2)))/(t1 + t2) + ...
        (1-exp(-l.*(t1 + t2)))/(t1 + t2)/d))*geneL(j);
    end
    
    plot(xi/geneL(j), LL_L, ['--',color(j)], 'LineWidth', 2);
    hold on;
    plot(xi_L, f_L, color(j), 'LineWidth', 2);
    hold on;
end

xlim([-0.1 1.1]);
legend('500', '500', '1000', '1000','2000','2000', '4000','4000', '8000', '8000');


subplot(2, 1, 2); 

l = linspace(100, 10000, 101);

if (Model ==4)
    for i=1:101
        
        if (2*h<l(i))
  
               norm(i) =  (exp(-2*h*(t1 + t2)) - exp(-l(i)*(t1 + t2)))/(t1 + t2) + (1-exp(-l(i)*(t1 + t2)))/(t1 + t2)/d;
        else
   
               norm(i) =  (1-exp(-l(i)*(t1 + t2)))/(t1 + t2)/d;
end

    end
    
end

if (Model ==3)
    for i=1:101
        
        if (2*h<l(i))
  
               norm(i) =  (exp(-2*h*t2-l(i)*t1) - exp(-l(i)*(t1 + t2)))/t2 + (exp(-l(i)*t1)-exp(-l(i)*(t1 + t2)))/t2/d;
        else
   
               norm(i) =  (exp(-l(i)*t1)-exp(-l(i)*(t1 + t2)))/t2/d;
end

    end
end

if (Model ==5)
    for i=1:101
        
        if (2*h<l(i))
  
               norm(i) =  (exp(-l(i)*t1 - 2*h*t2)*(t1 + t2)^2 - exp(-l(i)*(t1 + t2))*t1^2 + t1*t2*exp(-2*h*(t1 + t2))*(l(i)*t2 -2*h*t1 -2*h*t2+l(i)*t1 - t2/t1 - 2))/(t1 + t2)^2/t1^2/t2 + ...
         (l(i)-1/(t1 + t2) - 1/t1 - t1/t2/(t1+t2)*exp(-l(i)*(t1 + t2))+(t1 + t2)/t1/t2*exp(-l(i)*t1))/(t1 + t2)/t1/d;
        else
   
               norm(i) =  (l(i)-1/(t1 + t2) - 1/t1 - t1/t2/(t1+t2)*exp(-l(i)*(t1 + t2))+(t1 + t2)/t1/t2*exp(-l(i)*t1))/(t1 + t2)/t1/d;
end

    end
end

if (Model == 1)
    norm = (2*h<l).*(l-2*h) + l/d;
end
    
if (Model == 2)    
    norm = ((2*h<l).*(t1.*(exp(-2*h*(t1+t2))-exp(-l.*(t1+t2)))+t2*(t1+t2).*(l-2*h).*exp(-l.*(t1+t2)))/(t1 + t2)^2 + ...
                exp(-l.*(t1+t2)).*(l.*t2^2+t1*(exp(l.*(t1 + t2))+l.*t2-1))/(t1 + t2)^2/d);
end
    
if (Model == 6)
    norm = a*( (2*h<l).*(t1.*(exp(-2*h*(t1+t2))-exp(-l.*(t1+t2)))+t2*(t1+t2).*(l-2*h).*exp(-l.*(t1+t2)))/(t1 + t2)^2 + ...
        (exp(-l.*(t1+t2)).*(l.*t2^2+l.*t2*t1-t1)+t1)/(t1 + t2)^2/d) + ...
        (1-a)*( (2*h<l).*(exp(-2*h*(t1 + t2)) - exp(-l.*(t1 + t2)))/(t1 + t2) + ...
        (1-exp(-l.*(t1 + t2)))/(t1 + t2)/d);
end

plot(l, norm, 'LineWidth', 2);
xlabel('gene length');
ylabel('normalization factor');


end

%%
for i = 1:6
    n_nan(i) = length(find(isnan(SSchain(i, :))));
    n_inf(i) = length(find(isinf(SSchain(i, :))));
    vec = SSchain(i, find((~isnan(SSchain(i, :))&(~isinf(SSchain(i, :))))));
    LL_min(i) = min(vec);
    LL_median(i) = median(vec);
    LL_mad(i) = mad(vec);
end
LL_ratio = LL_median-min(LL_median);


  

save R25SS25R2.ExtremeTemp2.17Dec_all6 Chain SSchain param s_param LL_min LL_median LL_mad LL_ratio n_nan n_inf RejectionRate;