% clear 


MFS=facsseq();
MFS.load('MFS_combinedcounts.csv');
%% organize sequence counts into binning experiments
ids=1:24;
subrunnames={'theo1-','theo2-','theo1+','theo2+'};
subrunindices={ids(6:-1:1),ids(12:-1:7),ids(18:-1:13),ids(24:-1:19)};
MFS.assignsubrun(subrunnames,subrunindices,[1:5 6]);

% rescale sequence read counts to intended read counts
cellssorted=[249946 451904 257676 89239 31936 7610 288793 390949 259622 99416 33511 5426 335452 436197 202463 75533 38997 7776 352298 426368 192559 72889 36995 7119 382763 481121 312594 91994 32122 8693 429470 481218 288812 83868 29775 7914 590158 368955 152677 48145 21733 6739 579887 343241 142412 45584 21561 6860];
MFS.scalebincounts(cellssorted,[1:5 6.5]);

% find control ribozymes
MFS.findctrls();
%% find validated theophylline sequences
val.seqs={
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCCATAAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAAAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAGAAAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCCAGAAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAGGAAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAGAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCATTCAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCCTGAAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCCGGCAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCCTTGAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCTGGTAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCGAAAGGGACGAAACAGC',
                    };
             val.names={'CAUAA','AAAAA','AGAAA','CAGAA','AGGAA','AAAGA','AUUCA','CUGAA','CGGCA','CUUGA','UGGUA','GAAAG'};

MFS.findvalids(val);



%% analyze theophylline library
minreadcount=100;

hasTheoaptamer=~cellfun('isempty',regexp(MFS.seqs,'^GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCC[A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G][A|T|C|G]?GGACGAAACAGC$'));
gc=ones(length(MFS.origbincounts(:,1)),1);

for k=1:4
        enoughreads=(sum(MFS.subrun(k).origbincounts,2)>minreadcount);
    gc=gc & enoughreads;
    medcount(k)=median((sum(MFS.subrun(k).origbincounts(enoughreads'&hasTheoaptamer,:),2)));

end


goodrbz=~cellfun('isempty',regexp(MFS.seqs,'GCTGTCACCGG[A|T|C|G]+CCGGTCTGATGAGTC[A|T|C|G]+GACGAAACAGC'));
theofilter=gc'&hasTheoaptamer;


%% how many seqs
sum(hasTheoaptamer)

%% analyze theophylline library
alpha=1e-125;
filters=theofilter;
usescaledbincounts=false;
mapcontrols=true;
switches=false;
setfig('theo replicates');clf
subplot(1,2,1)
id1=1;
id2=2;
theo_minus_reps=MFS.pairplot(id1,id2,filters,usescaledbincounts,minreadcount,switches,mapcontrols,alpha);
xlabel('\mu replicate 1 (-ligand)')
ylabel('\mu replicate 2 (-ligand)')
theo_minus_reps=MFS.findVYBmusbyctrls(theo_minus_reps);

subplot(1,2,2)
hold on
id1=3;
id2=4;
theo_plus_reps=MFS.pairplot(id1,id2,filters,usescaledbincounts,minreadcount,switches,mapcontrols,alpha);
xlabel('\mu replicate 1 (+ligand)')
ylabel('\mu replicate 2 (+ligand)')
theo_plus_reps=MFS.findVYBmusbyctrls(theo_plus_reps);


%%
setfig('theo switches combined');clf
id1=[1 2];
id2=[3 4];
switches=true;

theo_mean_switch=MFS.pairplot(id1,id2,filters,usescaledbincounts,minreadcount,switches,mapcontrols,alpha);
xlabel('\mu (-ligand)')
ylabel('\mu (+ligand)')

theo_mean_switch.ctrlcounts1=[MFS.origbincounts(MFS.ctrls.ids,MFS.subrun(1).idx) MFS.origbincounts(MFS.ctrls.ids,MFS.subrun(2).idx)];
theo_mean_switch.ctrlcounts2=[MFS.origbincounts(MFS.ctrls.ids,MFS.subrun(3).idx) MFS.origbincounts(MFS.ctrls.ids,MFS.subrun(4).idx)];
theo_mean_switch.valcounts1=[MFS.origbincounts(MFS.val.ids,MFS.subrun(1).idx) MFS.origbincounts(MFS.val.ids,MFS.subrun(2).idx)];
theo_mean_switch.valcounts2=[MFS.origbincounts(MFS.val.ids,MFS.subrun(3).idx) MFS.origbincounts(MFS.val.ids,MFS.subrun(4).idx)];
theo_mean_switch.ctrlseqs=MFS.ctrls.seqs;
theo_mean_switch.valseqs=MFS.val.seqs;
% obtain VYB mus
theo_mean_switch=MFS.findVYBmusbyctrls(theo_mean_switch);


theo_mean_switch.minus=theo_minus_reps;
theo_mean_switch.plus=theo_plus_reps;

theo_mean_switch.fold=10.^(theo_mean_switch.VYBmus2-theo_mean_switch.VYBmus1);
theo_mean_switch.fold1=10.^(theo_mean_switch.plus.VYBmus1-theo_mean_switch.minus.VYBmus1);
theo_mean_switch.fold2=10.^(theo_mean_switch.plus.VYBmus2-theo_mean_switch.minus.VYBmus2);


setfig('VYB scaled theo');clf
scatter(theo_mean_switch.VYBmus1,theo_mean_switch.VYBmus2,'filled','markerfacealpha',0.2,'markerfacecolor',[0 0.4 1])
hold on
plot(theo_mean_switch.ctrlVYBmus1,theo_mean_switch.ctrlVYBmus2,'+','markersize',14,'linewidth',2)
plot(theo_mean_switch.valVYBmus1,theo_mean_switch.valVYBmus2,'o','markersize',14,'linewidth',2)
x=linspace(min(theo_mean_switch.VYBmus1),max(theo_mean_switch.VYBmus1));
plot(x,x,'k--','linewidth',2)

xlabel('mCherry/BFP (-ligand)')
ylabel('mCherry/BFP (+ligand)')

set(gca,'linewidth',2)
set(gca,'fontsize',22)
axis([-1 1.5 -1 1.5])


setfig('VYB scaled theo normalize to sTRSVctl');clf
scatter(theo_mean_switch.VYBmus1-theo_mean_switch.ctrlVYBmus1(2)+2,...
        theo_mean_switch.VYBmus2-theo_mean_switch.ctrlVYBmus2(2)+2,'filled','markerfacealpha',0.2,'markerfacecolor',[0 0.4 1])
hold on
plot(theo_mean_switch.ctrlVYBmus1-theo_mean_switch.ctrlVYBmus1(2)+2,theo_mean_switch.ctrlVYBmus2-theo_mean_switch.ctrlVYBmus2(2)+2,'+','markersize',14,'linewidth',2)
plot(theo_mean_switch.valVYBmus1-theo_mean_switch.ctrlVYBmus1(2)+2,theo_mean_switch.valVYBmus2-theo_mean_switch.ctrlVYBmus2(2)+2,'o','markersize',14,'linewidth',2)
x=linspace(min(theo_mean_switch.VYBmus1-theo_mean_switch.ctrlVYBmus1(2)+2),max(theo_mean_switch.VYBmus1-theo_mean_switch.ctrlVYBmus1(2)+2));
plot(x,x,'k--','linewidth',2)

xlabel('mCherry/BFP (-ligand)')
ylabel('mCherry/BFP (+ligand)')

set(gca,'linewidth',2)
set(gca,'fontsize',24)
axis([0 3 0 3])
box on
legend('Library sequence','Control ribozymes','Validation candidates','1:1','location','southeast')
%%
setfig('theo volcano');clf
MFS.volcanoplot(theo_mean_switch)
xlabel('AR (mCherry/BFP)')
ylabel('-log_1_0(p-value)')


switches_MFSVc=theo_mean_switch.seqs(find(theo_mean_switch.ttestH));
box on
set(gca,'YScale','log')

yL=get(gca,'YLim');
hold on
plot(log10(2)*ones(1,100),linspace(yL(1),yL(2)),'k:','linewidth',2)

legend('Library sequence','Putative switches','AR = 2','location','northwest')
set(gca,'fontsize',26)


%%
setfig('theo replicates concatenated normalize to sTRSVctl');clf
rep1=[theo_minus_reps.VYBmus1-theo_minus_reps.ctrlVYBmus1(2)+2 theo_plus_reps.VYBmus1-theo_plus_reps.ctrlVYBmus1(2)+2];
rep2=[theo_minus_reps.VYBmus2-theo_minus_reps.ctrlVYBmus2(2)+2 theo_plus_reps.VYBmus2-theo_plus_reps.ctrlVYBmus2(2)+2];
crep1=[theo_minus_reps.ctrlVYBmus1-theo_minus_reps.ctrlVYBmus1(2)+2 theo_plus_reps.ctrlVYBmus1-theo_plus_reps.ctrlVYBmus1(2)+2];
crep2=[theo_minus_reps.ctrlVYBmus2-theo_minus_reps.ctrlVYBmus2(2)+2 theo_plus_reps.ctrlVYBmus2-theo_plus_reps.ctrlVYBmus2(2)+2];
vrep1=[theo_minus_reps.valVYBmus1-theo_minus_reps.ctrlVYBmus1(2)+2 theo_plus_reps.valVYBmus1-theo_plus_reps.ctrlVYBmus1(2)+2];
vrep2=[theo_minus_reps.valVYBmus2-theo_minus_reps.ctrlVYBmus2(2)+2 theo_plus_reps.valVYBmus2-theo_plus_reps.ctrlVYBmus2(2)+2];

scatter(rep1,rep2,70,'filled','markerfacealpha',0.2,'markerfacecolor',[0 0.4 1])
hold on
plot(crep1,crep2,'+','linewidth',2,'markersize',14)
plot(vrep1,vrep2,'o','linewidth',2,'markersize',14)

mdl=fitlm(rep1,rep2);
rsq=mdl.Rsquared.Adjusted;
fitcoeffs=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2)];
x=linspace(min(rep1),max(rep1));
y=(fitcoeffs(1)+fitcoeffs(2)*x);
plot(x,y,'k--','linewidth',2)

text(0.5,2,sprintf('R^2=%0.2f',rsq),'fontsize',20)

xlabel('mCherry/BFP replicate 1')
ylabel('mCherry/BFP replicate 2')
set(gca,'fontsize',22)
set(gca,'linewidth',2)

box on
axis([0 3 0 3])

legend('Library sequence','Control ribozymes','Validation candidates','1:1','location','southeast')



%% analyze cdG

%% END












