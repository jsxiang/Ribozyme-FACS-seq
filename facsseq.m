classdef facsseq < handle

  properties

    seqs;
    seqlengths;
    origbincounts;
    scaledbincounts;
    subrun;
    motif;
    seqmotif
    switches;
    runID;
    ctrls;
    refs;
    rounds;
    val;
    
    
    rbzseqs;
      end

  methods
    function obj=facsseq()
    end

    function load(data,file)
        % Backward compatibility with old data.mat files
        if nargin<2
            file='combinedcounts.csv';
        end
        fprintf('Loading data from %s...',file);
        d=csvread(file,0,1);
        fid=fopen(file);

        s=textscan(fid,'%s');
        fclose(fid);

        seqs={};
        for j=1:length(s{1})
            o=regexp(s{1}{j},'(\W*)([A-Z]+)(\W*)','tokens');
            seqs{end+1}=o{1}{2};
        end
            data.origbincounts=d;
            data.seqs=seqs;
            data.seqlengths=cellfun('length',data.seqs);

        fprintf('done\n');
%         fn=setdiff(fieldnames(data),'pts');
%         for i=1:length(fn)
%             data.(fn{i})=data.(fn{i});
%         end
        %       obj.cleanup();
    end
    
    function assignsubrun(data,subrunnames,subrunindices,edges)
        if nargin<4
            edges=1:6;
        end
        for i=1:length(subrunnames)
            data.subrun(i).name=subrunnames{i};
            data.subrun(i).idx=subrunindices{i};
            data.subrun(i).origbincounts=data.origbincounts(:,data.subrun(i).idx);
            
            for j=1:length(data.subrun(i).origbincounts(:,1));
                c1=data.subrun(i).origbincounts(j,:);

                [m,s]=normfit(edges,[],[],c1);
                data.subrun(i).mu(j)=m;
                
            end
        end
    end
    
    function findrefs(data,refs,refnames)
        for i=1:length(refs)
            isref=~cellfun('isempty',regexp(data.seqs,strcat('^',refs{i},'$')));
            for j=1:length(data.subrun)
                data.refs(i).subrun(j).mu=data.subrun(j).mu(isref);
                data.refs(i).subrun(j).mus=data.subrun(j).mus(isref);
                data.refs(i).subrun(j).adjustedmus=data.subrun(j).adjustedmus(isref);
                if ~isnan(data.refs(i).subrun(j).mu)
                    fprintf('%16s\t%4s\t%0.4f\t%0.4f\t%0.4f\n',refnames{i},data.subrun(j).name,data.refs(i).subrun(j).mu,data.refs(i).subrun(j).mus,data.refs(i).subrun(j).adjustedmus);
                end
            end
        end
    end
    
    function findseqlengths(data)
        data.seqlengths=cellfun('length',data.seqs);
    end
    
    function findvalids(data,val)
        data.val=val;
        valids=[];
        for i=1:length(data.val.seqs)
            id=find(~cellfun('isempty',regexp(data.seqs,strcat('^',data.val.seqs{i},'$'))));
            if ~isempty(id)
                valids(i)=id;
            end
        end
        data.val.ids=valids;
    end    
    
    
    function scalebincounts(data,cellssorted,edges)
        if nargin<3
            edges=1:6;
        end
        csmat=zeros(length(data.subrun(1)),length(data.subrun(1).idx));
        sumbccount=csmat;
        for i=1:length(data.subrun)
            csmat(i,:)=cellssorted(data.subrun(i).idx);
            sumbccount(i,:)=sum(data.subrun(i).origbincounts);
            scalefac=repmat(1./sumbccount(i,:).*csmat(i,:),length(data.subrun(i).origbincounts(:,1)),1);
            data.subrun(i).scaledbincounts=data.subrun(i).origbincounts.*scalefac;
            for j=1:length(data.subrun(i).scaledbincounts(:,1));
                c1=data.subrun(i).scaledbincounts(j,:);

                [m,s]=normfit(edges,[],[],c1);
                data.subrun(i).mus(j)=m;
                
            end
            
            rsbc=round(data.subrun(i).scaledbincounts);
            
            ttestX={};
            for j=1:length(rsbc)
                c=ones(1,rsbc(j,1));
                for k=2:length(edges)
                   c=[c ones(1,rsbc(j,k))*k];
                end
                ttestX{end+1}=c;
            end
            
            data.subrun(i).ttestX=ttestX;
            
        end
        setfig('barcode count');clf
        plot(csmat',sumbccount','-o')        

    end
    
    function pairedruns=pairplot(data,id1,id2,filters,usescaledbincounts,minreadcount,switches,mapcontrols,alpha)
        if nargin < 7
            switches=false;
            mapcontrols=true;
            alpha=1e-50;
        end
        if usescaledbincounts
            if length(id1)>1
                mu1mat=[];
                ctrl1mat=[];
                val1mat=[];
                for i=1:length(id1)
                    mu1mat=[mu1mat; data.subrun(id1(i)).mus(filters)];
                    ctrl1mat=[ctrl1mat; data.subrun(id1(i)).mus(data.ctrls.ids)];
                    val1mat=[val1mat; data.subrun(id1(i)).mus(data.val.ids)];
                end
                mu1=mean(mu1mat);
                ctrl1=mean(ctrl1mat);
                val1=mean(val1mat);
            else
            mu1=data.subrun(id1).mus(filters);
            ctrl1=data.subrun(id1).mus(data.ctrls.ids);
            val1=data.subrun(id1).mus(data.val.ids);
            end
            if length(id2)>1
                mu2mat=[];
                ctrl2mat=[];
                val2mat=[];
                for i=1:length(id2)
                    mu2mat=[mu2mat; data.subrun(id2(i)).mus(filters)];
                    ctrl2mat=[ctrl2mat; data.subrun(id2(i)).mus(data.ctrls.ids)];
                    val2mat=[val2mat; data.subrun(id2(i)).mus(data.val.ids)];
                end
                mu2=mean(mu2mat);
                ctrl2=mean(ctrl2mat);
                val2=mean(val2mat);
            else
                
            mu2=data.subrun(id2).mus(filters);
            ctrl2=data.subrun(id2).mus(data.ctrls.ids);
            val2=data.subrun(id2).mus(data.val.ids);
            end
        else
            if length(id1)>1
                mu1mat=[];
                ctrl1mat=[];
                val1mat=[];
                for i=1:length(id1)
                    mu1mat=[mu1mat; data.subrun(id1(i)).mu(filters)];
                    ctrl1mat=[ctrl1mat; data.subrun(id1(i)).mu(data.ctrls.ids)];
                    val1mat=[val1mat; data.subrun(id1(i)).mu(data.val.ids)];
                end
                mu1=mean(mu1mat);
                ctrl1=mean(ctrl1mat);
                val1=mean(val1mat);
            else
            mu1=data.subrun(id1).mu(filters);
            ctrl1=data.subrun(id1).mu(data.ctrls.ids);
            val1=data.subrun(id1).mu(data.val.ids);
            end
            if length(id2)>1
                mu2mat=[];
                ctrl2mat=[];
                val2mat=[];
                for i=1:length(id2)
                    mu2mat=[mu2mat; data.subrun(id2(i)).mu(filters)];
                    ctrl2mat=[ctrl2mat; data.subrun(id2(i)).mu(data.ctrls.ids)];
                    val2mat=[val2mat; data.subrun(id2(i)).mu(data.val.ids)];
                end
                mu2=mean(mu2mat);
                ctrl2=mean(ctrl2mat);
                val2=mean(val2mat);
            else
                
            mu2=data.subrun(id2).mu(filters);
            ctrl2=data.subrun(id2).mu(data.ctrls.ids);
            val2=data.subrun(id2).mu(data.val.ids);
            end
        end
        
        if mapcontrols
            ctrlmdl=fitlm(ctrl2,ctrl1);
            fitcoeffs=[ctrlmdl.Coefficients.Estimate(1) ctrlmdl.Coefficients.Estimate(2)];
            adjustedmu2=fitcoeffs(1)+fitcoeffs(2)*mu2;
            adjustedctrl2=fitcoeffs(1)+fitcoeffs(2)*ctrl2;
            adjustedval2=fitcoeffs(1)+fitcoeffs(2)*val2;
        else
            adjustedmu2=mu2;
            adjustedctrl2=ctrl2;
            adjustedval2=val2;
        end
            
        scatter(mu1,adjustedmu2,'filled','markerfacealpha',0.2)

        
        hold on
        x=linspace(1,6);
        plot(x,x,'k--','linewidth',2)
        plot(ctrl1,adjustedctrl2,'+','markersize',12,'linewidth',2);
        plot(val1,adjustedval2,'o','markersize',12,'linewidth',2);

        legend('FACSseq data','1:1','ribozyme controls','validated hits','location','southeast')
        set(gca,'linewidth',1.5)
        set(gca,'fontsize',18)
        
        mdl=fitlm(mu1,adjustedmu2);
        rsq=mdl.Rsquared.Adjusted;
        if ~switches
            text(2,5,sprintf('R^2 = %0.2f',rsq),'fontsize',14)
        end
        title(sprintf('%d unique sequences (>%d total reads)',sum(filters),minreadcount),'fontsize',14)        
        
        filteri=find(filters);
        pairedruns.ttestH=[];
        pairedruns.ttestP=[];
        for j=1:length(filteri)
            if length(id1)>1
                counts1=[];
                counts2=[];
                for k=1:length(id1)
                    counts1=[counts1 data.subrun(id1(k)).ttestX{filteri(j)}];
                    counts2=[counts2 data.subrun(id2(k)).ttestX{filteri(j)}];
                end
            else
                counts1=data.subrun(id1).ttestX{filteri(j)};
                counts2=data.subrun(id2).ttestX{filteri(j)};
            end
            
            if mapcontrols
                adjustedcounts2=fitcoeffs(1)+fitcoeffs(2)*counts2;
            else
                adjustedcounts2=counts2;
            end
%         alpha=0.0001/sum(filters);
        [H,P]=ttest2(counts1,adjustedcounts2,'alpha',alpha,'tail','both','vartype','unequal');
        pairedruns.ttestH(j)=H;
        pairedruns.ttestP(j)=P;
        end
        
        pairedruns.seqs=data.seqs(filters);
        pairedruns.mu1=mu1;
        pairedruns.mu2=mu2;
        pairedruns.adjustedmu2=adjustedmu2;
        
        pairedruns.ctrl1=ctrl1;
        pairedruns.ctrl2=ctrl2;
        pairedruns.adjustedctrl2=adjustedctrl2;
        pairedruns.val1=val1;
        pairedruns.val2=val2;
        pairedruns.adjustedval2=adjustedval2;
        
        
    end
    function findrunID(data,motifID,motifname)
        for i=1:length(motifID)
            hasmotif=sum(data.origbincounts(:,motifID(i)),2)>0;
            data.runID(i).name=motifname{i};
            data.runID(i).seqs=data.seqs(hasmotif);
            data.runID(i).hasmotif=hasmotif;  

            for j=1:length(data.subrun)
                data.runID(i).subrun(j).name=data.subrun(j).name;
                data.runID(i).subrun(j).idx=data.subrun(j).idx;
                data.runID(i).subrun(j).origbincounts=data.subrun(j).origbincounts(hasmotif,:);
                data.runID(i).subrun(j).mus=data.subrun(j).mus(hasmotif);
%                 data.motif(i).subrun(j).VYBmus=data.subrun(j).VYBmus(hasmotif);
                data.runID(i).subrun(j).scaledbincounts=data.subrun(j).scaledbincounts(hasmotif,:);
                data.runID(i).subrun(j).ttestX=data.subrun(j).ttestX(hasmotif);
            end
        end
    end
    
    function findctrls(data,ctrls)
        if nargin<2
            ctrls.seqs={'GCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGC',
                        'GCTGTCACCGGATGTGCTTTCCGGTACGTGAGGTCCGTGAGGACGAAACAGC',
                        'GCTGTCACCGGATGTTTCCGGTCTGATGAGTCCATAAGGACGAAACAGC',
                        'GCTGTCACCGGATACTTCCGGTCTGATGAGTCCCAGAGGACGAAACAGC',
                        'GCTGTCACCGGATGCATCCGGTCTGATGAGTCCCGCGGGACGAAACAGC'};
            ctrls.names={'sTRSV','sTRSVctl','g814','g833','g862'};
        end
        data.ctrls=ctrls;
        ctrlids=[];
        for i=1:length(data.ctrls.seqs)
            id=find(~cellfun('isempty',regexp(data.seqs,strcat('^',data.ctrls.seqs{i},'$'))));
            if ~isempty(id)
                ctrlids(i)=id;
            end
        end
        data.ctrls.ids=ctrlids;
    end
    
    function findmotif(data,motif,motifname)
        allmotifs=zeros(1,length(data.seqs));
        for i=1:length(motif)
            hasmotif=~cellfun('isempty',regexp(data.seqs,motif{i}));
            
            if isempty(regexp(motif{i},'\W'))
                nohasyet=~hasmotif;
                nohasyeti=find(nohasyet);

                for p=1:length(nohasyeti)
                    seq=data.seqs{nohasyeti(p)};
                    [s,a]=swalign(seq,motif{i});
                    nummismatch=sum(~(a(2,:)=='|'));
                    if nummismatch<3
                        hasmotif(nohasyeti(p))=1;
                    end
                end
            end
            
            if i~=length(motif)
                allmotifs=hasmotif|allmotifs;
                data.motif(i).name=motifname{i};
                data.motif(i).seqs=data.seqs(hasmotif);
                data.motif(i).hasmotif=hasmotif;  
            else
                hasmotif=hasmotif&(~allmotifs);
                data.motif(i).name=motifname(i);
                data.motif(i).seqs=data.seqs(hasmotif);
                data.motif(i).hasmotif=hasmotif;
            end
            for j=1:length(data.subrun)
                data.motif(i).subrun(j).name=data.subrun(j).name;
                data.motif(i).subrun(j).idx=data.subrun(j).idx;
                data.motif(i).subrun(j).origbincounts=data.subrun(j).origbincounts(data.motif(i).hasmotif,:);
                data.motif(i).subrun(j).mus=data.subrun(j).mus(data.motif(i).hasmotif);
%                 data.motif(i).subrun(j).VYBmus=data.subrun(j).VYBmus(hasmotif);
                data.motif(i).subrun(j).scaledbincounts=data.subrun(j).scaledbincounts(data.motif(i).hasmotif,:);
                data.motif(i).subrun(j).ttestX=data.subrun(j).ttestX(data.motif(i).hasmotif);
            end
        end
    end
    
    function adjustMus(data,motifID,subrunsID,usescaledbincounts,method)
        for i=1:length(subrunsID(1,:))
            
            goodstats1=sum(data.motif(motifID).subrun(subrunsID(1,i)).origbincounts,2)>80 & sum(data.motif(motifID).subrun(subrunsID(1,i)).origbincounts,2)<500;
            goodstats2=sum(data.motif(motifID).subrun(subrunsID(2,i)).origbincounts,2)>80 & sum(data.motif(motifID).subrun(subrunsID(2,i)).origbincounts,2)<500;
%             goodstats1=sum(data.motif(motifID).subrun(subrunsID(1,i)).origbincounts,2)>80;
%             goodstats2=sum(data.motif(motifID).subrun(subrunsID(2,i)).origbincounts,2)>80;

            goodstats=goodstats1&goodstats2;
            if strcmp(method,'linear')
            if usescaledbincounts
                preadjustedmus=data.motif(motifID).subrun(subrunsID(1,i)).mus(goodstats);
                refmus=data.motif(motifID).subrun(subrunsID(2,i)).mus(goodstats);
                mdl=fitlm(preadjustedmus,refmus);
            else
                preadjustedmus=data.motif(motifID).subrun(subrunsID(1,i)).mu(goodstats);
                refmus=data.motif(motifID).subrun(subrunsID(2,i)).mu(goodstats);
                mdl=fitlm(preadjustedmus,refmus);
            end
            fitcoeffs=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2)];
            
            
            
            if usescaledbincounts
                adjustedmus=fitcoeffs(1)+fitcoeffs(2)*data.subrun(subrunsID(1,i)).mus;
                data.subrun(subrunsID(1,i)).adjustedmus=adjustedmus;
                data.subrun(subrunsID(2,i)).adjustedmus=data.subrun(subrunsID(2,i)).mus;
                adjustedrefmus= data.subrun(subrunsID(2,i)).adjustedmus;
                
            else
                adjustedmus=fitcoeffs(1)+fitcoeffs(2)*data.subrun(subrunsID(1,i)).mu;
                data.subrun(subrunsID(1,i)).adjustedmus=adjustedmus;
                data.subrun(subrunsID(2,i)).adjustedmus=data.subrun(subrunsID(2,i)).mu;
                adjustedrefmus=data.subrun(subrunsID(2,i)).adjustedmus;
                
            end
            elseif strcmp(method,'median')
                if usescaledbincounts
                preadjustedmus=data.motif(motifID).subrun(subrunsID(1,i)).mus(goodstats);
                refmus=data.motif(motifID).subrun(subrunsID(2,i)).mus(goodstats);
                medianscale=median(preadjustedmus)-median(refmus);
            else
                preadjustedmus=data.motif(motifID).subrun(subrunsID(1,i)).mu(goodstats);
                refmus=data.motif(motifID).subrun(subrunsID(2,i)).mu(goodstats);
                medianscale=median(preadjustedmus)-median(refmus);
            end            
            
            
            if usescaledbincounts
                adjustedmus=data.subrun(subrunsID(1,i)).mus-medianscale;
                data.subrun(subrunsID(1,i)).adjustedmus=adjustedmus;
                data.subrun(subrunsID(2,i)).adjustedmus=data.subrun(subrunsID(2,i)).mus;
                adjustedrefmus= data.subrun(subrunsID(2,i)).adjustedmus;
                
            else
                adjustedmus=data.subrun(subrunsID(1,i)).mu-medianscale;
                data.subrun(subrunsID(1,i)).adjustedmus=adjustedmus;
                data.subrun(subrunsID(2,i)).adjustedmus=data.subrun(subrunsID(2,i)).mu;
                adjustedrefmus=data.subrun(subrunsID(2,i)).adjustedmus;
                
            end
            end

            
            setfig(strcat(data.subrun(subrunsID(1,i)).name,' vs ',data.subrun(subrunsID(2,i)).name,' adjust mu')),clf
            subplot(1,2,1)
            plot(preadjustedmus,refmus,'.');
            x=linspace(min(refmus),max(refmus));
            hold on
            plot(x,x,':','linewidth',2)
            
            if strcmp(method,'linear')
            plot(x,fitcoeffs(1)+fitcoeffs(2)*x,':','linewidth',3)
            elseif strcmp(method,'median')
            plot(x,x-medianscale,':','linewidth',3)
            end            
            xlabel(data.subrun(subrunsID(1,i)).name);
            ylabel(data.subrun(subrunsID(2,i)).name);
            legend('data','1:1','fit','location','northwest')
            
            
            subplot(1,2,2)
            plot(adjustedmus,adjustedrefmus,'.');
            x=linspace(min(adjustedrefmus),max(adjustedrefmus));
            hold on
            plot(x,x,':','linewidth',2)
            if strcmp(method,'linear')
            plot(x,fitcoeffs(1)+fitcoeffs(2)*x,':','linewidth',3)
            elseif strcmp(method,'median')
            plot(x,x-medianscale,':','linewidth',3)
            end
            legend('data','1:1','fit','location','northwest')
            xlabel(data.subrun(subrunsID(1,i)).name)
            ylabel(data.subrun(subrunsID(2,i)).name)
                        
        end
        
        for i=1:length(data.motif)
            for j=1:length(data.motif(i).subrun)
                data.motif(i).subrun(j).adjustedmus=data.subrun(j).adjustedmus(data.motif(i).hasmotif);
            end
        end
    end
    
    function subrun=findVYBmusbyctrls(~,subrun)
%            name    |   mu-    |   mu+   
%            t3hRECK |  0.1210  |  0.1210  
%              sTRSV |  0.3068  |  0.1210  
%               g814 |  0.4380  |  0.1210  
%               g833 |  0.7446  |  0.1210  
%               g862 |  2.2635  |  0.1210  
%         t3hRECKctl |  3.6983  |  0.1210  
%           sTRSVctl |  4.4253  |  0.1210  
%               FA-b |  2.0271  |  2.9204  
%             L2b9a1 |  2.1709  |  2.7494  
%           Theo1041 |  0.9684  |  3.5404  
        grbzVYB=log10([0.3068 4.4253 0.4380 0.7446 2.2635]);
        mdl=fitlm(subrun.ctrl1,grbzVYB);
        fitcoeffs=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2)];
        subrun.VYBmus1=fitcoeffs(1)+fitcoeffs(2)*subrun.mu1;
        subrun.ctrlVYBmus1=fitcoeffs(1)+fitcoeffs(2)*subrun.ctrl1;
        subrun.valVYBmus1=fitcoeffs(1)+fitcoeffs(2)*subrun.val1;

        mdl=fitlm(subrun.adjustedctrl2,grbzVYB);
        fitcoeffs=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2)];
        subrun.VYBmus2=fitcoeffs(1)+fitcoeffs(2)*subrun.adjustedmu2;
        subrun.ctrlVYBmus2=fitcoeffs(1)+fitcoeffs(2)*subrun.adjustedctrl2;
        subrun.valVYBmus2=fitcoeffs(1)+fitcoeffs(2)*subrun.adjustedval2;           
        
    end
    

    function findVYBmus(data,modelparams)
        
        for i=1:length(data.subrun)
            data.subrun(i).VYBmus=modelparams(i,1)*data.subrun(i).adjustedmus+modelparams(i,2);
            for j=1:length(data.motif)
                data.motif(j).subrun(i).VYBmus=data.subrun(i).VYBmus(data.motif(j).hasmotif);
            end
        end
        
    end
    
    function compareMu(data,subrunsID,alpha,minfoldchange)
        for i=1:length(subrunsID(1,:))
%             n1=sum(data.subrun(subrunsID(1,i)).origbincounts,2);
%             n2=sum(data.subrun(subrunsID(2,i)).origbincounts,2);
%             
%             z=(mu1-mu2)./sqrt(mu1.*(1-mu1)./n1+mu2.*(1-mu2)./n2);
            Hstat=zeros(1,length(data.seqs));
            pval=ones(1,length(data.seqs))*nan;
            
            morethan5=find((cellfun('length',data.subrun(subrunsID(1,i)).ttestX)>10) &(cellfun('length',data.subrun(subrunsID(2,i)).ttestX)>10));
            for j=1:length(morethan5)
                [H,P]=ttest2(data.subrun(subrunsID(1,i)).ttestX{morethan5(j)},data.subrun(subrunsID(2,i)).ttestX{morethan5(j)},'alpha',alpha,'tail','both','vartype','unequal');
                Hstat(morethan5(j))=H;
                pval(morethan5(j))=P;
            end

            mu1=data.subrun(subrunsID(1,i)).VYBmus;
            mu2=data.subrun(subrunsID(2,i)).VYBmus;
                        
                        
            isdiff=(abs(mu2-mu1)>log10(minfoldchange));
            
            data.switches(i).name=sprintf('%s vs %s',data.subrun(subrunsID(1,i)).name,data.subrun(subrunsID(2,i)).name);
            data.switches(i).VYBmus1=mu1;
            data.switches(i).VYBmus2=mu2;
            data.switches(i).count1=sum(data.subrun(subrunsID(1,i)).origbincounts,2);
            data.switches(i).count2=sum(data.subrun(subrunsID(2,i)).origbincounts,2);
            
            data.switches(i).isdiff=isdiff;
            data.switches(i).ttestH=Hstat;
            data.switches(i).ttestP=pval;
            data.switches(i).seqs=data.seqs(isdiff);
            
            
            
            
        end
        
        for i=1:length(data.motif)
            for j=1:length(data.switches)
                data.motif(i).switches(j).name=data.switches(j).name;
                data.motif(i).switches(j).VYBmus1=data.switches(j).VYBmus1(data.motif(i).hasmotif);
                data.motif(i).switches(j).VYBmus2=data.switches(j).VYBmus2(data.motif(i).hasmotif);
                data.motif(i).switches(j).isdiff=data.switches(j).isdiff(data.motif(i).hasmotif);
                data.motif(i).switches(j).ttestH=data.switches(j).ttestH(data.motif(i).hasmotif);
                data.motif(i).switches(j).ttestP=data.switches(j).ttestP(data.motif(i).hasmotif);
                data.motif(i).switches(j).count1=data.switches(j).count1(data.motif(i).hasmotif);
                data.motif(i).switches(j).count2=data.switches(j).count2(data.motif(i).hasmotif);
                
                data.motif(i).switches(j).seqs=intersect(data.switches(j).seqs,data.motif(i).seqs,'stable');
            end
        end
        
    end
    
    function extractLoop(data,motifID,regexpmotif,regexpmotifname)
        % regexpmotif assume 2nd and 4th tokens are the loop positions
        for k=1:length(motifID)
            for j=1:length(regexpmotif)
                data.motif(motifID(k)).regexpmotif(j).name=regexpmotifname{j};
                o=regexp(data.motif(motifID(k)).seqs,regexpmotif(j),'tokens');
                hasmotif=find(~cellfun('isempty',o));
                data.motif(motifID(k)).regexpmotif(j).hasmotif=hasmotif;
                data.motif(motifID(k)).regexpmotif(j).loop1={};
                data.motif(motifID(k)).regexpmotif(j).loop2={};
                for i=1:length(hasmotif)
                    data.motif(motifID(k)).regexpmotif(j).loop1{hasmotif(i)}=o{hasmotif(i)}{1}{2};
                    data.motif(motifID(k)).regexpmotif(j).loop2{hasmotif(i)}=o{hasmotif(i)}{1}{4};
                end
                data.motif(motifID(k)).regexpmotif(j).loop1len=cellfun('length',data.motif(motifID(k)).regexpmotif(j).loop1);
                data.motif(motifID(k)).regexpmotif(j).loop2len=cellfun('length',data.motif(motifID(k)).regexpmotif(j).loop2);

            end
        end
        
    end
    
    function seqID=extractSeqbyMu(data,motifID,switchID,prctilemu,numseq)
        if nargin<4
            prctilemu=10;
            numseq=40;
        end
        seqID={};
        for i=1:length(motifID)
            for j=1:length(switchID)
                structOI=data.motif(motifID(i)).switches(switchID(j));
                hasmin=structOI.count2>numseq;
                hasprctilemu=structOI.VYBmus2<prctile(structOI.VYBmus2(hasmin),prctilemu);
                seqID{end+1}=hasmin' & hasprctilemu;
            end
        end
        
    end
    
    function [seqs,c]=findLoopSeq(data,motifID,switchID,regexpmotifID,prctilemu,numseqs,loop,loopsize,showseqlogo);
        seqID=data.extractSeqbyMu(motifID,switchID,prctilemu,numseqs);
        seq10prcmu=find(seqID{1} & seqID{2});
        if loop==1
            hasloopsize=find(data.motif(motifID).regexpmotif(regexpmotifID).loop1len==loopsize);
        elseif loop==2
            hasloopsize=find(data.motif(motifID).regexpmotif(regexpmotifID).loop2len==loopsize);
        end
        [c,ia,ib]=intersect(seq10prcmu,hasloopsize);
        if loop==1
            desiredloop=data.motif(motifID).regexpmotif(regexpmotifID).loop1(c)';
        else
            desiredloop=data.motif(motifID).regexpmotif(regexpmotifID).loop2(c)';
        end

        dloop={};
        for i=1:length(desiredloop) % convert T to U
            s=desiredloop{i};
            s(s=='T')='U';
            dloop{end+1}=s;
            fprintf('%s\n',s)
        end
        % 
        seqs=dloop;
        if showseqlogo
            seqlogo(dloop,'SSCORRECTION',true);
        end
    end
    
    function seqnumeric=convertATCGto1234(data,seqs)
        RNA={'A','U','C','G'};
        seqnumeric=zeros(length(seqs),length(seqs{1}));
        for i=1:length(seqs)
            for j=1:length(RNA)
                seqnumeric(i,seqs{i}==RNA{j})=j;
            end
        end
    end
    
    function [mI,pJoint,pointMI]=findMutualInformation(data,seqs,name)
        % function calculates mutual information 
        % input is an aligned array of sequences in 1,2,3,4
        % output is a 2D vector of mutual information 
        % MI=sum(Px,y(X,Y)*log2(Px,y(X,Y)/Px(X)Py(Y)))
        % and joint probability Px,y(X,Y)
        
        load('MyColormaps','mycmap')
        DNA={'A','U','C','G'};
        addpath ~/Document/MATLAB/cbrewer/cbrewer/cbrewer/
        mycmap=cbrewer('seq','Reds',28);
        
        l1l2seqs=data.convertATCGto1234(seqs);


        pJoint=zeros(length(DNA)*length(l1l2seqs(1,:)));
        for i=1:length(DNA)
            for j=1:length(DNA)
                for k=1:length(l1l2seqs(1,:))
                    for h=1:length(l1l2seqs(1,:))
                        p=sum(l1l2seqs(:,k)==i&l1l2seqs(:,h)==j)/length(l1l2seqs(:,1));
                        pJoint(4*(k-1)+i,4*(h-1)+j)=p;

                    end
                end
            end
        end


        mI=zeros(length(l1l2seqs(1,:)));
        pointMI=zeros(size(pJoint));
        try
        [ent,pAll]=findEnt(l1l2seqs,name);
        catch
            [ent,pAll]=findEnt(l1l2seqs);
        end
        for k=1:length(l1l2seqs(1,:))
            for h=1:length(l1l2seqs(1,:))
                all16mi=zeros(length(DNA));
                for i=1:length(DNA)
                    for j=1:length(DNA)
                        all16mi(i,j)=pJoint(4*(k-1)+i,4*(h-1)+j)*log2(pJoint(4*(k-1)+i,4*(h-1)+j)/pAll(i,k)/pAll(j,h));
                        pointMI(4*(k-1)+i,4*(h-1)+j)=log2(pJoint(4*(k-1)+i,4*(h-1)+j)/pAll(i,k)/pAll(j,h));
                    end
                end
                mI(k,h)=sum(sum(all16mi));
            end
        end
        
        try
            figtitle=strcat(name,': mutual information');
        catch
            figtitle='mutual information';
        end
        setfig(figtitle);clf
        imagesc(mI)
        c=colorbar;
        ylabel(c,'mutual information')
        set(c,'fontsize',28)
        set(c,'linewidth',2)

        colormap(mycmap)
        set(gca,'fontsize',48)
        set(gca,'linewidth',2)
        xlabel('position')
        ylabel('position')
        if nargin==3
        title(figtitle,'interpreter','none')
        end
    end
    
    function plotcompare(data,subrun1,subrun2,numseqs)
        setfig(strcat(data.subrun(subrun1).name,' vs ',data.subrun(subrun2).name));clf
        count1=sum(data.subrun(subrun1).origbincounts,2)>numseqs;
        count2=sum(data.subrun(subrun2).origbincounts,2)>numseqs;
        hasmin=count1&count2;
        plot(data.subrun(subrun1).VYBmus(hasmin),data.subrun(subrun2).VYBmus(hasmin),'.','MarkerSize',10,'color',[0.5 0.5 0.5])
        
        x=linspace(-1,1.5);
        hold on
        plot(x,x,'k:','linewidth',2)
        
        n=sum(hasmin);
        
        text(-0.8,1,sprintf('n=%d',n));
        
    end
    
    function plothandles=plotcomparemotif(data,motifID,subrun1,subrun2,numseqs)
        plothandles={};
        for i=1:length(motifID)
            
            id=motifID(i);
            try
                setfig(strcat(data.motif(id).name,' compare mu'));
            catch
                setfig(strcat(data.motif(id).name{1},' compare mu'));
            end
            
            clf
            c1=sum(data.motif(id).subrun(subrun1).origbincounts,2);
            c2=sum(data.motif(id).subrun(subrun2).origbincounts,2);
            count1=c1>numseqs;
            count2=c2>numseqs;
            hasmin=count1&count2;
            muminus=10.^data.motif(id).subrun(subrun1).VYBmus(hasmin);
            muplus=10.^data.motif(id).subrun(subrun2).VYBmus(hasmin);
            
            
            
%             plot(muminus,muplus,'.','MarkerSize',10,'color',[0.5 0.5 0.5])
            h=scatter(muminus,muplus,'filled');
            set(h,'MarkerFaceColor',[0.2 0.2 0.2])
            s=(c1(hasmin).*c2(hasmin));
            snorm=log10(s)/log10(max(s));
            hold on
            for j=1:length(s)
                set(h,'MarkerFaceAlpha',snorm(j))
                
            end
            plothandles{end+1}=h;

            x=logspace(-1,1.5);
            hold on
            plot(x,x,'k:','linewidth',2)

            n=sum(hasmin);

            text(10.^-0.8,10.^1,sprintf('n=%d',n));
            
            set(gca,'linewidth',2)
            set(gca,'fontsize',22)
            set(gca,'YScale','log')
            set(gca,'XScale','log')
            axis([1e-1 10^1.5 1e-1 10^1.5])
            
            box ON
            
            xlabel('\mu (-ligand)')
            ylabel('\mu (+ligand)')
        end
    end
    
    function plotswitches(data,switchID,numseqs)
        if nargin<2
            switchID=1:length(data.switches);
            numseqs=40;
        end
        
        for i=1:length(switchID)
        setfig(strcat(data.switches(switchID(i)).name,' switches'));clf
        
        hasmin=(data.switches(switchID(i)).count1>numseqs) & (data.switches(switchID(i)).count2>numseqs);
        fprintf('%s: n=%d\n',data.switches(switchID(i)).name,sum(hasmin))
        h1=scatter(10.^data.switches(switchID(i)).VYBmus2(hasmin),10.^data.switches(switchID(i)).VYBmus1(hasmin),'filled');
        set(h1,'MarkerFaceAlpha',0.4)
        box on
        x=logspace(min(data.switches(switchID(i)).VYBmus2(hasmin)),max(data.switches(switchID(i)).VYBmus2(hasmin)));
        hold on
        plot(x,x,'k:','linewidth',2)
        
        mdl=fitlm(data.switches(switchID(i)).VYBmus2(hasmin),data.switches(switchID(i)).VYBmus1(hasmin),'weights',sqrt(data.switches(switchID(i)).count1(hasmin).*data.switches(switchID(i)).count2(hasmin)));
        rsq=mdl.Rsquared.Adjusted;
        
        text(10,1e-1,sprintf('R^2 = %0.2f',rsq),'fontsize',20)
        
        xlabel('\mu (-ligand)')
        ylabel('\mu (+ligand)')
        set(gca,'linewidth',2)
        set(gca,'fontsize',22)
        set(gca,'YScale','log')
        set(gca,'XScale','log')
        
        setfig(strcat(data.switches(switchID(i)).name,' volcano'));clf
        h2=scatter(data.switches(switchID(i)).VYBmus1(hasmin)-data.switches(switchID(i)).VYBmus2(hasmin),-log10(data.switches(switchID(i)).ttestP(hasmin)),'filled');
        set(h2,'MarkerFaceAlpha',0.5)
%         set(gca,'YScale','log')
        xlim([-3,3])
        xlabelnames={};
        for i=-3:1:3
            if i<0
                xlabelnames{end+1}=sprintf('10^-^%d',abs(i));
            else
                xlabelnames{end+1}=sprintf('10^%d',i);
            end
        end
        set(gca,'XTickLabel',xlabelnames)
        end
        

    end
    
    function volcanoplot(~,subrun)
        subrun.fold=(subrun.VYBmus2-subrun.VYBmus1);
        
        scatter(subrun.fold,-log10(subrun.ttestP),'filled','markerfacealpha',0.7,'markerfacecolor',[0.4 0.4 0.4])
        
        hold on
        scatter(subrun.fold(find(subrun.ttestH)),-log10(subrun.ttestP(find(subrun.ttestH))),'filled')
        
        
        [c,ia,ib]=intersect(subrun.seqs,subrun.valseqs);
        
%         plot(subrun.fold(ia),-log10(subrun.ttestP(ia)),'o','linewidth',2,'markersize',14)
        
        xlim([-1,1])
        xlabel('fold change in \mu')
        ylabel('-log10(p-value)')
        set(gca,'fontsize',22)
        set(gca,'linewidth',2)
        
        L=get(gca,'Xticklabel');
        Lnum=[];
        L10powerstr={};
        for i=1:length(L)
            L10powerstr{i}=num2str(10^str2num(L{i}));
        end

        set(gca,'Xticklabel',L10powerstr)
        
    end
    
    function plotswitchesmotif(data,motifID,switchID,numseqs)
        
        for i=1:length(motifID)
            for j=1:length(switchID)
%                 fprintf('i=%d,j=%d; motifID=%d,switchID=%d; \n',i,j,motifID(i),switchID(j))
                try
                    title=sprintf('%s:%s switches',data.motif(motifID(i)).name,data.motif(motifID(i)).switches(switchID(j)).name);
                    setfig(title);clf
                catch
                    title=sprintf('%s:%s switches',data.motif(motifID(i)).name{1},data.motif(motifID(i)).switches(switchID(j)).name);
                    setfig(title);clf
                end
                structOI=data.motif(motifID(i)).switches(switchID(j));
                
                idm=((structOI.VYBmus1-structOI.VYBmus2)>0.3)&((structOI.VYBmus1-structOI.VYBmus2)<log10(50));
                idp=-log10(structOI.ttestP)>12;
                id=find(idm&idp);
                data.motif(motifID(i)).switches(switchID(j)).isswitch=id;
                hasmin=(structOI.count1>numseqs) & (structOI.count2>numseqs);
                hasp=(~isnan(structOI.ttestP))&hasmin';
                fprintf('n=%d\n',sum(hasp))
                h1=scatter(10.^structOI.VYBmus2(hasp),10.^structOI.VYBmus1(hasp),'filled');
                set(h1,'MarkerFaceAlpha',0.5)
                set(h1,'markerfacecolor',[0.4 0.4 0.4])
                hold on
                h12=scatter(10.^structOI.VYBmus2(id),10.^structOI.VYBmus1(id),'filled');
                set(h12,'MarkerFaceColor',[1 1 1])
                h13=scatter(10.^structOI.VYBmus2(id),10.^structOI.VYBmus1(id),'filled');
                set(h13,'MarkerFaceAlpha',0.75)
                set(h13,'MarkerFaceColor',[0.9 0.1 0.1])
                
                x=linspace(min(10.^structOI.VYBmus1(hasp)),max(10.^structOI.VYBmus2(hasp)));
                plot(x,x,'k:','linewidth',2)
                
                xlabel('\mu (-ligand)')
                ylabel('\mu (+ligand)')
                set(gca,'linewidth',2)
                set(gca,'fontsize',22)
                set(gca,'YScale','log')
                set(gca,'XScale','log')
                axis([10^-1,10^1.5,10^-1,10^1.5])
                box ON
                
                try
                    setfig(sprintf('%s:%s volcano',data.motif(motifID(i)).name,data.motif(motifID(i)).switches(switchID(j)).name));clf
                catch
                    setfig(sprintf('%s:%s volcano',data.motif(motifID(i)).name{1},data.motif(motifID(i)).switches(switchID(j)).name));clf
                end
                
                h2=scatter((structOI.VYBmus1-structOI.VYBmus2),-log10(structOI.ttestP),'filled');
                set(h2,'markerfacealpha',0.5)
                set(h2,'markerfacecolor',[0.4 0.4 0.4])
                
                hold on

                h3=scatter((structOI.VYBmus1(id)-structOI.VYBmus2(id)),-log10(structOI.ttestP(id)),'filled');
                set(h3,'MarkerFaceColor',[1 1 1])
                h4=scatter((structOI.VYBmus1(id)-structOI.VYBmus2(id)),-log10(structOI.ttestP(id)),'filled');
                set(h4,'MarkerFaceAlpha',0.75)
                set(h4,'MarkerFaceColor',[0.9 0.1 0.1])
%                 set(gca,'YScale','log')
                set(gca,'linewidth',2)
                set(gca,'fontsize',22)
                xlim([-3,3])
                xlabelnames={};
                for k=-3:1:3
                    if k<0
                        xlabelnames{end+1}=sprintf('10^-^%d',abs(k));
                    else
                        xlabelnames{end+1}=sprintf('10^%d',k);
                    end
                end
                set(gca,'XTickLabel',xlabelnames)
                xlabel('fold change in \mu')
                ylabel('-log10(p-value)')
                box ON
                
                try
                    setfig(sprintf('%s:%s enrichment vs count',data.motif(motifID(i)).name,data.motif(motifID(i)).switches(switchID(j)).name));clf
                catch
                    setfig(sprintf('%s:%s enrichment vs count',data.motif(motifID(i)).name{1},data.motif(motifID(i)).switches(switchID(j)).name));clf
                end                
                
                notseqerror=(structOI.count2>600)' & (structOI.VYBmus2<.67);
                outlier=find(abs(structOI.VYBmus2-0.5803)<0.0001);
                notseqerror(outlier)=0;
                notseqerror=(structOI.count2>10);
                
                
                h5=scatter(log10(structOI.count2(notseqerror)),structOI.VYBmus2(notseqerror),50,'filled');
                set(h5,'markerfacealpha',0.5)
                set(h5,'markerfacecolor',[0.4 0.4 0.4])
                
                mdl=fitlm(log10(structOI.count2(notseqerror)),structOI.VYBmus2(notseqerror));
                rsq=mdl.Rsquared.Adjusted;

                fitcoeffs=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2)];
                x=linspace(min(log10(structOI.count2(notseqerror))),max(log10(structOI.count2(notseqerror))));
                hold on
                plot(x,(fitcoeffs(1)+fitcoeffs(2)*(x)),'k:','linewidth',2)
                
                text(3.6,.5,sprintf('R^2=%0.2f',rsq),'fontsize',16)
                text(3.6,.45,sprintf('y=%0.2fx+%0.2f',fitcoeffs(2),fitcoeffs(1)),'fontsize',16)

                set(gca,'linewidth',2)
                set(gca,'fontsize',22)
%                 set(gca,'XScale','log')
                xlabel('log10(NGS read count)')
                ylabel('log10(\mu)')
                box on
            end
        end
    end
    
    
    function plotenrichmentVScountmotif(data,motifID,switchID,numseqs)
        
        
        
        
    end

  end
end










