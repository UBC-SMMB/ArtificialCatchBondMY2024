%% Load and process data RUN THIS FIRST
c1 = '#80CCEA';
c2 = '#C44596';
c3 = '#911136';
c1rgb = [hex2dec(c1(2:3))/256 hex2dec(c1(4:5))/256 hex2dec(c1(6:7))/256];
c2rgb = [hex2dec(c2(2:3))/256 hex2dec(c2(4:5))/256 hex2dec(c2(6:7))/256];
c3rgb = [hex2dec(c3(2:3))/256 hex2dec(c3(4:5))/256 hex2dec(c3(6:7))/256];


slowunwrapped = [slowhoard{:}];
slowrupF = [slowunwrapped.rupF];
slowLR = [slowunwrapped.loadingRate];
slowoff = [slowunwrapped.offsets];

medunwrapped = [medhoard{:}];
medrupF = [medunwrapped.rupF];
medLR = [medunwrapped.loadingRate];
medoff = [medunwrapped.offsets];

fastunwrapped = [fasthoard{:}];
fastrupF = [fastunwrapped.rupF];
fastLR = [fastunwrapped.loadingRate];
fastoff = [fastunwrapped.offsets];

mslowLR = median(slowLR);
mmedLR = median(medLR);
mfastLR = median(fastLR);

v = [51.67,206.7, 2067];
v2 = -2:0.5:6;

edges = 0:2:55;

T = 300;
F = 0:0.1:130;
l = [forceTruncatedIntersections.hookL];
xddagger = 0.238*2 + (l-2) .* (0.005.*l + 0.085);
b = physconst('Boltzmann');         %J/K

tjaw = 10.^([forceTruncatedIntersections.jawInt] + F'.* [forceTruncatedIntersections.jawSlope]);
thookunzip = 10.^([forceTruncatedIntersections.hookInt] + F'.* [forceTruncatedIntersections.hookSlope]);
thookshear = 10.^[forceTruncatedIntersections.hookInt] .* exp((-F'.*xddagger.*10.^-21)/(b*T));
pjaw = thookunzip./(tjaw + thookunzip);
phook = tjaw./(tjaw + thookunzip);
meant = pjaw .* thookshear + phook .* thookunzip;

dmeant = diff(meant,1,1);
sd = size(dmeant);
dmeant1 = [dmeant; zeros(1,sd(2))];
dmeant2 = [zeros(1,sd(2)); dmeant];
inflection = (dmeant1 > 0 & dmeant2<0) | (dmeant2 > 0 & dmeant1<0);
inflectiontau = meant.*inflection;
inflectionF = repmat(F',1,sd(2)).*inflection;

for i = 1:sd(2)
    its = (nonzeros(inflectiontau(:,i)));
    log10its = log10(its);
    ifs = nonzeros(inflectionF(:,i));
    dt = diff(its);
    df = diff(ifs);
    dlog10its = diff(log10its);
    if ~length(dt)
        dt = nan;
        df = nan;
        its = [nan nan];
        ifs = [nan nan];
        dlog10its = nan;
    end
    deltalog10tau(i) = dlog10its;
    deltatau(i) = dt;
    deltaF(i) = df;
    inflectiontautrunc(:,i) = its;
    inflectionFtrunc(:,i) = ifs;
end
catchtf = ~isnan(deltaF) & min(inflectiontautrunc) > 0.01 & min(inflectionFtrunc) < 10 & max(inflectiontautrunc) < 100 ;
catchbonds = forceTruncatedIntersections(catchtf);
catchbonddeltaF = deltaF(catchtf);
catchbonddeltalog10tau = deltalog10tau(catchtf);
catchbondmeant = meant(:,catchtf);

colortauF = [rescale([catchbonddeltalog10tau]','InputMax',max(catchbonddeltalog10tau),'InputMin',min(catchbonddeltalog10tau)), zeros(length(catchbonddeltalog10tau),1), rescale([catchbonddeltaF]','InputMax',max(catchbonddeltaF),'InputMin',min(catchbonddeltaF))];

pd = polymerParams('dsDNA','pNnm');
ps = polymerParams('ssDNA','pNnm');
pp = polymerParams('PEG','pNnm');
pd.S = 1200;

interval = 0.1;
plotF = 0.2:interval:35;

pegLength = pp.Lc*9;

dsDNALength1 = pd.Lc*2914;
ssDNALength1 = ps.Lc*4;

dsDNALength2 = pd.Lc*(2914 + 9);
ssDNALength2 = ps.Lc*42;

o = 265;

curve1 = XWLCContour(plotF,pd.Lp,pd.S,pd.T)*dsDNALength1 + XWLCContour(plotF,ps.Lp,ps.S,ps.T)*ssDNALength1 + XWLCContour(plotF,pp.Lp,pp.S,pp.T)*pegLength;
curve2 = XWLCContour(plotF,pd.Lp,pd.S,pd.T)*dsDNALength2 + XWLCContour(plotF,ps.Lp,ps.S,ps.T)*ssDNALength2+ XWLCContour(plotF,pp.Lp,pp.S,pp.T)*pegLength;
lc1 = dsDNALength1 + ssDNALength1 + pegLength;
lc2 = dsDNALength2 + ssDNALength2 + pegLength;

%% Figure 2d

figure;
idx = 9952;
hold on
ylabel('Lifetime (s)')
xlabel('Force (pN)')
yscale('log')
plot(F,tjaw(:,idx))
plot(F,thookunzip(:,idx))
plot(F,thookshear(:,idx))
plot(F,meant(:,idx))
xlim([0 20]);

%% Figure 2e-f

figure;
hold on;
plot(curve1,plotF,"Color",c1,"LineWidth",3);
plot(curve2,plotF,"Color",c3,"LineWidth",3);
xlabel('Extension (nm)')
ylabel('Force (pN)')

%% Figure 3a

figure;
tiledlayout('horizontal');
nexttile;
hold on;
for i = 1:length(FEC)
    line(FEC(i).x(1:5:end) + fecoffsets(i)/2 + o,repmat(i,size(FEC(i).F(1:5:end))),FEC(i).F(1:5:end),'Color',fecswitched(i,:))
end
view(-55,18);
xlabel('Extension (nm)')
ylabel('Trace #')
zlabel('Force (pN)')
zlim([-2, 40])

%% Figure 3b

figure;
tiledlayout('horizontal');
nexttile;
hold on;
for i = 1:length(FEC)
    line(FEC(i).x(1:5:end) + fecoffsets(i)/2 + o,FEC(i).F(1:5:end),'Color',fecswitched(i,:))
end
xlabel('Extension (nm)')
ylabel('Force (pN)')

%% Figure 3b (inset)

figure
histogram(medfit(4,medjawopen)-medfit(3,medjawopen),15:0.3:25,'FaceColor',c2,'FaceAlpha',1)
ylim([0 150])
xline(lc2-lc1,'k--','LineWidth',2)
xlabel('\Delta x (nm)')

%% Figure 3c

figure;
t1 = tiledlayout(3,2,"Padding","tight",'TileSpacing','compact');
nexttile(t1,1);
histogram(slowrupF,edges,'FaceColor',c1,'FaceAlpha',1,'EdgeColor','k');
title(['Loading Rate = ' num2str(v(1)) ' nm/s'])
legend(['n = ' num2str(length(slowrupF))],'Location','north')

nexttile(t1,3);
histogram(medrupF,edges,'FaceColor',c2,'FaceAlpha',1,'EdgeColor','k');
title(['Loading Rate = ' num2str(v(2)) ' nm/s'])
legend(['n = ' num2str(length(medrupF))],'Location','north')

nexttile(t1,5);
histogram(fastrupF,edges,'FaceColor',c3,'FaceAlpha',1,'EdgeColor','k');
title(['Loading Rate = ' num2str(v(3)) ' nm/s'])
legend(['n = ' num2str(length(fastrupF))],'Location','north')

a = {};

for i = 1:length(v)
    switched_sum = sum(nnz(rswitched_lr(i,:)));
    a{end+1} = nexttile(t1,i*2);
    histogram(rrupture_forces(i,:),edges);
    title(['Loading Rate, ' num2str(v(i)) ' nm/s '])
end

set(a{1}.Children(1), 'FaceColor',c1,'FaceAlpha',1,'EdgeColor','k');
set(a{2}.Children(1), 'FaceColor',c2,'FaceAlpha',1,'EdgeColor','k');
set(a{3}.Children(1), 'FaceColor',c3,'FaceAlpha',1,'EdgeColor','k');
xlabel(t1,'Rupture Force (pN)');

set(a{1}.Children(1), 'FaceColor',c1,'FaceAlpha',1,'EdgeColor','k');
set(a{2}.Children(1), 'FaceColor',c2,'FaceAlpha',1,'EdgeColor','k');
set(a{3}.Children(1), 'FaceColor',c3,'FaceAlpha',1,'EdgeColor','k');
xlabel(t1,'Rupture Force (pN)');

%% Figure 3d

figure;
hold on;
xscale(gca,"log");

scatter(slowLR(~slowjawopen), slowrupF(~slowjawopen),[],c1rgb,'.')
scatter(slowLR(slowjawopen), slowrupF(slowjawopen),[],c3rgb,'.')
scatter(medLR(~medjawopen), medrupF(~medjawopen),[],c1rgb,'.')
scatter(medLR(medjawopen), medrupF(medjawopen),[],c3rgb,'.')
scatter(fastLR(~fastjawopen), fastrupF(~fastjawopen),[],c1rgb,'.')
scatter(fastLR(fastjawopen), fastrupF(fastjawopen),[],c3rgb,'.')
xlabel('Loading Rate (pN/s)')
ylabel('Rupture Force (pN)')

scatter(rrupture_lr(1,:),rrupture_forces(1,:),'ko','filled','MarkerFaceAlpha',0.1,'SizeData',5)
scatter(rrupture_lr(2,:),rrupture_forces(2,:),'ko','filled','MarkerFaceAlpha',0.1,'SizeData',5)
scatter(rrupture_lr(3,:),rrupture_forces(3,:),'ko','filled','MarkerFaceAlpha',0.1,'SizeData',5)

%% Figure 3e

figure; xscale("log")
hold on
errorbar(10.^v2, sp, sqrt(sp.*(1-sp)./1000),'m.')
plot(10.^v2, sp,'k')


errorbar(v, p0809, sqrt(p0809.*(1-p0809)./[325 609 463])*1.96, 'o','Color',c1,'MarkerFaceColor','auto')
errorbar(v, p0811, sqrt(p0811.*(1-p0811)./[101 178 130])*1.96, 'o','Color',c2,'MarkerFaceColor','auto')
errorbar(v, p1009, sqrt(p1009.*(1-p1009)./[253 490 187])*1.96, 'o','Color',c3,'MarkerFaceColor','auto')
xlabel("Pulling Rate (nm/s)")
ylabel("P_F of Jaw Opening")
plot(10.^v2, sp,'b-')

%% Figure 4-part 2

test = [[catchbonds.jawCG]' [catchbonds.jawL]'];
[~,~,ic] = unique(test, 'rows', 'stable');
h = accumarray(ic, 1);
maph = h(ic);
colortest = [maph];

figure;
nexttile;
	iscJaw = nan(max([catchbonds.jawCG]) - min([catchbonds.jawCG])+ 1,max([catchbonds.jawL])-min([catchbonds.jawL])+ 1);
	for i = 1:length(colortest)
    	x = catchbonds(i).jawCG + 1;
    	y = catchbonds(i).jawL + 1;
   	iscJaw(x,y) = colortest(i);
	end

	catchbonds_fa = struct2fieldarray(catchbonds);
	jawCG_min = min(catchbonds_fa.jawCG);
	jawCG_max = max(catchbonds_fa.jawCG);
	jawL_min = min(catchbonds_fa.jawL);
	jawL_max = max(catchbonds_fa.jawL);

	imagesc([0 jawL_max],[jawCG_max 0],flipud(iscJaw));
	axis equal;
	axis tight;
	ylim([-1 12]);
	xlim([10 31]);

	set(gca,'YDir','normal')
	colormap(gca, 'jet');
	colormap([1 1 1; colormap]);
	ylabel('Jaw CG')
	xlabel('Jaw Length')
	colorbar

nexttile;
	test = [[catchbonds.hookCG]' [catchbonds.hookL]'];
	[~,~,ic] = unique(test, 'rows', 'stable');
	h = accumarray(ic, 1);
	maph = h(ic);
	colortest = [maph];
	
	iscHook = nan(max([catchbonds.hookCG]) - min([catchbonds.hookCG])+ 1,max([catchbonds.hookL])-min([catchbonds.hookL])+ 1);
	for i = 1:length(colortest)
    	x = catchbonds(i).hookCG + 1;
    	y = catchbonds(i).hookL + 1;
   	iscHook(x,y) = colortest(i);
	end

	catchbonds_fa = struct2fieldarray(catchbonds);
	hookCG_min = min(catchbonds_fa.hookCG);
	hookCG_max = max(catchbonds_fa.hookCG);
	hookL_min = min(catchbonds_fa.hookL);
	hookL_max = max(catchbonds_fa.hookL);
	
	imagesc([0 hookL_max],[hookCG_max 0],flipud(iscHook)); 
	axis equal;
	axis tight;
	ylim([-1 10]);
	xlim([6 22]);

	set(gca,'YDir','normal')
	colormap(gca, 'jet');
	colormap([1 1 1; colormap]);
	
	ylabel('Hook CG')
	xlabel('Hook Length')
	colorbar

fixedHook = [catchbonds.hookL] == 8 & [catchbonds.hookCG] == 8;
fixedJaw = [catchbonds.jawL] == 27 & round([catchbonds.jawCG],1) == 0;

nexttile;
	hold on;
	title('Jaw Length: 27; CGs: 0')
	scatter([catchbonds(fixedJaw).hookL],[catchbonds(fixedJaw).hookCG],[],colortauF(fixedJaw,:),'.');
	xlabel('Length (bp)')
	ylabel('CG Proportion')
	axis equal
	xlim([5 20])
	ylim([-1 10])
	

nexttile;
	hold on;
	title('Hook Length: 8; CGs: 8')
	scatter([catchbonds(fixedHook).jawL],[catchbonds(fixedHook).jawCG],[],colortauF(fixedHook,:),'.');
	xlabel('Length (bp)')
	ylabel('CG #')
	axis equal
	xlim([10 30])
	ylim([-1 10])

%% Figure 4-part 1

figure

ax1 = nexttile;
	inflection_xx = inflectionFtrunc(:,catchtf);
	inflection_yy = inflectiontautrunc(:,catchtf);
	slope = catchbonddeltalog10tau./catchbonddeltaF;
	[~,sort_order] = sort(slope);
	inflection_xx = inflection_xx(:,sort_order);
	inflection_yy = inflection_yy(:,sort_order);
	colortauF_sorted = colortauF(sort_order,:);

	plot_gap = 1;
	line(inflection_xx(:,1:plot_gap:end),inflection_yy(:,1:plot_gap:end));
	colororder(ax1,colortauF_sorted(1:plot_gap:end,:))
	yscale('log')
	ylabel('\tau (s)')
	xlabel('F (pN)')

nexttile 
	hold on;
	scatter(catchbonddeltaF,catchbonddeltalog10tau,[],colortauF,'.')
	ylabel('\Delta log(\tau)')
	xlabel('\Delta F (pN)')

nexttile
	hold on

	[~,idx] = max(catchbonddeltaF);
	plot(F,catchbondmeant(:,idx),"Color",colortauF(idx,:))
	yscale('log')
	xlim([0 50])
	ylim([0.001 1000])
	xlabel('Force (pN)')
	ylabel('\tau (s)')

nexttile
	hold on
	[~,idx] = max(catchbonddeltalog10tau);
	plot(F,catchbondmeant(:,idx),"Color",colortauF(idx,:))
	
	yscale('log')
	xlim([0 50])
	ylim([0.001 1000])
