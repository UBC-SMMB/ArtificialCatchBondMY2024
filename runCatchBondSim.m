%% 1: Generate Sequences and Sequence tauF
%PBS:   Na: 0.153
%       K: 0.0047
%       Mg: 0
%T50M5: Na: 0.050
%       K: 0
%       Mg: 0.005

nSeq = 1000;            %number of sequences per condition
sLength = 7:30;         %conditions: length
salt.Na=0.050;          %salt concentrations [M]
salt.K=0;                       
salt.Mg=0.005;
tic()
for L = 1:length(sLength)
    for CG = 1:sLength(L) + 1
        [num2str(sLength(L)) ':' num2str(CG-1)]
        [seqProps(L,CG), seq{L,CG}, tau{L,CG}] = catchBondSim(nSeq, sLength(L), CG - 1,salt,false);
    end
end
toc()
F = seqProps(1,1).F;

for L = 1:length(sLength)
    for CG = 1:sLength(L) + 1
        seqProps(L,CG).CG = CG - 1;
        seqProps(L,CG).L = sLength(L);
        seqProps(L,CG).median = median(tau{L,CG},2);
        seqProps(L,CG).line = fitlm(F, log10(seqProps(L,CG).median));
    end
end


%% 2: find valid intersections

%Produces intersections:

%F              = Fc            (pN)
%tau            = log(tauc)     (log(s)) 
%jawL           = jaw length    (bp) 
%jawCG          = jaw cg        (bp)
%jawSlope       = slope of log(tau) versus F
%jawInt         = intercept of log(tau) versus F
%jawR           = R value 
%hookL          = hook length (bp)
%hookCG         = hook cg (bp)
%hookSlope      = slope of log(tau) versus F
%hookInt        = intercept of log(tau) versus F
%hookR          = R value

intersections = [];
tic()
for Ls = 1:length(sLength)
    for CGs = 1:sLength(Ls) + 1
        slopeJaw = seqProps(Ls,CGs).line.Coefficients{2,1};
        for Lb = 1:length(sLength)
            for CGb = 1:sLength(Lb) + 1
                slopeHook = seqProps(Lb,CGb).line.Coefficients{2,1};
                if slopeJaw < slopeHook % this is good, means jaw will have a longer lifetime than hook before crossing and higher after
                    intJaw = seqProps(Ls,CGs).line.Coefficients{1,1};
                    intHook = seqProps(Lb,CGb).line.Coefficients{1,1};
                    RJaw = seqProps(Ls,CGs).line.Rsquared.Ordinary;
                    RHook = seqProps(Lb,CGb).line.Rsquared.Ordinary;
                    [x,y] = line_intersection([slopeJaw, intJaw],[slopeHook,intHook]);
                    intersections(end + 1).F = x;
                    intersections(end).tau = y;
                    intersections(end).jawL =  sLength(Ls);
                    intersections(end).jawCG = CGs - 1;
                    intersections(end).jawSlope = slopeJaw;
                    intersections(end).jawInt = intJaw;
                    intersections(end).jawR = RJaw;
                    intersections(end).hookL = sLength(Lb);
                    intersections(end).hookCG = CGb - 1;
                    intersections(end).hookSlope = slopeHook;
                    intersections(end).hookInt = intHook;
                    intersections(end).hookR = RHook;
                    intersections(end).jawRaw = seqProps(Ls,CGs).median;
                    intersections(end).hookRaw = seqProps(Lb,CGb).median;
                end
            end
        end
    end
end
toc()
%% 3. Monte Carlo Simulation to get catch bond rupture times

%% Mode 1: non-linear (due to DNA handle elasticity) force ramp experiment
idx = 11587;                        %row number in intersections
T = 300;                            %temp (K)
Fparam = 51.7;                      %ramp speed (nm/s)
mode = 'ramps2f';                   %mode
experiments = 1000;                 %n
p = polymerParams('dsDNA','pNnm');  %polymer handle parameters
nbp = 2914;                         %number of base pairs in handles

[rupture_times, switched_array,switched_times, rupture_forces,switched_forces,lr,slr] = getCatchBondRuptureTimes(intersections(idx),Fparam,T,mode,experiments,p,nbp);

figure;
tiledlayout("horizontal")
nexttile;
histogram(rupture_forces)
nexttile; 
plot(lr,rupture_forces,'o')
%% Mode 2: linear force ramp experiment
idx = 11587;                        %row number in intersections
T = 300;                            %temp (K)
Fparam = 70;                        %loading rate (pN/s)
mode = 'ramp';                      %mode
experiments = 1000;                 %n

[rupture_times, switched_array,switched_times, rupture_forces,switched_forces] = getCatchBondRuptureTimes(intersections(idx),Fparam,T,mode,experiments);
figure;
histogram(rupture_forces)
%% Mode 3:  for clamp experiment
idx = 11587;                        %row number in intersections
T = 300;                            %temp (K)
Fparam = 0:2:30;                    %Forces (pN)
mode = 'clamp';                     %mode
experiments = 1000;                 %n

[rupture_times, switched_array,switched_times] = getCatchBondRuptureTimes(intersections(idx),Fparam,T,mode,experiments);
figure;
semilogy(Fparam, mean(rupture_times,1))