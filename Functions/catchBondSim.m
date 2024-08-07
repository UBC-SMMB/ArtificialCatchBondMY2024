%Catch bond simulator:
%First, model force dependent lifetime of one unzipping dna strand
%then, model another one
%If they intersect, you got a potential catch bond candidate!
%If they don't intersect, one will always rupture before the other and you
%have boring slip bond behaviour.
%The DNASeqFreeEnergy function is sequence dependent, not just CG content
%dependent, which is really nice!


%Things that affect force dependent lifetime (Murad & Li 2019):
%Slope of log(tau(F)) is affected by barrier width, which is mostly
%controlled by length of sequence
%Vertical shift is affected by the ratio between barrier width and height;
%barrier height is affected by CG content/actual sequence.

%so let's test:
% percent CG content
% length
%Things that may be important but we didn't test:
% number of CG
% order of CG (which I don't think yousif did, but might be important here)

%input: 
% nSeq: number of sequences to generate (1 x 1 double)
% L: length of sequence (1 x 1 double)
% nCG: number of CG (1 x 1 double, must be <= L)
% salt: salt molarity (struct containing .Na, .Mg, .K [M])
% plotresults: plot option (1 x 1 boolean)
function [seqProps, seq, tau] = catchBondSim(nSeq, L, nCG, salt, plotresults)
%% Step 1/4: Initial parameters
loop = '';                      %if it's empty there's no loop
T=300;                          %[K]

%% Step 2/4: Sequence and lifetime generation

seq = randomDNAseq(L,nCG,nSeq);

%% Step 3/4: Lifetime generation

F = 0:2e-12:20e-12;                       %[N]

tau = zeros(length(seq),length(F))';      %[s]
dGf1 = zeros(length(seq),length(F))';      %[kbT]


    for i = 1:length(F)
        for j = 1:length(seq)
            [tau(i,j),Fmax,~,~,dGf1(i,j)] = tauDNASeqUnzip(char(seq(j)),loop,salt,F(i),T);
        end
    end

%% Step 4/4: Getting parameters for output
 
 FpN = F*1e12;                   %convert F to pN from N

 for i = 1:length(FpN)
     seqProps.median(i) = median(tau(i,:));
     seqProps.mean(i) = mean(tau(i,:));
     seqProps.std(i) = std(tau(i,:));
     seqProps.error(i) = seqProps.std(i)/sqrt(nSeq);
 end
 seqProps.F = FpN;
 seqProps.Fmax = Fmax;
 seqProps.Fmin = 1e-2;

%%
if plotresults == true

    for i = 1:length(FpN)
        figure;
        histogram(tau(i,:));
        xline(seqProps.mean(i),'b')
        xline(seqProps.median(i),'r');
        xline(seqProps.mean(i) + seqProps.std(i),'m');
        xline(seqProps.mean(i) - seqProps.std(i),'m');
    end

end