function seq1 = readFeaturesFromDAtA(P, QRS, T, flat, nbeats, Nsig)
%READFEATURESFROMDATA put together the features of the different states
%
%   input:  Pwave, QRScomp, Twave, idle
%           nbBeats     number of beats
%           nbSig       number of signals
%   output: seq         cell with the 

seq1 = [];

if length(Nsig) == 1
    sSig = 1;
    nSig = Nsig;
else
    sSig = Nsig(1);
    nSig = Nsig(2);
end

if length(nbeats) == 1
    sBeat = 1;
    eBeat= nbeats;
else
    sBeat = nbeats(1);
    eBeat = nbeats(2);
end

a=1;
for ll = sSig : nSig
    for ii = sBeat : eBeat
                 
            type = [T{1,ll}.type{1,ii};...
                  P{1,ll}.type{1,ii};...
                  QRS{1,ll}.type{1,ii};...
                  flat{1,ll}.type{1,ii}];
            dS = [T{1,ll}.dS{1,ii};...
                  P{1,ll}.dS{1,ii};...
                  QRS{1,ll}.dS{1,ii};...
                  flat{1,ll}.dS{1,ii}];
            deltaDS = [T{1,ll}.deltaDS{1,ii};...
                  P{1,ll}.deltaDS{1,ii};...
                  QRS{1,ll}.deltaDS{1,ii};...
                  flat{1,ll}.deltaDS{1,ii}];
            summe = [T{1,ll}.summe{1,ii};...
                  P{1,ll}.summe{1,ii};...
                  QRS{1,ll}.summe{1,ii};...
                  flat{1,ll}.summe{1,ii}];
            count_dt = [T{1,ll}.count_dt{1,ii};...
                  P{1,ll}.count_dt{1,ii};...
                  QRS{1,ll}.count_dt{1,ii};...
                  flat{1,ll}.count_dt{1,ii}];
            du_dt = [T{1,ll}.du_dt{1,ii};...
                  P{1,ll}.du_dt{1,ii};...
                  QRS{1,ll}.du_dt{1,ii};...
                  flat{1,ll}.du_dt{1,ii}];

            
              
            seq{a} = ([type, dS, deltaDS, summe, count_dt, du_dt])';
            a = a + 1;
    end

    seq1 = [seq1, seq];
end
end