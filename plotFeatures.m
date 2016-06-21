function plotFeatures(idle, Pwave, QRScomp,Twave, feat1, feat2, Q, beats, marker)

switch feat1
    case 1
        idlePNFeat1 = idle{1,1}.typePN;
        PwaveFeat1 = Pwave{1,1}.type;
        idleSTFeat1 = idle{1,1}.typeST;
        QRSFeat1 = QRScomp{1,1}.type;
        idleTPFeat1 = idle{1,1}.typeTP;
        TwaveFeat1 = Twave{1,1}.type;
        idleFeat1 = idle{1,1}.type;
    case 2
        idlePNFeat1 = idle{1,1}.dS_PN;
        PwaveFeat1 = Pwave{1,1}.dS;
        idleSTFeat1 = idle{1,1}.dS_ST;
        QRSFeat1 = QRScomp{1,1}.dS;
        idleTPFeat1 = idle{1,1}.dS_TP;
        TwaveFeat1 = Twave{1,1}.dS;
        idleFeat1 = idle{1,1}.dS;
    case 3
        idlePNFeat1 = idle{1,1}.deltaDS_PN;
        PwaveFeat1 = Pwave{1,1}.deltaDS;
        idleSTFeat1 = idle{1,1}.deltaDS_ST;
        QRSFeat1 = QRScomp{1,1}.deltaDS;
        idleTPFeat1 = idle{1,1}.deltaDS_TP;
        TwaveFeat1 = Twave{1,1}.deltaDS;
        idleFeat1 = idle{1,1}.deltaDS;
    case 4
        idlePNFeat1 = idle{1,1}.summePN;
        PwaveFeat1 = Pwave{1,1}.summe;
        idleSTFeat1 = idle{1,1}.summeST;
        QRSFeat1 = QRScomp{1,1}.summe;
        idleTPFeat1 = idle{1,1}.summeTP;
        TwaveFeat1 = Twave{1,1}.summe;
        idleFeat1 = idle{1,1}.summe;
    case 5
        idlePNFeat1 = idle{1,1}.count_dtPN;
        PwaveFeat1 = Pwave{1,1}.count_dt;
        idleSTFeat1 = idle{1,1}.count_dtST;
        QRSFeat1 = QRScomp{1,1}.count_dt;
        idleTPFeat1 = idle{1,1}.count_dtTP;
        TwaveFeat1 = Twave{1,1}.count_dt;
        idleFeat1 = idle{1,1}.count_dt;
    case 6
        idlePNFeat1 = idle{1,1}.du_dtPN;
        PwaveFeat1 = Pwave{1,1}.du_dt;
        idleSTFeat1 = idle{1,1}.du_dtST;
        QRSFeat1 = QRScomp{1,1}.du_dt;
        idleTPFeat1 = idle{1,1}.du_dtTP;
        TwaveFeat1 = Twave{1,1}.du_dt;
        idleFeat1 = idle{1,1}.du_dt;
end
        

switch feat2
    case 1
        idlePNFeat2 = idle{1,1}.typePN;
        PwaveFeat2 = Pwave{1,1}.type;
        idleSTFeat2 = idle{1,1}.typeST;
        QRSFeat2 = QRScomp{1,1}.type;
        idleTPFeat2 = idle{1,1}.typeTP;
        TwaveFeat2 = Twave{1,1}.type;
        idleFeat2 = idle{1,1}.type;
    case 2
        idlePNFeat2 = idle{1,1}.dS_PN;
        PwaveFeat2 = Pwave{1,1}.dS;
        idleSTFeat2 = idle{1,1}.dS_ST;
        QRSFeat2 = QRScomp{1,1}.dS;
        idleTPFeat2 = idle{1,1}.dS_TP;
        TwaveFeat2 = Twave{1,1}.dS;
        idleFeat2 = idle{1,1}.dS;
    case 3
        idlePNFeat2 = idle{1,1}.deltaDS_PN;
        PwaveFeat2 = Pwave{1,1}.deltaDS;
        idleSTFeat2 = idle{1,1}.deltaDS_ST;
        QRSFeat2 = QRScomp{1,1}.deltaDS;
        idleTPFeat2 = idle{1,1}.deltaDS_TP;
        TwaveFeat2 = Twave{1,1}.deltaDS;
        idleFeat2 = idle{1,1}.deltaDS;
    case 4
        idlePNFeat2 = idle{1,1}.summePN;
        PwaveFeat2 = Pwave{1,1}.summe;
        idleSTFeat2 = idle{1,1}.summeST;
        QRSFeat2 = QRScomp{1,1}.summe;
        idleTPFeat2 = idle{1,1}.summeTP;
        TwaveFeat2 = Twave{1,1}.summe;
        idleFeat2 = idle{1,1}.summe;
    case 5
        idlePNFeat2 = idle{1,1}.count_dtPN;
        PwaveFeat2 = Pwave{1,1}.count_dt;
        idleSTFeat2 = idle{1,1}.count_dtST;
        QRSFeat2 = QRScomp{1,1}.count_dt;
        idleTPFeat2 = idle{1,1}.count_dtTP;
        TwaveFeat2 = Twave{1,1}.count_dt;
        idleFeat2 = idle{1,1}.count_dt;
    case 6
        idlePNFeat2 = idle{1,1}.du_dtPN;
        PwaveFeat2 = Pwave{1,1}.du_dt;
        idleSTFeat2 = idle{1,1}.du_dtST;
        QRSFeat2 = QRScomp{1,1}.du_dt;
        idleTPFeat2 = idle{1,1}.du_dtTP;
        TwaveFeat2 = Twave{1,1}.du_dt;
        idleFeat2 = idle{1,1}.du_dt;
end


for ll = 1 : beats
    if Q == 6
        plot(idlePNFeat1{1, ll}, idlePNFeat2{1, ll}, marker{1})
        plot(PwaveFeat1{1, ll}, PwaveFeat2{1, ll}, marker{2})
        plot(idleSTFeat1{1, ll}, idleSTFeat2{1, ll}, marker{3})               
        plot(QRSFeat1{1, ll}, QRSFeat2{1, ll}, marker{4})
        plot(idleTPFeat1{1, ll}, idleTPFeat2{1, ll}, marker{5})                
        plot(TwaveFeat1{1, ll}, TwaveFeat2{1, ll}, marker{6})
    elseif Q == 4
        plot(idleFeat1{1, ll}, idleFeat2{1, ll}, marker{1})
        plot(PwaveFeat1{1, ll}, PwaveFeat2{1, ll}, marker{2})
        plot(QRSFeat1{1, ll}, QRSFeat2{1, ll}, marker{3})
        plot(TwaveFeat1{1, ll}, TwaveFeat2{1, ll}, marker{4})
    end
end
end