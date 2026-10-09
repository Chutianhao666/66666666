function result = prachDetection(factor, snr, pdp, prachConfig)
NIFFT = prachConfig.bwpInfo.NFFT;
L_RA = prachConfig.L_RA;
N_CS = prachConfig.N_CS;

overSampling = NIFFT/L_RA;
windowSize = ceil(overSampling*N_CS);

thr = thresholdCalculation(factor,snr);

numDet = 0;
valid = [];
position = [];

figure; subplot(211); hold on; plot(pdp(:,1),'LineWidth',1); plot(thr,'LineWidth',1); subplot(212); hold on; plot(pdp(:,2),'LineWidth',1); plot(thr,'LineWidth',1);


% % used for detection probability
% if size(pdp,2) > 1
%     pdp = pdp(:,1);
% end



for i = 1:windowSize
    if pdp(i) > thr(i)
        numDet = numDet + 1;
        position(numDet) = i;
        valid(numDet) = 1;
    end
end

if numDet == 0
    result.detection = 0;
    result.miss_detection = 1;

elseif numDet == 1

    result.detection = 1;
    result.miss_detection = 0;

else
    for n = 2:numDet
        if position(n) - position(n-1) < overSampling
            if pdp(position(n)) >= pdp(position(n-1))
                valid(n-1) = 0;
            else
                valid(n) = 0;
            end
            numDet = numDet - 1;
        end
    end

    if numDet == 0
  
        result.detection = 0;
        result.miss_detection = 1;

    elseif numDet == 1
     
        result.detection = 1;
        result.miss_detection = 0;
    else
   
        result.detection = 1;
        result.miss_detection = 0;
    end

end

% figure; plot(pdp); hold on; plot(thr);



end

%% 
function thr = thresholdCalculation(factor, snr)

    peak_snr = [0.1374906761928795
        0.11139271280810693
        0.09024263064526394
        0.07244720262852593
        0.05812212341301452
        0.04683069276082391
        0.03758086483009969
        0.03020266908073748
        0.02405773785498541
        0.019427935442184744
        0.01569200494080036
        0.012624200087737512
        0.010373297366500256
        0.008212314081922013
        0.006742220947265425
        0.00549597416771544
        0.0046305707113795955
        0.004008088439249774
        0.003640208872663226
        0.0033964225010390376
        0.0032184072016652535
        ];
    mean_snr = [3.7577007062559596E-4
        3.4634806832516227E-4
        3.229624375608037E-4
        3.0271529852042057E-4
        2.86554145676448E-4
        2.7408153042108023E-4
        2.6391706639443695E-4
        2.551761676872991E-4
        2.4863205791924507E-4
        2.4351809295339405E-4
        2.3949480146124278E-4
        2.360946997763212E-4
        2.331038388339646E-4
        2.3066505981281915E-4
        2.2908254914590965E-4
        2.273503376455066E-4
        2.266558070280672E-4
        2.2569419375110646E-4
        2.2556361546327496E-4
        2.2450445853618433E-4
        2.2425237558780766E-4
        ];
    peak2ave_snr = [365.665707396225
        321.3527511859909
        279.18046837030295
        239.04229184259242
        202.60274566367198
        170.61234918106118
        142.17397997432533
        118.18779757131351
        96.5588791686003
        79.59936845943984
        65.41138673828127
        53.373483019959366
        44.40115062123802
        35.51891472012163
        29.380565984057355
        24.130226666538444
        20.40494669125373
        17.725792754416258
        16.125112598075045
        15.11906714416662
        14.34324439313929
        ];

    

    snr_seq = -5:-1:-25;

    peak_currentSNR = peak_snr(snr_seq==snr);
    mean_currentSNR = mean_snr(snr_seq==snr);
    peak2ave_currentSNR = peak2ave_snr(snr_seq==snr);

    thr = factor*peak2ave_currentSNR*mean_currentSNR.*ones(4096,1);
   

    % figure; plot(pdp); hold on; plot(thr); 



end
