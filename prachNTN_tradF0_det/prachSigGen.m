function [txSig, prachConfig] = prachSigGen(simConfig)
mu = log2(simConfig.puschSCS/15);
format = simConfig.prachFormat;
zeroCorrelationZoneConfig = simConfig.prachZeroConfigZone;
rootSequenceIndex = simConfig.prachRootSeqIdx;
scs_BWP = simConfig.puschSCS;
nTx = simConfig.nTx;

switch format
    case 'F0'
        L_RA = 839;
        scs_RA = 1.25;
        Nu = 24576;
        Nrep = 1;
        NCP_RA = 3168;
    case 'F1'
        L_RA = 839;
        scs_RA = 1.25;
        Nu = 24576;
        Nrep = 2;
        NCP_RA = 21024;
    case 'F2'
        L_RA = 839;
        scs_RA = 1.25;
        Nu = 24576;
        Nrep = 4;
        NCP_RA = 4688;
    case 'F3'
        L_RA = 839;
        scs_RA = 5;
        Nu = 6144;
        Nrep = 4;
        NCP_RA = 3168;
    case 'A1'
        L_RA = 139;
        scs_RA = 15*2^mu;
        Nu = 2048/2^mu;
        Nrep = 2;
        NCP_RA = 288/2^mu;
    case 'A2'
        L_RA = 139;
        scs_RA = 15*2^mu;
        Nu = 2048/2^mu;
        Nrep = 4;
        NCP_RA = 576/2^mu;
    case 'A3'
        L_RA = 139;
        scs_RA = 15*2^mu;
        Nu = 2048/2^mu;
        Nrep = 6;
        NCP_RA = 864/2^mu;
    case 'B1'
        L_RA = 139;
        scs_RA = 15*2^mu;
        Nu = 2048/2^mu;
        Nrep = 2;
        NCP_RA = 216/2^mu;
    case 'B2'
        L_RA = 139;
        scs_RA = 15*2^mu;
        Nu = 2048/2^mu;
        Nrep = 4;
        NCP_RA = 360/2^mu;
    case 'B3'
        L_RA = 139;
        scs_RA = 15*2^mu;
        Nu = 2048/2^mu;
        Nrep = 6;
        NCP_RA = 504/2^mu;
    case 'B4'
        L_RA = 139;
        scs_RA = 15*2^mu;
        Nu = 2048/2^mu;
        Nrep = 12;
        NCP_RA = 936/2^mu;
    case 'C0'
        L_RA = 139;
        scs_RA = 15*2^mu;
        Nu = 2048/2^mu;
        Nrep = 1;
        NCP_RA = 1240/2^mu;
    case 'C2'
        L_RA = 139;
        scs_RA = 15*2^mu;
        Nu = 2048/2^mu;
        Nrep = 4;
        NCP_RA = 2048/2^mu;
end

gZCconfig_Ncs_139 = [0, 2, 4, 6, 8, 10, 12, 13, 15, 17, 19, 23, 27, 34, 46, 69];
gZCconfig_Ncs_839_Fmt012 = [0, 13, 15, 18, 22, 26, 32, 38, 46, 59, 76, 93, 119, 167, 279, 419];
gZCconfig_Ncs_839_Fmt3 = [0, 13, 26, 33, 38, 41, 49, 55, 64, 76, 93, 119, 139, 209, 279, 419];

switch scs_RA
    case 1.25
        N_CS = gZCconfig_Ncs_839_Fmt012(zeroCorrelationZoneConfig+1);
    case 5
        N_CS = gZCconfig_Ncs_839_Fmt3(zeroCorrelationZoneConfig+1);
    case 15*2^mu
        N_CS = gZCconfig_Ncs_139(zeroCorrelationZoneConfig+1);
end

numCv = floor(L_RA/N_CS);
nRootSeq = ceil(64/numCv);

d_PSS = [...
      1  -1  -1   1  -1  -1  -1  -1   1   1  -1  -1  -1   1   1  -1   1 ...
     -1   1  -1  -1   1   1  -1  -1   1   1   1   1   1  -1  -1   1  -1 ...
     -1   1  -1   1  -1  -1  -1   1  -1   1   1   1  -1  -1   1   1  -1 ...
      1   1   1  -1   1   1   1   1   1   1  -1   1   1  -1   1   1  -1 ...
     -1   1  -1   1   1  -1  -1  -1  -1   1  -1  -1  -1   1   1   1   1 ...
     -1  -1  -1  -1  -1  -1  -1   1   1   1  -1  -1  -1   1  -1  -1   1 ...
      1   1  -1   1  -1   1   1  -1   1  -1  -1  -1  -1  -1   1  -1   1 ...
     -1   1  -1   1   1   1   1  -1].';

n2 = 2;
seq = d_PSS(1 + mod(43*n2 + (0:126),127));
seq = fft(seq,127)/sqrt(127);
seq = seq./abs(seq);

u = logical_idx_to_seq_number_839(rootSequenceIndex+1);
x = exp(-1i*pi*u*(0:L_RA-1).*(1:L_RA)/L_RA);
y = fft(x,L_RA)/sqrt(L_RA);
y = y./abs(y);


u_ = L_RA - u;
x_ = exp(-1i*pi*u_*(0:L_RA-1).*(1:L_RA)/L_RA);
y_ = fft(x_,L_RA)/sqrt(L_RA);
y_ = y_./abs(y_);



[k_bar, N_RB_RA] = PRACH_subcarrier_offset(L_RA, scs_RA, scs_BWP);

K = scs_BWP/scs_RA;
bwpInfo = cal_bwp_info(simConfig);
nFFT = bwpInfo.NFFT*K;
os_rate = nFFT/Nu;
Nu = Nu*os_rate;
NCP_RA = NCP_RA*os_rate;

N_BWP_start = 0;
n_RA_start = 0;
n_RA = 0;

k0 = 0;
k1 = k0 + 12*(N_BWP_start - 0) + 12*n_RA_start + 12*n_RA*N_RB_RA - bwpInfo.NRB*12/2;
firstSC = nFFT/2 + K*k1 + k_bar;

if simConfig.simAlgo
    y_all = y + y_;
    ifftin = zeros(1,nFFT);
    ifftin(firstSC + (1:L_RA)) = y_all;
    ifftout = ifft(fftshift(ifftin))*sqrt(nFFT);
else
    ifftin = zeros(1,nFFT);
    ifftin(firstSC + (1:L_RA)) = y;
    ifftout = ifft(fftshift(ifftin))*sqrt(nFFT);
end


prach_sig = [ifftout(end-NCP_RA+1:end), repmat(ifftout,1,Nrep)];
if nTx > 1
    prachSig =  zeros(size(prach_sig, 2), nTx);
    for ii=1:nTx
        prachSig(:,ii) = prach_sig.';
    end
else
    prachSig =  prach_sig.';
end

Tc = 1/(480e3*4096);                                    % Maximum sampling rate in NR
Ts = 1/(15e3*2048);                                     % Baseline sampling rate in LTE
Kappa = Ts/Tc;                                          % Ratio between Tc and Ts = 64
SamplingRate = bwpInfo.SamplingRate;
cpMode = bwpInfo.CyclicPrefix;
Ds = (1/Tc)/SamplingRate;                               % Down-sampling rate compared to maximum sampling rate
if strcmpi(cpMode,'Normal')                             % CP length
    Ncp = 144*Kappa/2^mu/Ds*ones(1,14);
    Ncp(1:7*2^mu:end) = Ncp(1:7*2^mu:end) + 16*Kappa/Ds;
else
    Ncp = 512*Kappa/2^mu/Ds*ones(1,12);
end
len_subframe = (sum(Ncp)+bwpInfo.NFFT*14)*2^mu;
startidx = 0;
len=len_subframe*10;
if startidx+length(prachSig) >len
    txSig = zeros(startidx+length(prachSig), nTx);
else
    txSig = zeros(len, nTx);
end

txSig(startidx+1:startidx+size(prachSig,1),:) = prachSig;  % data within a radio frame

prachConfig.mu = mu;
prachConfig.bwpInfo = bwpInfo;
prachConfig.nRootSeq = nRootSeq;
prachConfig.fdRootSeq = y;
prachConfig.fdRootSeq_ = y_;
prachConfig.L_RA = L_RA;
prachConfig.scs_RA = scs_RA;
prachConfig.Nu = Nu;
prachConfig.Nrep = Nrep;
prachConfig.NCP_RA = NCP_RA;
prachConfig.N_CS = N_CS;
prachConfig.format = format;
prachConfig.nFFT = nFFT;
prachConfig.K = K;
prachConfig.k1 = k1;
prachConfig.k0 = k0;
prachConfig.k_bar = k_bar;
prachConfig.firstSC = firstSC;
prachConfig.numCv = numCv;


end

%% Auxilary Functions

function u = logical_idx_to_seq_number_839(index)
%
% TS38.211 Table 6.3.3.1-3: Mapping L_RA = 839
%
table = [
    129, 710, 140, 699, 120, 719, 210, 629, 168, 671,  84, 755, 105, 734,  93, 746,  70, 769,  60, 779 ...
    2, 837,   1, 838,  56, 783, 112, 727, 148, 691,  80, 759,  42, 797,  40, 799,  35, 804,  73, 766         ...
    146, 693,  31, 808,  28, 811,  30, 809,  27, 812,  29, 810,  24, 815,  48, 791,  68, 771,  74, 765         ...
    178, 661, 136, 703,  86, 753,  78, 761,  43, 796,  39, 800,  20, 819,  21, 818,  95, 744, 202, 637         ...
    190, 649, 181, 658, 137, 702, 125, 714, 151, 688, 217, 622, 128, 711, 142, 697, 122, 717, 203, 636         ...
    118, 721, 110, 729,  89, 750, 103, 736,  61, 778,  55, 784,  15, 824,  14, 825,  12, 827,  23, 816         ...
    34, 805,  37, 802,  46, 793, 207, 632, 179, 660, 145, 694, 130, 709, 223, 616, 228, 611, 227, 612         ...
    132, 707, 133, 706, 143, 696, 135, 704, 161, 678, 201, 638, 173, 666, 106, 733,  83, 756,  91, 748         ...
    66, 773,  53, 786,  10, 829,   9, 830,   7, 832,   8, 831,  16, 823,  47, 792,  64, 775,  57, 782         ...
    104, 735, 101, 738, 108, 731, 208, 631, 184, 655, 197, 642, 191, 648, 121, 718, 141, 698, 149, 690         ...
    216, 623, 218, 621, 152, 687, 144, 695, 134, 705, 138, 701, 199, 640, 162, 677, 176, 663, 119, 720         ...
    158, 681, 164, 675, 174, 665, 171, 668, 170, 669,  87, 752, 169, 670,  88, 751, 107, 732,  81, 758         ...
    82, 757, 100, 739,  98, 741,  71, 768,  59, 780,  65, 774,  50, 789,  49, 790,  26, 813,  17, 822         ...
    13, 826,   6, 833,   5, 834,  33, 806,  51, 788,  75, 764,  99, 740,  96, 743,  97, 742, 166, 673         ...
    172, 667, 175, 664, 187, 652, 163, 676, 185, 654, 200, 639, 114, 725, 189, 650, 115, 724, 194, 645         ...
    195, 644, 192, 647, 182, 657, 157, 682, 156, 683, 211, 628, 154, 685, 123, 716, 139, 700, 212, 627         ...
    153, 686, 213, 626, 215, 624, 150, 689, 225, 614, 224, 615, 221, 618, 220, 619, 127, 712, 147, 692         ...
    124, 715, 193, 646, 205, 634, 206, 633, 116, 723, 160, 679, 186, 653, 167, 672,  79, 760,  85, 754         ...
    77, 762,  92, 747,  58, 781,  62, 777,  69, 770,  54, 785,  36, 803,  32, 807,  25, 814,  18, 821         ...
    11, 828,   4, 835,   3, 836,  19, 820,  22, 817,  41, 798,  38, 801,  44, 795,  52, 787,  45, 794         ...
    63, 776,  67, 772,  72, 767,  76, 763,  94, 745, 102, 737,  90, 749, 109, 730, 165, 674, 111, 728         ...
    209, 630, 204, 635, 117, 722, 188, 651, 159, 680, 198, 641, 113, 726, 183, 656, 180, 659, 177, 662         ...
    196, 643, 155, 684, 214, 625, 126, 713, 131, 708, 219, 620, 222, 617, 226, 613, 230, 609, 232, 607         ...
    262, 577, 252, 587, 418, 421, 416, 423, 413, 426, 411, 428, 376, 463, 395, 444, 283, 556, 285, 554         ...
    379, 460, 390, 449, 363, 476, 384, 455, 388, 451, 386, 453, 361, 478, 387, 452, 360, 479, 310, 529         ...
    354, 485, 328, 511, 315, 524, 337, 502, 349, 490, 335, 504, 324, 515, 323, 516, 320, 519, 334, 505         ...
    359, 480, 295, 544, 385, 454, 292, 547, 291, 548, 381, 458, 399, 440, 380, 459, 397, 442, 369, 470         ...
    377, 462, 410, 429, 407, 432, 281, 558, 414, 425, 247, 592, 277, 562, 271, 568, 272, 567, 264, 575         ...
    259, 580, 237, 602, 239, 600, 244, 595, 243, 596, 275, 564, 278, 561, 250, 589, 246, 593, 417, 422         ...
    248, 591, 394, 445, 393, 446, 370, 469, 365, 474, 300, 539, 299, 540, 364, 475, 362, 477, 298, 541         ...
    312, 527, 313, 526, 314, 525, 353, 486, 352, 487, 343, 496, 327, 512, 350, 489, 326, 513, 319, 520         ...
    332, 507, 333, 506, 348, 491, 347, 492, 322, 517, 330, 509, 338, 501, 341, 498, 340, 499, 342, 497         ...
    301, 538, 366, 473, 401, 438, 371, 468, 408, 431, 375, 464, 249, 590, 269, 570, 238, 601, 234, 605         ...
    257, 582, 273, 566, 255, 584, 254, 585, 245, 594, 251, 588, 412, 427, 372, 467, 282, 557, 403, 436         ...
    396, 443, 392, 447, 391, 448, 382, 457, 389, 450, 294, 545, 297, 542, 311, 528, 344, 495, 345, 494         ...
    318, 521, 331, 508, 325, 514, 321, 518, 346, 493, 339, 500, 351, 488, 306, 533, 289, 550, 400, 439         ...
    378, 461, 374, 465, 415, 424, 270, 569, 241, 598, 231, 608, 260, 579, 268, 571, 276, 563, 409, 430         ...
    398, 441, 290, 549, 304, 535, 308, 531, 358, 481, 316, 523, 293, 546, 288, 551, 284, 555, 368, 471         ...
    253, 586, 256, 583, 263, 576, 242, 597, 274, 565, 402, 437, 383, 456, 357, 482, 329, 510, 317, 522         ...
    307, 532, 286, 553, 287, 552, 266, 573, 261, 578, 236, 603, 303, 536, 356, 483, 355, 484, 405, 434         ...
    404, 435, 406, 433, 235, 604, 267, 572, 302, 537, 309, 530, 265, 574, 233, 606, 367, 472, 296, 543         ...
    336, 503, 305, 534, 373, 466, 280, 559, 279, 560, 419, 420, 240, 599, 258, 581, 229, 610
    ];
u = table(index);
end

function [k_bar, N_RB_RA] = PRACH_subcarrier_offset(L_RA, scs_RA, scs_PUSCH)
% 
% Get PRACH subcarrier offset, TS 38.211 Table 6.3.3.2-1
% 
k_bar_table = [
    % L_RA/scs_RA/scs_PUSCH/N_RB_RA/k_bar
    839  1.25    15    6     7;
    839  1.25    30    3     1;
    839  1.25    60    2   133;
    839     5    15   24    12;
    839     5    30   12    10;
    839     5    60    6     7;
    139    15    15   12     2;
    139    15    30    6     2;
    139    15    60    3     2;
    139    30    15   24     2;
    139    30    30   12     2;
    139    30    60    6     2;
    139    60    60   12     2;
    139    60   120    6     2;
    139   120    60   24     2;
    139   120   120   12     2
    ];

k_bar = NaN;
N_RB_RA = NaN;
warning('off');
for n = 1:size(k_bar_table)
   if k_bar_table(n,1) == L_RA && k_bar_table(n,2) == scs_RA && k_bar_table(n,3) == scs_PUSCH
       N_RB_RA = k_bar_table(n,4);
       k_bar = k_bar_table(n,5);
       break;
   end
end

end

