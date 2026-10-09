# 原 MATLAB PRACH 工程的铁路配对 ZC 改造与实现路线

算法的完整数学推导、独立 12 文件工程说明见 [paired_zc_5gr/README.md](../paired_zc_5gr/README.md)。
本文件解释如何逐步移植到**现有波形**，不把序列域验证结果当作 NR 系统级验证。
历史文件尚未覆盖；下面的修改需要在 MATLAB 中逐层验证后应用。

## 一、选择正确的原工程基线

以 `prachNTN_time_freq_recovery` 为移植基线，其 `simAlgo=1` 是同资源配对根叠加。
`prachNTN_pairF0_det` 可保留作原配对检测基线。
不要误选 `PRACH/main.m` 的默认 `simAlgo=2`：那个分支是 PSS 加单根，不是同资源配对 ZC。
`railway_conjugate` 是此前另一个时分重复原型，本次选型不使用它。

推荐改造顺序：固定 F0 和当前采样配置 → 平坦无噪声铁路信道 → AWGN → 块平坦 Rician →
验证抽取滤波器 → 地面 TDL → 分数时延 → 先验失配、多用户及公平性能比较。
在当前设计中，严格公式成立于整数循环时延与块静态平坦路径；TDL 阶段属于后续扩展。

## 二、原系统的三个采样域必须分开

默认 PUSCH 带宽 100 MHz、SCS 30 kHz、F0：

| 量 | 值 | 含义 |
|---|---|---|
| BWP NFFT | 4096 | 不是 PRACH 原始 IFFT 长度 |
| Fs0 | 122.88 MHz | `bwpInfo.SamplingRate` |
| PRACH Δf | 1.25 kHz | `prachConfig.scs_RA*1e3` |
| N0 | 98304 | `prachConfig.nFFT=Fs0/Δf` |
| CP0 | 12672 | F0 原始 CP 3168 经 ×4 过采样 |
| 抽取因子 D | 24 | 当前前端级联 4×3×2 |
| FsRX | 5.12 MHz | `Fs0/D` |
| M | 4096 | `N0/D`，一次有效 PRACH 接收块 |
| Fseq | 1.04875 MHz | `839*Δf`，独立算法的等效序列域 |

原始、接收、序列域时延样点分别满足：

`tau=dRaw/Fs0=dRx/FsRX=dSeq/Fseq`。

因此 `dRaw=D*dRx`、`dSeq=(839/M)*dRx`。
839 点样点约 0.9535 us，4096 点样点约 0.1953 us，原始样点约 8.138 ns。
绝对不能把 `time_comp` 的 839 点下标直接用于 `rxSig(time_comp+1:end,:)`。
所有接口同时返回 `delaySamples`、`delaySeconds`、`sampleRateHz`，明确估计所在域。

## 三、逐文件修改

### 1. `main.m`：场景、真值、总流程

删除铁路入口中的 `SatelliteAltitude`、`MobileAltitude`、`ElevationAngle`、
`SatelliteDopplerShift` 和 `NTN-TDL-C`；不再调用轨道多普勒函数。
速度以 `SpeedKmh` 存储并显式转 m/s，加入基站横向距、沿轨位置及振荡器频偏。
默认 1.9 GHz、350 km/h 是实验配置，最终要与课题采用的 5G-R 要求逐项对应。

使用 `railwayDoppler` 计算真实主径 CFO 及接收机独立先验，不能给接收机传真实 CFO。
时延采用明确的随机接入往返约定 `2R/c`，以秒为唯一场景输入，再按实际 Fs0 转样点。
移除 `delay_samples=200` 的覆盖。保存完整延迟缓存：

```matlab
delayRaw = round(delaySeconds*prachConfig.bwpInfo.SamplingRate);
delayed = [zeros(delayRaw,size(txSig,2)); txSig];
```

不要再 `txSig(1:end-delayRaw,:)` 截断尾部。留出搜索窗口、最大多径时延、滤波延迟和噪声保护区。
AWGN 用基础随机数实现，明确噪声方差的定义，不依赖 `awgn` 默认功率假设。
若测量有效突发功率，不能把 10 ms 空白帧计入信号平均功率。

第一阶段只用平坦信道，避免把正确性与 TDL 前端问题混在一起；
地面 TDL 后续使用官方支持的地面 `nrTDLChannel` 配置，并把主径运动 Doppler、
散射径 Doppler 扩展和公共振荡器 CFO 区分，不能重复加同一个主径偏移。

新接收顺序：运动先验原始波形补偿 → 前端及精确参考 → 配对候选粗定时 →
完整配对相位残差 → CFO 补偿 → 一维二次定时 → 归一化检测。
停止调用旧二维 `time_freq_refinement`。各次 trial 在 true 结构中保存实际时延和 CFO，
接收函数只读 receiver 配置；统计时才将估计与 true 比较。

### 2. `prachSigGen.m`：保持配对、补齐元数据、等能量

保留已有逻辑根表与 `u_ = L_RA-u`。逻辑索引 22 实际映射到 u=1，不能写成根 22 或 23。
保留 `y_=fft(x_)/sqrt(L_RA)`，不能替换成逐点 `conj(y)`。
将叠加分支改成精确归一化：

```matlab
pairFD = y+y_;
pairScale = sqrt(sum(abs(pairFD).^2)/L_RA);
pairFD = pairFD/pairScale;
ifftin = zeros(1,nFFT);
ifftin(firstSC+(1:L_RA)) = pairFD;
ifftout = ifft(fftshift(ifftin))*sqrt(nFFT);
```

这样有用块能量保持 L_RA，与单根分支一致；不能对双根增加 3 dB 后归因于算法收益。
把元数据加入输出：

```matlab
prachConfig.rootU = u;
prachConfig.rootV = u_;
prachConfig.pairFD = pairFD;
prachConfig.pairScale = pairScale;
prachConfig.deltaHz = scs_RA*1e3;
prachConfig.rawFsHz = bwpInfo.SamplingRate;
prachConfig.rawUsefulLength = nFFT;
prachConfig.rawReference = ifftout(:);
```

精估频可以先直接使用 `rawReference` 在**原始采样域**验证，以绕开现有抽取器的误差。
这一验证计算量较大，但参考的 CP、频带映射与发射完全一致，适合短无噪声回归。
单根/配对的公平对照只改变归一化波形，时间、采样率、噪声与天线设置保持一致。

### 3. `cal_bwp_info.m`：保留配置计算，加入一致性检查

保留 BWP NRB/NFFT/Fs 计算，检查带宽和 SCS 确实在支持集合内，拒绝空 `bwIdx`。
进入 PRACH 时检查 `Fs0/nFFT == scs_RA*1e3`、nFFT/抽取因子为整数、有效块长度一致。
先固定当前 100 MHz/30 kHz/F0；不能把未实现的 61.44 MHz 分支当作支持。
本次不自动把 F0 改成短格式，因为这会改变 L、CP、重复数和论文资源开销。

### 4. `prachFrontendProcessing.m`：频带、采样率、群延迟是必修项

当前前端去 CP、把频带下沿移到 DC，然后做实低通抽取。
F0 的 839 根序列频带跨度约 1.04875 MHz，但原通带 Fpass=.54 MHz 是半带宽量；
若不把频带居中，不能假设低通保留完整单边 PRACH 频带。
F3 的频带跨度约 4.195 MHz，不能沿用 F0 通带。先仅验证 F0。

推荐将频带**中心**移到 DC：

```matlab
bandCenterHz = (K*k1+k_bar+(L_RA-1)/2)*prachConfig.deltaHz;
n = (0:size(rxSig,1)-1).';
rxBaseband = rxSig.*exp(-1i*2*pi*bandCenterHz*n/prachConfig.rawFsHz);
```

运动先验补偿是另外一个 CFO 项，不能把已补偿的先验重复扣除。
滤波器通带至少覆盖根频带半宽和**残余** CFO 裕量，阻带满足每级抽取抗混叠要求。
采用可复现的系数与实际频率验证，检查两支路频谱都被完整保留。

现有群延迟按 `nFFT_RA=2048` 转样点，但默认实际块 M=4096，存在因子 2 不一致；
FIR 群延迟也应是 `(length(h)-1)/2`，不是 `length(h)/2`。
若第 i 级输入采样率为 Fs_i、FIR 长度为 Ni：

`tauFilter=sum_i (Ni-1)/(2*Fs_i)`，`delayRx=tauFilter*FsRX`。

保留分数部分，使用匹配参考或分数延迟补偿；不要每级取整后再 `circshift(-2)`。
推荐发射本地参考也通过**相同的下变频、滤波、抽取和窗口**，使接收参考与实际前端一致。
这可以自动携带已知滤波器影响，但绝对传播时延仍须根据窗口原点和已知群延迟解释。

让输出包含：

```matlab
frontend.rawFsHz = Fs0;
frontend.decimation = D;
frontend.rxFsHz = Fs0/D;
frontend.usefulLength = prachConfig.nFFT/D;
frontend.groupDelaySeconds = tauFilter;
frontend.bandOrigin = 'centered';
```

原 `timeData` 为 Nrx×M×Nrep，可用 `blocks=permute(timeData,[2 1 3])` 得 M×Nrx×Nrep。
不要依赖 `squeeze` 在单天线时保持同样维度。

### 5. `prachPowerDelayProfileGen.m`：统一参考和相关功率

停止逐样点 AGC，它可能改变完整配对参考的幅度权重；改用明确整体能量归一化。
不要把天线复相关先加再取平方；用 `P[k]=sum_a abs(C_a[k]).^2`。
删除第二支路的单独 `fftshift(c_)`，否则 c_ 与原 `lag_` 不再对应。
删除 `circshift(-2)` 以及 54458400 等固定功率除数。

若前端采用频带居中，构造一致的 M 点理想参考：

```matlab
bins = mod(-(L_RA-1)/2:(L_RA-1)/2,M)+1;
plusBins = zeros(M,1); minusBins = zeros(M,1); pairBins = zeros(M,1);
plusBins(bins) = prachConfig.fdRootSeq(:);
minusBins(bins) = prachConfig.fdRootSeq_(:);
pairBins(bins) = prachConfig.pairFD(:);
refPlus = ifft(plusBins)*sqrt(M);
refMinus = ifft(minusBins)*sqrt(M);
refPair = ifft(pairBins)*sqrt(M);
```

实际滤波后的参考优先由同一前端获得，上例仅在理想无滤波/已知补偿条件下成立。
如果保留频带下沿到 DC，频率索引应为 `1:L_RA`，参考与接收机必须选择同一个原点。
M 点参考有用能量为 L_RA，平均功率 L_RA/M；不要无声地归一化成 M 能量。

PDP 之外返回 `refPair`、零基相关 lag、FsRX 和所用参考能量，供恢复和检测一致使用。
循环相关适用于合法 CP 内的完整有效块；超 CP、缺失样本或窗口截断时改用线性缓冲模型。

### 6. `time_freq_recovery.m`：受先验支持的整数类消歧

建议替换为新的 `railwayPairedRecovery` 接口，旧入口可以转发，避免名字含义混乱。
839 点模型可直接使用模逆峰方程；M=4096 域则不能照抄整数索引公式。

原因：M/L≈4.882，理想模峰映射是插值周期中的 `kRx≈(M/L)*kSeq`，
传播时延 dRx 未必对应整数 dSeq。简单 `round(kRx*L/M)` 会丢失接收域时延精度。
应在真实接收采样域保留小数峰位、周期折返和滤波延迟，生成有界候选，再用完整参考做验证。
可以先把 M 点数据正确重采样/资源还原到 L 点验证整数模关系，但该阶段只能给粗同步，
最终必须回到 M 点或原始 Fs0 做定时修正，不能宣称 839 点结果具有原始样点精度。

更保守的移植第一版：运动先验将残余限制在远小于 Δf/2 的范围，
在 M 点合法时延窗口比较两支路匹配功率的同位置证据，并输出粗 dRx；
随后完整参考相位估频和一维二次定时。此版不承诺宽 CFO 的精确整数模消歧。
宽先验和任意根的 M 点插值模解另设回归目标，不混入首个可验证版本。

粗恢复输出 `coarseDelayRx`、`priorHz`、`coarseCfoHz`、`valid`，禁止返回无单位的整数再用于切片。

### 7. `time_freq_refinement.m`：停止二维搜索，换相位估频

删除旧调用中的 `fRange=-500:10:500`、`tRange=-15:1:15` 双循环和与算法无关的固定散点图。
可新增 `phaseCFOEstimator.m`，或保留旧文件名仅包装新的函数，参数必须包含完整配对参考和实际 Fs。

M 点通用实现无需知道 ZC 峰斜率：

```matlab
reference = circshift(refPair,coarseDelayRx);
phase = phaseCFOEstimator(rxBlock,reference,frontend.rxFsHz,H);
```

H 建议约 M/4，无模糊范围 `FsRX/(2H)`；M=4096、H=1024 时约 ±2500 Hz。
rxBlock 必须已补偿先验及所选整数 CFO，保持同一个时间原点。
复共轭去扩频自动移除本地参考包含的已知频带相位，不要另扣一次频带偏移。
参考与接收维度为 M×1、M×Nrx。累加逐天线非循环滞后，不跨序列边界。

估计残余 CFO 后只重算一维匹配相关，在合法接收时延窗口选择位置；
可用峰附近插值求分数时延，但插值结果只是估计，需要分数延迟测试，不是精度保证。
如果多径下相位偏差显著，记录为方法局限；不要用真实信道参数校正后宣称算法有效。

### 8. `prachDetection.m`：完整流程标定代替真值 SNR 查表

改接口为 `det=prachDetection5GR(recovery,threshold)`。
删除 `snr_seq`、真实 SNR 对应峰值表和固定 4096 点门限。
不要只取 `pdp(:,1)` 丢弃第二根，判决必须使用与候选解算一致的配对/完整参考统计量。
采用归一化匹配质量 gamma，并用同一前端、候选筛选、相位估计、一维精定时和先验门控的 H0 标定。
系统级标定还需经过真实抗混叠滤波器，滤波有色噪声不能直接用未滤波白噪声分位数。
训练 H0、验证 H0 和有信号试验随机样本分开，报告目标与实测 Pfa、试验数和置信区间。

### 9. `dopplerShiftCircularOrbit.m`：铁路入口不再调用

历史卫星对照可以保留文件。铁路入口新增并使用 `railwayDoppler.m`；
不要仅将卫星高度设零后继续使用轨道模型。
基站附近 CFO 变号，沿轨位置 x 和有符号 v 均需要跨越最近点的场景回归。

### 10. 统计、测试和已有结果文件

替换 `cdfplot`/`ksdensity` 依赖为排序经验 CDF，取消自动覆盖 `tErrors.mat`、`fErrors.mat`。
有信号 trial 必须同时记检测、捕获、无效候选和条件 RMSE；不能只画成功样本误差图。
TDL 中最强径与首径要分别定义，时间误差真值不能含糊地混用传播时延和相关最大峰。

## 四、每一层必须通过的检查

| 层 | 必需检查 | 未通过时不能声称 |
|---|---|---|
| 839 等效模型 | 根 1/23/129，正负 CFO，0/1/边界时延，半整数 CFO，相反天线相位 | 模关系和估频实现正确 |
| 原始 NR 映射 | 本地发射/reference 全链路匹配；pairFD 能量等于单根；CP 与窗口一致 | 原波形已经接通 |
| 抽取前端 | 两根频谱完整保留；滤波群延迟测量与理论一致；无噪声 CFO 正负号 | 4096 接收域和真实时延单位一致 |
| AWGN/Rician | 独立 H0、正确检测/捕获率、条件误差和拒检率 | 低 SNR 性能改善 |
| 地面多径 | 首径/最强径区分，不同路径 Doppler、分数时延、超 CP 场景 | 适用于完整铁路传播环境 |
| 公平对比 | 同总能量/时长/天线/先验/搜索支持/目标 Pfa，并报告实测 Pfa | 优于单根、原配对或 JTFR |

运动先验消融应先保持额外信息预算公平；如果只有新方法获知精确速度，不能把收益全部算作波形创新。
方法潜在研究贡献是铁路先验误差下的模消歧、配对参考相位估计和低复杂度定时链路权衡，
不是把已有的配对 ZC、相位差估频重新命名为全新理论。需要检索具体 JTFR 文献后才能做逐项比较。

## 五、现阶段完成情况

已实现可独立调用的 `paired_zc_5gr` 12 文件工程，提供确定性测试、H0 标定、
平坦 AWGN/Rician 仿真、统计和绘图；已写完整公式及本逐文件迁移说明。
MATLAB 文件通过静态语法解析；独立数值核对通过 72 个完整恢复组合。
当前未执行 MATLAB/Octave，也未改写或验证历史 NR 前端，因此不能宣称系统级移植完成。
遵照此前要求，未配置 MATLAB。
