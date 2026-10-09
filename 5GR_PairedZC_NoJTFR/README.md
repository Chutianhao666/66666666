# 5G-R 配对 ZC：运动先验、模约束粗定时与去扩频相位估频（无 JTFR）

本目录是可独立运行的 **839 点序列域等效 MATLAB 工程**，共 11 个 `.m` 文件与本说明。
保持原代码 `simAlgo=1` 的**同资源双根叠加**，不采用之前 `railway_conjugate` 目录的时分重复结构。
不调用 JTFR，不进行连续时频二维网格搜索。只需要基础 MATLAB R2016b 或更新版本。
这不是完整的标准 PRACH 波形：理想 CP、资源映射、抽取和群延迟均未在第一层模拟。
共轭双根叠加本身也是研究扩展，不能声称无需协议修改即可部署。

## 运行

在 MATLAB 中进入本目录：

```matlab
test_5GR
results = main_5GR;
```

无绘图的小规模诊断：

```matlab
cfg = defaultConfig();
cfg.makePlots = false;
cfg.epochs = 20;
results = main_5GR(cfg);
```

结果保存在返回结构，不自动写文件。要保存请自行指定仓库外的路径，防止覆盖旧结果。
不要 `addpath(genpath(...))` 混入其他场景；新旧目录均有 `railwayChannel.m`，MATLAB 路径顺序会影响解析。

## 铁路几何与接收机先验

令基站投影到轨道的位置为原点，列车沿轨道坐标为 x，横向间距为 b，沿轨速度为 v：

`R=sqrt(x^2+b^2)`，`fD=-(fc/c)*v*x/R`。

x>0、v<0 表示接近，产生正频偏。速度配置为 km/h，必须先除以 3.6。
往返随机接入到达偏移为 `tau=2R/c`；普通单程信道时延则为 `R/c`，两者不能混用。
默认 `fc=1.9 GHz`、350 km/h、x=3 km、b=50 m：主径多普勒约 616 Hz。
默认振荡器频偏 +150 Hz，总 CFO 约 766 Hz。这些是课题实验假设，不是强制铁路规范值。

接收机 `fPrior` 使用**另行配置**的 `priorMotion`（默认 330 km/h）和 `priorOscillatorHz`（默认 0），
从不读取 `trueMotion`、`oscillatorHz` 或 `truth`。默认先验误差约 185 Hz，设不确定度 B=300 Hz。
测速/位置先验可来自已有车载定位、列控测速或前一时刻的无线测量；可用性、测量时延及误差分布须单独研究。
不能假设初次随机接入时基站天然知道每辆列车的精确位置或速度。

消融“无运动先验”：将 `priorMotion.speedKmh=0`、`priorOscillatorHz=0`，并扩大
`priorUncertaintyHz` 到覆盖实际频偏的值（例如 1000 Hz）。必须重新标定门限。
违反 B 时允许拒检，不能偷偷用真实频偏回退。先验误差下的拒检率由 `invalidRate` 单独输出。

## 1. 配对 ZC 及功率归一化

奇数 L=839，实际根 u 与 L 互素，n=0,...,L-1：

`x_u[n]=exp(-j*pi*u*n*(n+1)/L)`，`x_(L-u)[n]=conj(x_u[n])`。

保持叠加结构：`s[n]=(x_u[n]+x_(L-u)[n])/A`，
`A=sqrt(mean(abs(x_u+x_(L-u)).^2))`，使每样点平均发射功率为 1。
不能简单把两个根相加后与单根作等功率比较，也不能默认交叉项恰好为零而只除以 sqrt(2)。

原代码逻辑索引 22 对应**实际根 1**，新工程默认 root=1。23 和 129 仅用于一般根的回归检查。
频谱关系为 `Y_(L-u)[k]=conj(Y_u[(-k) mod L])`，不是逐点 `conj(Y_u[k])`。

## 2. 接收等效模型与先验补偿

静态平坦路径、整数循环时延 d 下，第 a 根接收天线：

`r_a[n]=h_a*s[(n-d) mod L]*exp(j*2*pi*f*n/Fseq)+w_a[n]`。

`Fseq=L*Delta_fRA`，默认 PRACH 间隔 `Delta_fRA=1250 Hz`，所以 `Fseq=1.04875 MHz`。
一个序列样点约 0.9535 us，最大 6 km 往返时延对应搜索 d=0,...,42。
该循环延迟假设等价于合法 CP 后完整有效块，**不是**时域延迟接收缓存的直接实现。

先用运动先验：`r0_a[n]=r_a[n]*exp(-j*2*pi*fPrior*n/Fseq)`。
残余归一化频偏 `epsilon=(f-fPrior)/Delta_fRA=m+delta`，其中 m 为整数类。

## 3. 配对相关峰与模解算

采用与 `ifft(fft(r).*conj(fft(x)))` 一致的相关定义：

`C_u[k]=sum_n r0[n]*conj(x_u[(n-k) mod L])`。

先忽略另一根交叉项，展开两个二次相位：

`C_u[k]=h*exp(j*phi(k,d))*D_L(epsilon-u*(k-d))`，
`D_L(z)=sum_n exp(j*2*pi*z*n/L)`。

其闭式是 `exp(j*pi*z*(L-1)/L)*sin(pi*z)/sin(pi*z/L)`（零分母用极限）。
整数 CFO m 时主自相关峰满足：

`k_u=d+inv(u)*m (mod L)`，`k_(L-u)=d-inv(u)*m (mod L)`。

因 L 为奇数，`inv(2)=(L+1)/2=420`：

`d=inv(2)*(k_u+k_(L-u)) (mod L)`，
`m=u*inv(2)*(k_u-k_(L-u)) (mod L)`。

m 需转换成有符号代表并按物理频偏界筛选。一般根必须使用模逆，不能把根号直接当作线性斜率。
非整数 CFO 下存在谱泄漏，配对叠加还存在互根 Gauss 交叉项，因此这些是**整数类的候选关系**，
不是非整数 CFO、多径或采样率变化下的连续精确等式。

直接各取一个全局最大峰有无噪声反例：u=1、d=42、epsilon=0.49 时，两个峰可能为 42 与 41，
直接模半和给错误时延 461。新实现不盲目平均全局峰，而按物理支持反向生成候选。

具体做法：`|m| <= floor(B/Delta_fRA+1/2)`，`0<=d<=dmax`；
由 d,m 算出合法 `k_u,k_(L-u)`，比较该配对的功率证据：

`J(d,m)=sqrt(P_u[k_u]*P_(L-u)[k_(L-u)])/(1+(m*Delta_fRA/B)^2)`，
`P_u[k]=sum_a abs(C_(u,a)[k])^2`。

分母是显式先验偏好启发式，不宣称严格 MAP。默认 B=300 Hz 时只有 m=0，退化成一个受约束的
一维时延曲线；B 扩大时才枚举有限整数模类。这是解析产生的离散消歧候选，不是 JTFR 的连续频率网格。
若先验过宽，整数类增多，复杂度和误选概率也会上升。

## 4. 完整配对去扩频与相位差估计

选粗候选 d0,m0 后先补偿 `f0=fPrior+m0*Delta_fRA`，得到 r1。
本地参考必须是**完整叠加配对** `s_d[n]=s[(n-d0) mod L]`。

`z_a[n]=r1_a[n]*conj(s_d[n])`。

正确 d0、静态平坦路径和无噪声时：

`z_a[n]=h_a*abs(s_d[n])^2*exp(j*2*pi*fres*n/Fseq)`。

取 H=209，用**非循环**时间滞后累加：

`Q_H=sum_a sum_(n=0)^(L-H-1) z_a[n+H]*conj(z_a[n])`，
`fres_hat=Fseq/(2*pi*H)*angle(Q_H)`。

幅度系数为 `sum_a |h_a|^2 sum_n |s_d[n+H]|^2 |s_d[n]|^2 >= 0`，
所以相位只保留残余频偏。频偏无模糊范围为 `|fres|<Fseq/(2H)`，默认约 ±2509 Hz。
H 增大可提升相位对频偏的敏感度，但缩小无模糊范围且增加对时变信道的敏感性。

实现等价地直接累加
`r1[n+H]*conj(r1[n])*conj(s_d[n+H])*s_d[n]`，不除以 s 的近零样点。
先逐天线作乘积，再相加；未知接收天线初相位被消去。
不能将 n+H 越过序列末端循环回绕，因为非整数 CFO 的相位斜坡不是 L 周期。
只用 x_u 做参考会保留另一根干扰项，这也是原精估计单根参考不适用于叠加配对的原因。

## 5. 二次定时、判决与复杂度

补偿总频偏后，用完整配对参考做一次 FFT 循环相关，仅在合法时延窗口选最大值：

`gamma[d]=sum_a |sum_n rComp_a[n]*conj(s[(n-d) mod L])|^2 / (Es*Er)`，
`Es=sum_n |s[n]|^2`，`Er=sum_a sum_n |r_a[n]|^2`。

由 Cauchy-Schwarz，理想数学值 gamma 在 [0,1]。正确无噪声平坦信道为 1。
再按新的 d 做一轮相位估计与一维定时，默认总共两轮，固定迭代次数保证标定和检测流程一致。
估计总 CFO 超出先验支持则判无效，最终检测条件是有效且 `max(gamma)>threshold`。
此时延指标是平坦路径到达位置；扩展到多径后最大峰可能为最强径，不能自动称为首径估计。

FFT 相关约 `O(Nrx*L*log L)`，模候选约 `O((2M+1)*(dmax+1))`，
每次相位估计 `O(Nrx*L)`。没有 `Ntime*Nfrequency` 次重复相关的矩形搜索。
效率提升不等于性能提升；低 SNR 错粗定时会污染相位，先验失配会增加拒检。

## 6. 门限与统计口径

`calibrateThreshold` 对**完整**接收流程做 N=500 次独立无信号试验，保存每次最终峰值。
排序后取秩 `ceil((N+1)*(1-targetPfa))`。随后使用另外 500 次 H0 测试虚警率并报告 Wilson95 区间。
这只是白噪声下的经验标定，不是解析 CFAR；配置、搜索支持或迭代次数改变后须重新标定。
论文研究建议至少 10000 次 H0，目标很低时还需更多；不能把 0 次虚警写成真实 Pfa=0。

信道噪声定义为名义输入**每复样点 SNR**，平均发射样点功率 1、`sigma2=10^(-SNR/10)`。
它不是整个前导能量 Es/N0；长度 839 的总前导能量额外带来约 29.24 dB。
Rician 系数每个块内恒定，统计平均功率为 1，不对每次衰落重新缩放噪声。

输出 Pd、正确时延捕获率、已检测样本的定时/CFO 条件 RMSE、无效候选率及已检测样本数。
条件 RMSE 不计漏检，必须和捕获率一起报告。无人检测时 RMSE 为 NaN，不应记成 0。
绘图只用基础 MATLAB 的排序、折线，不依赖 `cdfplot` 或 `ksdensity` 工具箱。

## 文件及接口

| 文件 | 输入 → 输出 | 职责 |
|---|---|---|
| `defaultConfig.m` | 无 → cfg | 真值、先验、物理单位、统计参数集中配置 |
| `pairedZCGenerate.m` | cfg → s,ref | 同资源双根、完整参考、精确功率归一化 |
| `railwayDoppler.m` | motion,fc,c → fHz,R | 轨旁几何，接近/远离频偏符号 |
| `railwayChannel.m` | s,cfg,SNR → L×Nrx 的 rx,truth | AWGN 或块平坦 Rician，理想循环时延 |
| `pairedZCCorrelation.m` | rx,ref → 双支路复相关及 L×2 功率 | 统一零基 lag，非相干天线合并 |
| `railwayPairedRecovery.m` | rx,cfg,ref → recovery | 先验补偿、模候选、相位估频、一维二次定时 |
| `phaseCFOEstimator.m` | rx,完整参考,Fs,H → 相位频偏、相干度、有效性 | 非循环加权滞后，不做频率搜索 |
| `prachDetection5GR.m` | recovery,门限 → detected,时延,CFO | 保留无效/未检测结果 NaN |
| `calibrateThreshold.m` | cfg,ref → threshold,H0分数 | 完整流程噪声标定 |
| `main_5GR.m` | 可选 cfg → results | 独立 H0、Monte Carlo、统计及绘图 |
| `test_5GR.m` | 无 | 确定性回归检查，不代替统计性能评测 |

## 验证边界与下一步

当前 **没有安装或执行 MATLAB/Octave**。11 个 MATLAB 文件通过 MISS_HIT 语法解析。
独立 NumPy 实现通过 60 个整数循环时延/正负残余 CFO/三根/相反天线相位组合；
另通过 12 个 ±0.49、±0.51 子载波间隔附近的半整数 CFO 案例。
这些不证明 MATLAB 首次执行兼容性，也不证明 NR 前端、多径或分数时延性能。

独立 NumPy 小样本诊断：500 次 H0 标定，另 500 次 H0 出现 7 次虚警（约 1.4%）；
每 SNR 100 次 Rician 诊断在 -20/-15/-10/-5/0 dB 检出 20/51/81/93/100 次，
这只是另一实现的诊断数据，**不是 MATLAB 输出，不是与 JTFR 的公平性能比较**。
这些结果暴露低 SNR 的先验拒检和相位估计问题，尚需进一步研究。

扩展到铁路多径、快速衰落、分数延迟、隧道、非白噪声或多用户之前，应先完成本地 `test_5GR`。
本推导的严格相位关系只保证静态平坦路径与正确整数时延；多径相位混合及首径问题尚未解决。
原始工程逐文件移植、采样率和滤波检查见仓库的 `PAIRED_ZC_5GR_DESIGN.md`。

配对根、运动补偿、去扩频相位估频分别都有既有研究基础，不能直接宣称全新理论。
建议课题定位为“铁路运动先验约束下的配对 ZC 分阶段同步”，用无先验、单根、
无相位修正等消融研究贡献；各方案保持总发射能量、时长、天线数、搜索支持与虚警率一致。
