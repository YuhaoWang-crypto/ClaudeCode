# AlphaZero 国际象棋（浏览器版）

把 [zeta36/chess-alpha-zero](https://github.com/zeta36/chess-alpha-zero) 的方法移植到浏览器里，做成一个可以直接对弈的网页。
它用的是原项目训练好的神经网络权重，所有计算都在你的浏览器里完成（有 GPU 时用 WebGL）。

## 在线玩

打开 **<https://raw.githack.com/YuhaoWang-crypto/ClaudeCode/claude/pensive-turing-r25ie6/alphazero_chess/web/index.html>** 即可（这是从 GitHub 仓库直接托管的页面，不需要登录）。第一次打开要下载约 22 MB 的权重，之后浏览器会缓存。

页面不依赖任何第三方 CDN：TensorFlow.js 和 chess.js 都放在 `web/vendor/` 里，权重存成文本文件，字体加载失败也会用系统字体照常显示。

## 本地运行

```bash
cd alphazero_chess/web
```

```bash
python3 -m http.server 8000
```

然后打开 <http://localhost:8000>。不能直接双击 `index.html`，因为浏览器不允许从 `file://` 读取权重文件。

- 点击或拖动棋子走棋；兵到底线会弹出升变选择
- **执白 / 执黑**：开新局并选择你的颜色
- **AI 思考深度**：直觉（不搜索）/ 快速约 1 秒 / 标准约 3 秒 / 深思约 8 秒。按时间计算，所以电脑越快，同样时间里搜索得越多、下得越强；页面会显示本机每秒能搜索多少次
- **后台推演**：轮到你思考时，AI 也在继续搜索；你走完后，它会直接沿用这部分结果
- **提示**：让 AI 替你搜索一次，用箭头标出建议走法
- **悔棋**：退回到你上一次走棋之前
- 右侧面板实时显示 AI 的搜索：每个候选步的搜索次数 N、网络直觉概率 P、搜索得到的胜率 Q
- 对局自动保存在浏览器里，刷新页面可以继续

## 学到并实现的技术

| 部分 | 原项目 | 这里的实现 |
|---|---|---|
| 局面编码 | `env/chess_env.py` 的 `canon_input_planes` | `web/engine.js` 的 `inputPlanes`：18 个 8×8 平面（12 层棋子、4 层易位权、50 步计数、吃过路兵格）。轮到黑方时整盘翻转，网络总是从走棋方视角看棋 |
| 走法编码 | `config.py` 的 `create_uci_labels` | 同样的 1968 个走法标签；黑方走棋时把网络输出按棋盘上下镜像映射回真实走法 |
| 神经网络 | `agent/model_chess.py` | 7 个残差块 × 256 个 3×3 卷积核；策略头输出 1968 个走法概率，价值头输出 [−1, 1] 的局面评估。用 TensorFlow.js 运行 |
| 搜索 | `agent/player_chess.py` | PUCT 蒙特卡洛树搜索：选择 `Q + c·P·√N/(1+n)` 最大的分支，叶子交给网络评估，价值逐层取反回传，最后走访问次数最多的一步。原项目用 16 个线程 + 虚拟损失并行；这里改成每批 8 个叶子一起送进 GPU，虚拟损失为 1（实测与逐个搜索的访问分布几乎一致，速度快一倍） |

权重转换：`tools/convert_weights.py` 读取原项目的 `model_best_weight.h5`，把 BatchNorm 合并进卷积，存成 float16，再拆成三个 base64 文本文件（每个约 8 MB，任何静态托管都能直接提供）。

## 在原项目之上加强了什么

每一项都先用 `tests/match.mjs` 做引擎对引擎的对局（同一网络、每步 100 次模拟、12 个随机开局各下黑白两盘），确认确实变强才保留：

| 改动 | 作用 | 对局结果 |
|---|---|---|
| 首次探索惩罚（FPU = 0.3） | 原项目把没搜过的走法当成“平局”（Q = 0），局面差时会浪费大量搜索在烂棋上。改为“比当前局面估值低一点”，搜索集中在好棋上 | 对不加 FPU 的同一引擎 FPU_RESULT |
| 搜索树复用 | 上一步搜过的分支留着继续用，相当于免费多搜 | 与杀棋求解器合并测试 |
| 杀棋求解器 | 搜索中一旦证明某步是将杀就直接走，证明会被将杀的走法就避开 | 与搜索树复用合并测试 |
| 按时间思考 | 难度按秒数而不是固定次数，GPU 快的电脑同样时间能搜几百上千次 | — |
| 后台推演 | 你思考时 AI 继续搜索，你走完后直接沿用 | — |

全部合在一起，对原项目的搜索方式（同样每步 100 次模拟）：PACKAGE_RESULT。

另外试过把 FPU 调到 0.6，反而略差（对 0.3 是 1 胜 7 和 4 负），所以保留 0.3。CPUCT_LINE

## 有多强

原作者用约 1 万盘人类对局做监督学习，估计它在每步 1200 次模拟时约 1200 等级分。棋力主要取决于每步能搜索多少次：有独立显卡的电脑在“标准”（3 秒）下通常能搜几百到上千次，比只有集成显卡或用 CPU 的设备强得多。
它能看出简单的杀棋（例如 Qxf7#、底线杀），但中残局里仍会犯错，适合初学者到中级玩家练习。网络本身没有重新训练，想要质的提升需要更强的网络（例如 Leela Chess Zero 的权重），那是另一项更大的工作。

## 校验

`tools/az_reference.py` 用 numpy 逐层复现了原 Keras 模型。`tests/vectors.json` 是它在 8 个局面上的输出，`tests/check_engine.mjs` 检查浏览器引擎与之一致（价值误差 < 0.001，概率误差 < 0.001），并检查搜索能找到白方和黑方的一步杀。

```bash
npm install chess.js@1.4.0 @tensorflow/tfjs@4.22.0
```

```bash
node alphazero_chess/tests/check_engine.mjs
```

重新生成权重（需要原项目的 h5 文件和 `pip install h5py numpy`）：

```bash
python3 alphazero_chess/tools/convert_weights.py chess-alpha-zero/data/model/model_best_weight.h5 alphazero_chess/web/model/weights --vectors alphazero_chess/tests/vectors.json
```

## 文件

```
alphazero_chess/
  web/index.html        页面和样式
  web/app.js            棋盘、走棋交互、对局流程、分析面板
  web/engine.js         编码 + 神经网络 + MCTS（可在浏览器和 Node 中运行）
  web/model/            转换后的权重（float16，base64 文本，约 22 MB）
  web/vendor/           TensorFlow.js 4.22.0 和 chess.js 1.4.0（与 npm 发布版逐字节相同）
  tools/az_reference.py numpy 参考实现
  tools/convert_weights.py
  tests/check_engine.mjs, tests/vectors.json
  tests/match.mjs       引擎对引擎的对局，用来检验某项改动是否真的变强
```

## 致谢与许可

模型结构、编码方式和权重来自 zeta36/chess-alpha-zero（MIT 许可，Copyright (c) 2017 Samuel Gravan）。
走法规则用 [chess.js](https://github.com/jhlywa/chess.js)，神经网络推理用 [TensorFlow.js](https://www.tensorflow.org/js)，两者都随页面一起打包在 `web/vendor/`。
