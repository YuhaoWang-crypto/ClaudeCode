# AlphaZero 国际象棋（浏览器版）

把 [zeta36/chess-alpha-zero](https://github.com/zeta36/chess-alpha-zero) 的方法移植到浏览器里，做成一个可以直接对弈的网页。
它用的是原项目训练好的神经网络权重，所有计算都在你的浏览器里完成（有 GPU 时用 WebGL）。

## 怎么玩

```bash
cd alphazero_chess/web
```

```bash
python3 -m http.server 8000
```

然后打开 <http://localhost:8000>。不能直接双击 `index.html`，因为浏览器不允许从 `file://` 读取权重文件。

- 点击或拖动棋子走棋；兵到底线会弹出升变选择
- **执白 / 执黑**：开新局并选择你的颜色
- **AI 思考深度**：每步的搜索模拟次数（1 / 64 / 200 / 800）。次数越多越强，也越慢；页面会测出本机速度，设备较慢时自动调低默认值
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

权重转换：`tools/convert_weights.py` 读取原项目的 `model_best_weight.h5`，把 BatchNorm 合并进卷积，存成 float16，并拆成两个不超过 9 MB 的文件（`web/model/`）。

## 有多强

原作者用约 1 万盘人类对局做监督学习，估计它在每步 1200 次模拟时约 1200 等级分。网页里默认 200 次模拟，所以会弱一些；在本地测试中，200 次模拟的版本对只用网络直觉（1 次模拟）的版本 2 战 2 胜。
它能看出简单的杀棋（例如 Qxf7#、底线杀），但中残局里仍会犯错，适合初学者到中级玩家练习。

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
  web/model/            转换后的权重（float16，约 17 MB）
  tools/az_reference.py numpy 参考实现
  tools/convert_weights.py
  tests/check_engine.mjs, tests/vectors.json
```

## 致谢与许可

模型结构、编码方式和权重来自 zeta36/chess-alpha-zero（MIT 许可，Copyright (c) 2017 Samuel Gravan）。
走法规则用 [chess.js](https://github.com/jhlywa/chess.js)，神经网络推理用 [TensorFlow.js](https://www.tensorflow.org/js)，均从 CDN 加载。
