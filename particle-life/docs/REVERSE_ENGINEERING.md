# 逆向工程记录 · What the platform actually computes

本文件记录如何从 `https://sandbox-science.com/particle-life` 的线上产物中，
把它真正执行的物理内核**逐字**取出来，以及取到了什么。
所有下游的理论（`THEORY.md`）和复现（`web/`, `theory/`）都以此为准，
而不是以"我猜它大概是这样"为准。

---

## 1. 提取路径

站点是一个 Nuxt (Vue 3) SPA，页面 HTML 里只有一个 loading 壳，因此：

```
GET /particle-life                       -> 只有 <div id="__nuxt"> 和 entry chunk
GET /_nuxt/<entry>.js                    -> 路由 -> 页面 chunk 的映射
GET /_nuxt/<page-chunk>.js               -> pinia store：所有默认参数
   __vite__mapDeps 里出现 ParticleLifeGpu / ParticleLifeGpu3D
GET /_nuxt/<ParticleLifeGpu>.js          -> 16 个 @compute WGSL 内核（明文模板字符串）
```

关键点：**模拟内核没有被编译掉**——WGSL 以模板字符串形式原样躺在 bundle 里，
所以拿到的是执行代码本身，不是复刻品。

---

## 2. 力的定律（逐字，来自 `forceCompute` 内核）

```wgsl
let dist = sqrt(dx * dx + dy * dy);
let index = typeA * options.numTypes + typeB;
let params = get_interaction(index, options.numTypes);   // (rule, minR, maxR)
let maxR = params.z;
if (dist > 0.0 && dist < maxR) {
    let rule = params.x;
    let minR = params.y;
    var force = 0.0;
    if (dist < minR) {
        force = (dist / minR - 1.0) * options.repel;
    } else {
        let mid = (minR + maxR) / 2.0;
        let slope = rule / (mid - minR);
        force = -(slope * abs(dist - mid)) + rule;
    }
    if (force != 0.0) {
        velocitySum.x += dx * (force / dist);
        velocitySum.y += dy * (force / dist);
    }
}
...
let forceFactor = options.forceFactor * deltaTime * 60.0;
particle.vx += velocitySum.x * forceFactor;
particle.vy += velocitySum.y * forceFactor;
```

推进内核 `particleAdvance`：

```wgsl
particle.vx *= dt.friction;          // friction = pow(1 - frictionFactor, deltaTime*60)  (CPU 端预算)
particle.vy *= dt.friction;
particle.x += particle.vx * dt.deltaTime;
particle.y += particle.vy * dt.deltaTime;
// isWallRepel: 反射；isWallWrap: 环面回绕；两者皆假 = 开放边界
```

于是精确的离散映射是

```
v <- ( v + 60 κ Δt · Σ_j f_{s_i s_j}(r_ij) r̂_ij ) · (1-λ)^{60Δt}
x <- x + v Δt
```

三条容易被忽略、但对理论至关重要的事实：

1. **力矩阵是有序对矩阵**：索引 `typeA * numTypes + typeB`。
   `A[a][b] ≠ A[b][a]` 完全合法 —— 牛顿第三定律在这里是**可选项**。
2. **min/max 半径同样是完整的 N×N 矩阵**，也可以不对称
   （UI 上就是 `Forces | Min. Radius | Max. Radius` 三个页签）。
3. 力直接加到速度上，质量恒为 1；`repel` 是唯一的全局标量（核心排斥强度）。

矩阵在 GPU 上被打包进一个 `u32`：

```wgsl
let rule = (f32((word >> 0u)  & 0xFFu) / 255.0) * 2.0 - 1.0;   // 8 bit,  [-1, 1]
let minR =  f32((word >> 8u)  & 0xFFu);                        // 8 bit,  [0, 255]
let maxR =  f32((word >> 16u) & 0xFFFFu);                      // 16 bit
```

即力系数只有 **8 bit 量化**（步长 2/255 ≈ 0.0078）。

---

## 3. 邻居搜索

三段式 GPU 空间哈希，格子边长 = 全局 `maxRadius`：

| 内核 | 作用 |
|---|---|
| `clearBinSize` / `fillBinSize` | 原子计数每个 bin 的粒子数 |
| `prefixSumStep` | 并行前缀和 -> bin 偏移 |
| `sortParticles` | 按 bin 重排粒子数组（保证访存局部性） |

非回绕模式下用 `extendedGrid* + gridOffset*` 把网格向外扩一圈，
这样开放边界不需要特殊分支。

---

## 4. 默认参数（来自 pinia store）

| 参数 | GPU 引擎 | CPU 引擎 |
|---|---|---|
| `numParticles` | 64 000 | 6 000 |
| `numColors` (物种数 S) | 7 | 7 |
| `particleSize` | 2 | 8 |
| `forceFactor` κ | 1 | 1 |
| `frictionFactor` λ | 0.3 | 0.3 |
| `repel` R | 1 | 1 |
| `minRadiusRange` | [12, 24] | [30, 60] |
| `maxRadiusRange` | [32, 64] | [90, 150] |
| `manualDeltaTime` | 0.0166 (=1/60) | — |
| 墙 | 开放 (wrap=repel=false) | 反射 |

λ = 0.3、Δt = 1/60 ⇒ 每帧速度乘 0.7，连续时间阻尼率 γ = −60 ln 0.7 = **21.4 s⁻¹**。

---

## 5. 预设生成器

站点带 31 个矩阵生成器、37 个配色、28 个二维分布（另有 3D 分布）。
`theory/plife/matrices.py` 与 `web/js/model.js` 里 1:1 移植了其中
结构上最有代表性的 19 个矩阵生成器，例如：

```js
// id 6 - Rock–Paper–Scissors
r[l][s] = (s===l) ? -0.1 : (s===(l+1)%e) ? 0.9 : (s===(l+e-1)%e) ? -0.7 : 0

// id 16 - Wavefield  (纯反对称：sin 是奇函数)
o[r][l] = (r===l) ? -0.05 : Math.sin(2*Math.PI*(l-r)/e) * 0.9

// id 2 - Snake
t[n][o] = (o===n) ? 1 : (o===(n+1)%e) ? 0.2 : 0
```

这些生成器的**图结构**（谁连谁、连的方向对不对称）正是形态学的控制量，
见 `THEORY.md` 第 8 节。

---

## 6. 复现与原版的差异（诚实清单）

| 项 | 本复现 | 原站 |
|---|---|---|
| 计算后端 | CPU（typed array + 均匀网格），JS | WebGPU compute（另有 CPU/3D 引擎） |
| 规模 | 约 10⁴ 粒子 @ 60 fps | 2×10⁵ 粒子 |
| 力系数精度 | float64 / float32 | 8 bit 量化 |
| 3D 模式 | 无 | 有 |
| 相机 / GIF 录制 / 笔刷类型 / 预设云端 | 仅基础笔刷与导入导出 | 完整 |
| 物理内核 | **相同** | — |
| 积分器与更新顺序 | **相同** | — |
| 边界模式 | 相同（wrap / bounce / open） | — |

即：**物理是等价的，工程规模不是**。
本复现额外提供原站没有的东西——右侧 "Physics X-ray" 面板，
把连续介质线性稳定性理论实时算在浏览器里，与模拟测量并排显示。
浏览器里的这套计算已与 `theory/plife` 的 Python 实现逐点对齐到 1.2×10⁻⁴ 相对误差
（见 `docs/VALIDATION.md`）。
