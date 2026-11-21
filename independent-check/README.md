# 检验材料包

- `检验报告.md`：完整检验报告
- `scripts/`：独立验算脚本（Python 3 + gcc，无第三方依赖）
  - `hbp.py`：通用 H-自举渗流模拟器
  - `sim.c` / `simr.c`：位集模拟器（`simr` 可限定轮数）。编译：`gcc -O2 -o sim sim.c`、`gcc -O2 -o simr simr.c`
  - `t_theta.py`、`t_theta2.py`：θ 图自闭包分类（`t_theta2.py` 需要先编译 `sim`）
  - `t_typed.py`：12 顶点带类型规则，E→E⁺ 与 ⌊lcm/2⌋
  - `t_rule.py`：重建 166,464 顶点规则，检查阶数、边数、2-连通、树宽
  - `t_host.py`：实际宿主图上的局部化不等式、度数，以及全规模提升见证（依赖 `t_rule.py`）
  - `t_lemma22.py`：Lemma II.2 的 2h 轮上界随机测试（需要先编译 `simr`）
  - 运行方式：在 `scripts/` 目录内执行 `python3 <脚本名>`
- `logs/`：`lake build` 与内核重放审计的完整日志，以及 θ 闭包模拟的输出
- 所检验的预印本 v1.0 见仓库 README 中的发布链接（此处不重复存放 PDF）
