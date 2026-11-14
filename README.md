# VCF 变异核对（FastAPI）

同一参考、同一样本下提供两项核查：

1. **双 VCF 归一化比对**（原有）：核对两套变异识别程序输出的单样本 VCF，按**归一化后的等位变异**
   列出两侧一致、仅单侧出现与剂量不符的结果，并保留每条结论的 VCF 原始行号；
2. **单样本冲突核查**（新增）：上传一份单样本 VCF，检查全部调用能否共同安排到**两条染色体拷贝**上，
   无解时给出一组导致无解的归一变异及其 VCF 原始行号。

## 单样本冲突核查规则

- **输入限制与原流程一致**：同样的文件大小上限（FASTA 1 MiB、VCF 256 KiB）、REF 校验、
  单条染色体纯 ACGT 参考、单样本 10 列 VCF、**非缺失未分相二倍体 GT**、左对齐归一与剂量聚合规则。
- **仅接受单碱基替换（SNV）与纯缺失（DEL）**；遇到插入或 MNV 等复杂变异直接判整次失败
  （返回 422，只显示错误，不展示局部结论）。
- **冲突定义**：两个变异若替换或删除了**同一个参考碱基**，就不能占用同一条染色体拷贝；
  缺失的**锚定碱基不计入**被删碱基（锚点上的 SNV 与该缺失不冲突）。
- **剂量语义**：剂量 1 的变异占一条拷贝，剂量 2（纯合）占满两条拷贝。
- **判定算法**：以归一变异为顶点、共享受影响参考碱基为边建冲突图——
  - 剂量 2 的顶点只要有相邻变异即无解（成对证据）；
  - 其余顶点均为剂量 1，冲突图**可二染色（无奇环）**才可安排；奇环即“任意两条可安排、
    三条合起来无解”的情形，例如同一参考碱基上的三个剂量 1 SNV。
- 页面展示**可行结论**，或一组导致无解的**归一变异、剂量、受影响碱基与原始行号 / 原始写法**；
  JSON 接口为 `POST /api/conflict`，页面表单为 `POST /conflict`。

## 双 VCF 比对规则

- **输入限制**
  - 单条染色体参考 FASTA，纯 A/C/G/T 序列，不超过 **1 MiB**；
  - 两份单样本 VCF（恰好 10 列、1 个样本），各不超过 **256 KiB**。
- **接受的变异**
  - A/C/G/T 构成的**单碱基替换**（SNV）；
  - 事件长度不超过 **20 bp** 的**纯插入 / 纯缺失**；
  - GT 必须是**非缺失、未分相的二倍体**（`0/1`、`1/1`，多 ALT 可为 `1/2`）。
  - 其他形态（MNV、符号化 ALT、`*`、分相 `|`、缺失 `.`、非二倍体等）一律判失败。
- **归一化**：对照参考序列对 indel 做左对齐（与 `bcftools norm` / `vt normalize` 等价），
  重复区域中位置写法不同的同一插入 / 缺失归为同一等位变异；SNV 位置不变。
- **剂量计算**：多 ALT 记录按 GT 分别计算每个 ALT 的 0/1/2 剂量，剂量为 0 的变异忽略；
  同一文件内同一归一变异的多条记录剂量累加，超过二倍体上限 2 即判定为**互相矛盾**。
- **整次失败**（返回 422，只显示错误、不展示任何局部结论）：
  - 任一记录 REF 与参考序列不符（含坐标越界、CHROM 不一致）；
  - 同一文件中同一归一变异存在互相矛盾的剂量；
  - 以及任何不满足上述输入规则的情形。
- **结果分组**（均按归一变异排序，并标注两侧原始行号与原始 `POS:REF>ALT` 写法）：
  两侧一致 / 仅 VCF A 或 B 出现 / 剂量不符。

## 目录结构

```
app/            FastAPI 应用与核心比对逻辑（bio.py 仅用标准库）
  main.py       路由与上传校验（/compare、/conflict 及对应 /api 接口）
  bio.py        FASTA/VCF 解析、左对齐归一化、剂量分组、拷贝冲突图判定
  templates/    上传页、双 VCF 结果页、冲突核查结果页（Jinja2）
tests/          单元 / 端到端测试
samples/        演示数据：demo.fa、双 VCF 比对样例与三份冲突核查样例
scripts/        演示 VCF 生成器（坐标由参考自动计算）
Dockerfile, docker-compose.yml, requirements.txt
```

## 用 Docker Compose 启动（推荐）

```bash
docker compose up --build
```

浏览器入口（启动后可直接在页面上操作两项核查）：

- 上传页面：<http://localhost:8000/>
- 交互式 API 文档（Swagger UI）：<http://localhost:8000/docs>

停止：`docker compose down`。容器以只读根文件系统、`no-new-privileges` 运行，
内置 `/` 健康检查。

## 本地运行（无 Docker 时）

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 试用

打开 <http://localhost:8000/>：

- **双 VCF 比对**：参考 `samples/demo.fa`，VCF A/B 分别选 `samples/caller_a.vcf`、
  `samples/caller_b.vcf`；
- **单样本冲突核查**：参考仍选 `samples/demo.fa`，待核查 VCF 可选
  - `samples/conflict_feasible.vcf`（可行：同碱基两个剂量 1 分占两拷贝、纯合 SNV、
    缺失与锚点 SNV 共存）；
  - `samples/conflict_odd_cycle.vcf`（无解：同一碱基三个剂量 1 SNV 的奇环，附一个无关变异）；
  - `samples/conflict_hom.vcf`（无解：剂量 2 SNV 与同碱基另一剂量 1 SNV）。

也可以走 JSON 接口：

```bash
# 双 VCF 比对
curl -F 'fasta=@samples/demo.fa' \
     -F 'vcf_a=@samples/caller_a.vcf' \
     -F 'vcf_b=@samples/caller_b.vcf' \
     http://localhost:8000/api/compare

# 单样本冲突核查
curl -F 'fasta=@samples/demo.fa' \
     -F 'vcf=@samples/conflict_odd_cycle.vcf' \
     http://localhost:8000/api/conflict
```

演示数据由生成器自动产出（坐标由参考计算，保证 REF 合法）：

```bash
python3 scripts/make_demo.py            # 双 VCF 样例
python3 -I scripts/make_conflict_demo.py  # 冲突核查三份样例
```

## 测试

```bash
# 核心逻辑（仅标准库）
python3 -I -m unittest discover -s tests -t . -p 'test_bio.py' -v

# 全部测试（Web 层需要先安装 requirements.txt）
python3 -I -m unittest discover -s tests -t . -v
```
