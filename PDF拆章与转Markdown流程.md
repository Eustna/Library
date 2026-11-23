# PDF 拆章与转 Markdown 流程

记录把整本 PDF（含扫描版）拆成按章 PDF、再转成 Markdown 的完整做法。适用于教材、讲义等需要按章阅读或导入知识库的书。

## 工具

- Python 3.13（anaconda 环境）
- `pypdf`：读书签、拆页、写分章 PDF
- `pymupdf`：把 PDF 页面渲染成 PNG
- Playwright（Python 包）：操作 DeepSeek 网页
- Chrome：Playwright 使用 `channel="chrome"` 启动
- 可选：DeepSeek 网页账号（网页版能直接读 PDF 并输出 Markdown）

装包：

```powershell
pip install pypdf pymupdf playwright
```

---

## 一、拆章

### 1. 先看书签结构

```python
import pypdf

r = pypdf.PdfReader('book.pdf')
def walk(items, depth=0):
    for it in items:
        if isinstance(it, list):
            walk(it, depth + 1)
        else:
            print('  ' * depth + it.title)

walk(r.outline)
print('pages', len(r.pages))
```

一级书签通常是章，二级书签是 §x.y 节，三级是 §x.y.z。拆章按一级书签；没有书签时按目录页人工确定页码范围。

### 2. 按页码拆

页码一律用原 PDF 的 1-based 页码，起止页都包含。

```python
import pypdf, os

src = r'C:\Users\17869\Desktop\book.pdf'
out = r'C:\Users\17869\Desktop\book-分章'
os.makedirs(out, exist_ok=True)

chapters = [
    (1, '绪言', 6, 11),
    (2, '聚变能', 12, 19),
    # (章号, 标题, 起始页, 结束页)
]

r = pypdf.PdfReader(src)
for num, title, a, b in chapters:
    w = pypdf.PdfWriter()
    for p in range(a - 1, b):
        w.add_page(r.pages[p])
    name = f'第{num:02d}章 {title}.pdf'
    with open(os.path.join(out, name), 'wb') as fh:
        w.write(fh)
    print(name, b - a + 1, '页')
```

前置、附录单独拆：

```python
(0,  '前置（封面·书名·版权·前言）', 1, 5),
(13, '附录（符号索引·专有名词索引）', 340, 356),
```

### 3. 目录清单

在分章目录里写 `目录清单.txt`，记录：

- 原文件名、总页数
- 拆分规则（按一级书签，页码为原 PDF 1-based）
- 每个分章文件对应的页码范围和页数

命名统一为 `第01章 标题.pdf`、`第02章 标题.pdf`，两位数补零，排序稳定。

---

## 二、转 Markdown

两种方案。优先用方案 A（识别质量好），方案 A 被限流时用方案 B（本地视觉）。

### 方案 A：DeepSeek 网页（老办法）

原理：把分章 PDF 上传到 DeepSeek 网页，让它输出 Markdown 代码块，再点“下载”得到 md 文件。

#### 2A.1 建立浏览器上下文

```python
from playwright.sync_api import sync_playwright

PROFILE = r'C:\Users\17869\AppData\Local\ms-playwright-mcp\mcp-chrome-f2e1067'

with sync_playwright() as p:
    ctx = p.chromium.launch_persistent_context(
        PROFILE,
        channel='chrome',
        headless=False,
        args=['--start-maximized'],
        no_viewport=True,
    )
    page = ctx.pages[0] if ctx.pages else ctx.new_page()
    page.goto('https://chat.deepseek.com', wait_until='domcontentloaded')
```

`PROFILE` 用已有的持久化目录，能保留登录态。没有登录态时先手动登录一次。

#### 2A.2 上传、发提示词、等待、下载

```python
import time

pdf = r'...\第03章 磁约束聚变原理.pdf'
prompt = (
    '请把上传的 PDF 完整转换为 Markdown。'
    '只输出一个 markdown 代码块，第一行写 <!-- 文件名：第03章 磁约束聚变原理.md -->。'
    '只转换正文，跳过页眉页脚页码。'
    '公式用 LaTeX：行内 \\(...\\)，独立 \\[...\\]。'
    '图注用 [图x.x 占位] 标出。'
    '不要总结、不要删减、不要输出解释，不要写“未完”。'
)

page.get_by_role('textbox', name='给 DeepSeek 发送消息')
page.locator('input[type=file]').last.set_input_files(pdf)
time.sleep(4)
box = page.get_by_role('textbox', name='给 DeepSeek 发送消息')
box.fill(prompt)
box.press('Enter')

# 等代码块出现，再等长度稳定
page.wait_for_selector('.md-code-block', timeout=120000)
time.sleep(15)
prev = -1
for _ in range(40):
    cur = page.evaluate(
        "() => {const b=document.querySelectorAll('.md-code-block');"
        "return b.length ? b[b.length-1].innerText.length : 0}"
    )
    if cur > 300 and cur == prev:
        break
    prev = cur
    time.sleep(6)

block = page.locator('.md-code-block').last
with page.expect_download(timeout=120000) as dl:
    block.get_by_role('button', name='下载').click()
dl.value.save_as(r'...\第03章 磁约束聚变原理.md')
```

下载文件第一行通常有 `<!-- 文件名：... -->`，合并前去掉。

#### 2A.3 长章节分块

DeepSeek 单次输出约 50 KB，超过会截断。经验阈值：

- 不超过 24 页：整章一次转
- 超过 24 页：按 15 页一块拆分，逐块转，最后按顺序合并

```python
import pypdf, os

r = pypdf.PdfReader(chapter_pdf)
pages = len(r.pages)
a = 1
while a <= pages:
    b = min(a + 14, pages)          # 15 页一块
    w = pypdf.PdfWriter()
    for i in range(a - 1, b):
        w.add_page(r.pages[i])
    with open(f'chunk_p{a}-{b}.pdf', 'wb') as fh:
        w.write(fh)
    a = b + 1
```

每块的提示词写“这是第 x/n 部分”，文件名用 `第03章_p1-15.md` 之类，最后按页码顺序拼接。

#### 2A.4 截断续写

输出到一半停止时，取最后一句做锚点，让 DeepSeek 接着写：

```text
上一个输出在 §9.4.3 中途被截断了，最后停在“也是 Krus”。
请从“kal 坐标……”继续，把本章剩余正文完整输出为另一个可下载的 .md 文件：
第一行写 <!-- 文件名：第09章 施瓦西时空_续.md -->。
不要重复前文，不要总结、不要删减，不要写“未完”。
```

#### 2A.5 限流处理

出现“服务器繁忙”“请删除异常文件再发送”时，是账号侧上传服务受限，与文件无关。处理：

- 等待 10–30 分钟重试
- 失败重试间隔逐步拉长：300 s、600 s、900 s
- 长时间不恢复，转方案 B

### 方案 B：本地视觉逐页

原理：把 PDF 每页渲染成 PNG，用支持视觉的模型逐页读图并转写，最后合并。

#### 2B.1 渲染页面

```python
import pymupdf, os

doc = pymupdf.open('第09章 托卡马克实验研究.pdf')
out = r'...\pages\第09章'
os.makedirs(out, exist_ok=True)
for i, page in enumerate(doc):
    pix = page.get_pixmap(dpi=150)
    pix.save(os.path.join(out, f'第09章_p{i+1:03d}.png'))
```

150 DPI 足够清晰；扫描版不要用 `pypdf` 提取文本，因为页面没有文本层。

#### 2B.2 逐页转写

每页单独写一个 md，避免一次输出过长：

```text
用 Read 工具读取第09章_p001.png，转写为 Markdown：
- 标题用 #/##/###，保留原书编号
- 正文完整，不总结、不删减
- 公式用 LaTeX：行内 \(...\)，独立 \[...\]，保留公式编号
- 图注写 [图x.x 占位]
- 去掉页眉页脚页码
- 无法确认的字保留原样，不要编造
输出到 第09章_p001.md
```

一批 8 页交给一个子代理，页数多时并行多组。

#### 2B.3 合并

```python
import io, os, re

d = r'...\vision\第09章'
pages = sorted(
    (int(m.group(1)), os.path.join(d, n))
    for n in os.listdir(d)
    if (m := re.match(r'第09章_p(\d+)\.md$', n))
)

parts = []
for _, p in pages:
    t = io.open(p, encoding='utf-8').read()
    t = re.sub(r'^\s*<!--\s*文件名[^\n]*-->\s*\n?', '', t)
    parts.append(t.strip('\n'))

io.open(r'...\第09章 托卡马克实验研究.md', 'w', encoding='utf-8').write('\n\n'.join(parts) + '\n')
```

---

## 三、质量检查

### 1. 书签覆盖

用 PDF 最后一两个书签核对 md 末尾是否包含，防止整章截断。

```python
# 比较最后几个书签标题与 md 文本
```

书签名和 md 标题常有措辞差异（“主要特性”对“主要特征”），模糊匹配即可。

### 2. 格式检查

- 页眉页脚页码：不应有单独成行的纯数字、书名行
- 公式：`\(` 与 `\)` 数量相等，`\[` 与 `\]` 数量相等
- 图注：保留 `[图x.x 占位]`
- 空文件、缺页：按页号连续检查

### 3. 常见 OCR 问题

- 公式编号缺失或重号：保留缺号，不臆造
- 上下标错位、符号串行：按上下文最小修正，拿不准的保留原样并标注
- 人名、装置名拼写：按原书保留，不擅自“更正”
- 扫描版没有文本层：一律走视觉方案

---

## 四、图片补入知识库（可选）

如果图片已按 `章.图号.png` 命名（如 `9.22.png` 对应图 9.22）：

```powershell
# 1. 复制图片
New-Item -ItemType Directory -Path "references\物理\等离子体物理\assets" -Force
Copy-Item "pic\*.png" "references\物理\等离子体物理\assets\"

# 2. 替换占位符
# [图9.22 占位] -> ![图9.22](assets/9.22.png)
```

正则：

```text
\[图\s*([0-9]+\.[0-9]+)\s*占位\]
```

替换为：

```text
![图$1](assets/$1.png)
```

替换后检查所有 `assets/*.png` 链接都存在，且没有残留占位符。

---

## 五、本次实际脚本

参考脚本在：

```text
C:\Users\17869\AppData\Local\Temp\opencode\
├─ split_pdfs.py       按页码拆章
├─ shi_md.py           DeepSeek 网页转 md（上传/等待/下载/分块/重试）
├─ render_pages.py     渲染页面 PNG
├─ merge_vision.py     视觉转写结果合并
└─ ds_probe.py         诊断 DeepSeek 页面状态
```

这些是一次性脚本，路径和章节表按书修改后复用。
