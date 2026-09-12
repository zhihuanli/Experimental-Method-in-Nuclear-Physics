# L1—L4 可编辑讲义源文件

本目录保存四册独立中文讲义的正文、排版配置、图件及导出辅助文件。已发布的四份 PDF 位于上一级目录，本次归档不改变讲义内容。

## 正文与排版文件

| 讲义 | Markdown 正文 | LaTeX 排版文件 | 排版配置 |
|---|---|---|---|
| L1 | [L1.md](independent/L1.md) | [L1.tex](independent/L1.tex) | [L1.yaml](independent/L1.yaml) |
| L2 | [L2.md](independent/L2.md) | [L2.tex](independent/L2.tex) | [L2.yaml](independent/L2.yaml) |
| L3 | [L3.md](independent/L3.md) | [L3.tex](independent/L3.tex) | [L3.yaml](independent/L3.yaml) |
| L4 | [L4.md](independent/L4.md) | [L4.tex](independent/L4.tex) | [L4.yaml](independent/L4.yaml) |

以 Markdown 为正文主文件：标题、中文段落、公式、表格及图注均可编辑。LaTeX 文件是由对应 Markdown 生成的版本；后续修改后应重新生成，避免同时维护两套不同的正文。

`l1_curated_figures` 至 `l4_curated_figures` 保存配套 PDF 图件，包括矢量线条及嵌入的照片、测量图。图像内部文字不属于 Markdown 正文。请保持这些文件夹与 `independent` 为同级目录，正文中的相对路径依赖这一结构。

## 重新导出

原导出环境使用 Pandoc、XeLaTeX、Python 与 `pypdf`，中文字体为 Fandol。辅助文件保留原工作环境的配置：

- `independent_compile.py`：把 Markdown 转成 LaTeX，连续编译三次以更新目录与交叉引用，再生成补齐中文字符映射的 PDF。
- `independent_pdf.py`：利用 Poppler 的 Adobe-GB1 字符表补齐缺失的 Unicode 映射，不改变页面绘制内容。

在本目录执行以下命令可导出全部四册；把末尾数字改为 `1` 可仅导出 L1：

```sh
python3 independent_compile.py 1 2 3 4
```

脚本及 YAML 保留了原机器上的 Pandoc、XeLaTeX、字体和 Poppler 字符表的绝对路径。在其他电脑上使用前，请改为该电脑的实际安装路径，并确认执行脚本的 Python 已有 `pypdf`。也可以在 `independent` 目录中手动使用 Pandoc 和 XeLaTeX 导出。

生成的 `L1-portable.pdf` 等文件位于 `independent` 目录；脚本不会自动覆盖上一级目录的正式讲义。重新导出后需检查中文、公式、插图边界、目录和交叉引用，再更新正式版本。

本目录不包含旧稿、原始幻灯片、临时页面预览、检查日志或逐页映射表。
