"""Keep legacy learning URLs as entrances to their already-integrated lessons."""
from pathlib import Path

root = Path(__file__).resolve().parents[1]
pages = {
    'binned ML-decay.html': ('likelihood_decay.html', 'Exponential decay: unbinned and binned likelihood'),
    'unbinned ML-decay.html': ('likelihood_decay.html', 'Exponential decay: unbinned and binned likelihood'),
    'binned ML-gaus+bkg.html': ('likelihood_signal_background.html', 'Signal plus background: unbinned, binned, and extended likelihood'),
    'unbinned ML-gaus+bkg.html': ('likelihood_signal_background.html', 'Signal plus background: unbinned, binned, and extended likelihood'),
}
for name, (target, title) in pages.items():
    assert (root/'chapt2'/target).is_file()
    page = f'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title><link rel="canonical" href="{target}">
<style>body{{max-width:850px;margin:50px auto;padding:0 24px;font:17px/1.8 system-ui,sans-serif}}a{{color:#145787}}</style>
</head><body><nav><a href="../reference.html">学习资料</a></nav>
<h1>{title}</h1><p>原 binned 与 unbinned 示例已整合，统一比较模型、参数估计与结果。</p>
<p><a href="{target}">阅读整合后的实例（可切换 PyROOT / ROOT C++）</a></p>
</body></html>'''
    (root/'chapt2'/name).write_text(page)

central = root/'chapt2/2.centerlimit.html'
source = central.read_text()
source = source.replace('中心极值定理', '中心极限定理').replace('Guassian', 'Gaussian').replace('Guass', 'Gaussian')
central.write_text(source)
