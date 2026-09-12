from pathlib import Path
import sys
from pypdf import PdfReader, PdfWriter
from pypdf.generic import NameObject, DecodedStreamObject

def portable(src,dst):
    mapping=Path('/Users/zhli/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/share/poppler/cidToUnicode/Adobe-GB1').read_text().splitlines()
    entries=[]
    for cid,v in enumerate(mapping):
        cp=int(v,16)
        if cp:entries.append(f'<{cid:04X}> <{chr(cp).encode("utf-16-be").hex().upper()}>')
    parts=['/CIDInit /ProcSet findresource begin','12 dict begin','begincmap','/CIDSystemInfo << /Registry (Adobe) /Ordering (UCS) /Supplement 0 >> def','/CMapName /Adobe-GB1-ToUnicode def','/CMapType 2 def','1 begincodespacerange','<0000> <FFFF>','endcodespacerange']
    for i in range(0,len(entries),100):
        batch=entries[i:i+100];parts.extend([f'{len(batch)} beginbfchar',*batch,'endbfchar'])
    parts+=['endcmap','CMapName currentdict /CMap defineresource pop','end','end']
    w=PdfWriter(clone_from=src);s=DecodedStreamObject();s.set_data(('\n'.join(parts)+'\n').encode('ascii'));ref=w._add_object(s.flate_encode())
    count=0
    for o in w._objects:
        if not hasattr(o,'get') or o.get('/Subtype')!='/Type0' or '/ToUnicode' in o:continue
        ds=o.get('/DescendantFonts',[])
        if not ds:continue
        ci=ds[0].get_object().get('/CIDSystemInfo',{})
        if ci.get('/Registry')=='Adobe' and ci.get('/Ordering')=='GB1':
            assert o.get('/Encoding')=='/Identity-H';o[NameObject('/ToUnicode')]=ref;count+=1
    w.write(dst)
    a,b=PdfReader(src),PdfReader(dst)
    assert len(a.pages)==len(b.pages)
    for x,y in zip(a.pages,b.pages):assert x.get_contents().get_data()==y.get_contents().get_data()
    print('Portable fonts',count,'pages',len(b.pages),dst)

if __name__=='__main__':portable(Path(sys.argv[1]),Path(sys.argv[2]))
