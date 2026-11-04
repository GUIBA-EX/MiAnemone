"""Apply Avenir Next to all editable scientific SVG figures, preserving data geometry."""
from pathlib import Path
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
ET.register_namespace('','http://www.w3.org/2000/svg')
count=0
for folder in [ROOT/'figures',ROOT/'publications/figures']:
 for path in folder.rglob('*.svg'):
  tree=ET.parse(path)
  for el in tree.iter():
   if el.tag.endswith('text') or el.tag.endswith('tspan'):
    el.set('font-family','Avenir Next');el.set('style',el.get('style','')+';font-family: Avenir Next;')
   elif el.get('font-family') is not None:el.set('font-family','Avenir Next')
  tree.write(path,encoding='utf-8',xml_declaration=True)
  content=path.read_text().replace('′',chr(39))
  if path.name.endswith('_binding.svg'):
   for a,b in {'#13aaba':'#159fa8','#288064':'#3d7993','#c88312':'#ce962d','#2779ef':'#4779b8','#922bde':'#9063a8','#697885':'#82888c','#f7fbff':'#ffffff'}.items():content=content.replace(a,b)
  path.write_text(content);count+=1
print(count,'scientific SVGs standardized')
