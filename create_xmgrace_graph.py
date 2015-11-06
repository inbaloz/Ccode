template = open('xmgrace_template.agr','rb').read()
import sys
coords = open(sys.argv[1],'rb').read()

out_text = template.format(coords=coords)
out_filename = sys.argv[1]+'.agr'
open(out_filename,'wb').write(out_text)
