#!/usr/bin/env bash
# Build the Open Journal of Astrophysics PDF (openjournal.cls) for Paper VIII from APDM_galactic.md.
set -e
cd "$(dirname "$0")"
mkdir -p tex
BASE=APDM_galactic
python3 tools/glyphs.py $BASE.md .g.md
python3 tools/makedoc.py .g.md .build.md
trap 'rm -f .build.md .g.md' EXIT
python3 - <<'PY'
import re
t=open('.build.md',encoding='utf-8').read()
# split off the YAML front matter: math-protection must NOT touch it (a \( inside a
# double-quoted YAML title is an invalid escape). pandoc handles metadata math itself.
mfm=re.match(r'^---\n.*?\n---\n', t, re.S)
head = mfm.group(0) if mfm else ''
body = t[len(head):]
body=re.sub(r'\n##\s+References\s*\n+:::\s*\{#refs\}\s*\n:::\s*\n','\n',body)
body=re.sub(r'(?<!\\)(?<!\$)\$(?!\$)((?:\\.|[^$\\\n]|\\\n)+?)\$(?!\$)',
         lambda m: '`\\('+m.group(1)+'\\)`{=latex}', body)
open('.build.md','w',encoding='utf-8').write(head+body)
PY
pandoc -f markdown-superscript-subscript .build.md -o tex/$BASE.tex \
  --standalone --shift-heading-level-by=-1 --natbib --template=tools/template_oja.tex
perl -0pi -e 's/\\\[/\\begin{equation}/g; s/\\\]/\\end{equation}/g' tex/$BASE.tex
( cd tex && \
  xelatex -interaction=nonstopmode $BASE.tex >$BASE.build.log 2>&1 ; \
  BIBINPUTS="..:$BIBINPUTS" bibtex $BASE      >>$BASE.build.log 2>&1 ; \
  xelatex -interaction=nonstopmode $BASE.tex >>$BASE.build.log 2>&1 ; \
  xelatex -interaction=nonstopmode $BASE.tex >>$BASE.build.log 2>&1 ) || true
[ -f tex/$BASE.pdf ] && cp tex/$BASE.pdf $BASE.pdf && echo "built paper/$BASE.pdf" || echo "FAILED"
