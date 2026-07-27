#!/usr/bin/env bash
# Build the MNRAS PDF (mnras.cls, one-column, author-year natbib) for Paper VIII
# from APDM_galactic.md. Pipeline identical to build_oja.sh; only the template,
# keyword separator and bibliography style differ. MNRAS is author-year, so the
# \citet/\citep distinction is preserved (no normalisation).
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
  --standalone --shift-heading-level-by=-1 --natbib --template=tools/template_mnras.tex
perl -0pi -e 's#(\\begin\{keywords\}\s*\n)([^\\]*?)(\n\s*\\end\{keywords\})#my ($a,$k,$z)=($1,$2,$3); $k =~ s{[;,]\s*}{ -- }g; "$a$k$z"#se' tex/$BASE.tex
# make long monospace paths breakable (prevents right-margin overflow)
perl -0pi -e 's{\\texttt\{([^{}]*)\}}{"\\texttt{".($1=~s#([/._])#$1\\allowbreak #gr)."}"}ge' tex/$BASE.tex
perl -0pi -e 's/\\\[/\\begin{equation}/g; s/\\\]/\\end{equation}/g' tex/$BASE.tex
( cd tex && \
  xelatex -interaction=nonstopmode $BASE.tex >$BASE.build.log 2>&1 ; \
  BIBINPUTS="..:$BIBINPUTS" bibtex $BASE      >>$BASE.build.log 2>&1 ; \
  xelatex -interaction=nonstopmode $BASE.tex >>$BASE.build.log 2>&1 ; \
  xelatex -interaction=nonstopmode $BASE.tex >>$BASE.build.log 2>&1 ) || true
if [ -f tex/$BASE.pdf ]; then
  cp tex/$BASE.pdf $BASE.pdf
  echo "built paper/$BASE.pdf  (MNRAS)"
  grep -c "^!" tex/$BASE.build.log | awk '{print $1" LaTeX errors (see tex/'$BASE'.build.log)"}'
  grep -c "Citation.*undefined" tex/$BASE.build.log 2>/dev/null | awk '{print $1" undefined citations"}'
else
  echo "BUILD FAILED — see tex/$BASE.build.log"; grep -A2 '^!' tex/$BASE.build.log | head -20
fi
