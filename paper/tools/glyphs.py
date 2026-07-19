#!/usr/bin/env python3
"""Shared glyph pass: map prose Unicode (Greek, sub/superscripts, math operators) that
STIX Two Text lacks in text mode into inline $...$ LaTeX, WITHOUT touching characters that
are already inside $...$ math or `code` spans. Reads APDM_galactic.md, writes .build.md."""
import re, sys, os

GREEK = {'α':'\\alpha','β':'\\beta','γ':'\\gamma','δ':'\\delta','ε':'\\epsilon','ζ':'\\zeta',
 'η':'\\eta','θ':'\\theta','ι':'\\iota','κ':'\\kappa','λ':'\\lambda','μ':'\\mu','ν':'\\nu',
 'ξ':'\\xi','π':'\\pi','ρ':'\\rho','σ':'\\sigma','ς':'\\varsigma','τ':'\\tau','υ':'\\upsilon',
 'φ':'\\phi','ϕ':'\\phi','χ':'\\chi','ψ':'\\psi','ω':'\\omega','ϵ':'\\epsilon','ϑ':'\\vartheta',
 'Γ':'\\Gamma','Δ':'\\Delta','Θ':'\\Theta','Λ':'\\Lambda','Ξ':'\\Xi','Π':'\\Pi','Σ':'\\Sigma',
 'Φ':'\\Phi','Ψ':'\\Psi','Ω':'\\Omega','𝒟':'\\mathcal{D}','𝒬':'\\mathcal{Q}','𝒥':'\\mathcal{J}',
 'ℤ':'\\mathbb{Z}','ℓ':'\\ell','𝜋':'\\pi','𝜏':'\\tau','𝜇':'\\mu','𝜎':'\\sigma','𝜌':'\\rho'}
SUP = {'⁰':'0','¹':'1','²':'2','³':'3','⁴':'4','⁵':'5','⁶':'6','⁷':'7','⁸':'8','⁹':'9','⁺':'+','⁻':'-','ⁿ':'n'}
SUB = {'₀':'0','₁':'1','₂':'2','₃':'3','₄':'4','₅':'5','₆':'6','₇':'7','₈':'8','₉':'9','₊':'+','₋':'-'}
OPS = {'→':'\\to','←':'\\leftarrow','↔':'\\leftrightarrow','⟺':'\\iff','⟹':'\\implies','⇒':'\\Rightarrow',
 '⇔':'\\Leftrightarrow','≲':'\\lesssim','≳':'\\gtrsim','≈':'\\approx','≃':'\\simeq','≅':'\\cong',
 '≤':'\\leq','≥':'\\geq','≪':'\\ll','≫':'\\gg','≠':'\\neq','≡':'\\equiv','∼':'\\sim','∝':'\\propto',
 '⟨':'\\langle','⟩':'\\rangle','∈':'\\in','∉':'\\notin','⊂':'\\subset','⊃':'\\supset','∪':'\\cup',
 '∩':'\\cap','·':'\\cdot','∓':'\\mp','±':'\\pm','√':'\\surd','∞':'\\infty','∂':'\\partial',
 '∇':'\\nabla','∫':'\\int','∑':'\\sum','∏':'\\prod','⊗':'\\otimes','⊕':'\\oplus','⊙':'\\odot',
 '∀':'\\forall','∃':'\\exists','∅':'\\emptyset','∎':'\\blacksquare','†':'\\dagger','⋅':'\\cdot'}

def convert(text):
    # protect $$...$$, $...$, and `code`
    parts = re.split(r'(\$\$.*?\$\$|\$[^$\n]*\$|`[^`\n]*`)', text, flags=re.S)
    for i in range(0, len(parts), 2):
        s = parts[i]
        # combining tilde over a base letter (e.g. m̃) -> $\tilde{m}$
        s = re.sub(r'(\w)̃', lambda m: '$\\tilde{' + m.group(1) + '}$', s)
        # runs of sub/superscript digits -> single $_{...}$ / $^{...}$
        s = re.sub('([' + ''.join(SUP) + ']+)', lambda m: '$^{' + ''.join(SUP[c] for c in m.group(1)) + '}$', s)
        s = re.sub('([' + ''.join(SUB) + ']+)', lambda m: '$_{' + ''.join(SUB[c] for c in m.group(1)) + '}$', s)
        for k, v in GREEK.items(): s = s.replace(k, f'${v}$')
        for k, v in OPS.items():   s = s.replace(k, f'${v}$')
        s = s.replace('−', '-')
        parts[i] = s
    # code spans are protected from $-wrapping, but Unicode super/subscripts inside `\texttt`
    # choke under unicode-math (STIX Two Math) as stray superscript operators -> map to ASCII
    # caret/underscore, which \texttt renders literally (e.g. `T⁴`->`T^4`, `10⁻⁴`->`10^-4`).
    for i in range(1, len(parts), 2):
        if parts[i].startswith('`'):
            parts[i] = re.sub('([' + ''.join(SUP) + ']+)',
                              lambda m: '^' + ''.join(SUP[c] for c in m.group(1)), parts[i])
            parts[i] = re.sub('([' + ''.join(SUB) + ']+)',
                              lambda m: '_' + ''.join(SUB[c] for c in m.group(1)), parts[i])
    t = ''.join(parts)
    # merge adjacent $...$$...$ (e.g. Greek then subscript) into one math span
    t = re.sub(r'\$([^$]*)\$\$([^$]*)\$', r'$\1\2$', t)
    t = re.sub(r'\$([^$]*)\$\$([^$]*)\$', r'$\1\2$', t)   # twice for triples
    # final sweep: unicode-math chars the agent wrote INSIDE $...$ spans (which we protected)
    # -> their LaTeX commands; valid in math, and text-mode ones are already $-wrapped above.
    for k, v in {'√': '\\surd ', '𝜋': '\\pi ', '𝜏': '\\tau ', '𝜇': '\\mu ', '𝜎': '\\sigma ',
                 '𝜌': '\\rho ', '𝒥': '\\mathcal{J} ', '𝒟': '\\mathcal{D} '}.items():
        t = t.replace(k, v)
    return t

if __name__ == "__main__":
    inp = sys.argv[1] if len(sys.argv) > 1 else "APDM_galactic.md"
    out = sys.argv[2] if len(sys.argv) > 2 else ".build.md"
    open(out, "w", encoding="utf-8").write(convert(open(inp, encoding="utf-8").read()))
