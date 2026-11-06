"""Rough APS (PRL) length estimate for main.tex.

Counts: body words (text between \\maketitle and the AI-disclosure paragraph,
captions included, inline math = 1 word), displayed-equation rows
(16 words per row in a single column), and figures by the APS formula:
single column: 150/aspect + 20; two-column: 300/(0.5*aspect) + 40.
"""
import re, sys

src = open(sys.argv[1] if len(sys.argv) > 1 else "main.tex").read()
src = re.sub(r"(?<!\\)%.*", "", src)
body = src.split(r"\maketitle", 1)[1]
body = body.split(r"\bibliography", 1)[0]
# acknowledgments are excluded by APS; data availability is kept (conservative)
body = re.sub(r"\\runin\{Acknowledgments\.\}.*?(?=\\runin\{Data availability)", " ", body, flags=re.S)

# figures: pull out captions, remember type
figs = []
def grab_fig(m):
    env, content = m.group(1), m.group(2)
    cap = re.search(r"\\caption\{(.*)\}\s*\\label", content, re.S)
    figs.append((env, cap.group(1) if cap else ""))
    return " "
body = re.sub(r"\\begin\{(figure\*?)\}(.*?)\\end\{\1\}", grab_fig, body, flags=re.S)

# displayed equations: count rows
rows = 0
def grab_eq(m):
    global rows
    content = m.group(2)
    rows += content.count(r"\\") + 1
    return " EQ "
body = re.sub(r"\\begin\{(equation\*?|align\*?)\}(.*?)\\end\{\1\}", grab_eq, body, flags=re.S)

def words(t):
    t = re.sub(r"\$[^$]*\$", " M ", t)            # inline math -> 1 word
    t = re.sub(r"\\(cite|ref|eqref|label)\{[^}]*\}", " ", t)
    t = re.sub(r"\\SMsec\{[^}]*\}", " ", t)
    t = t.replace(r"\SM", " ")
    t = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", t)
    t = re.sub(r"[{}~\\]", " ", t)
    return len([w for w in t.split() if re.search(r"[A-Za-z0-9]", w) and w != "EQ"])

import subprocess
def _aspect(name):
    out = subprocess.run(["pdfinfo", f"figs/{name}.pdf"], capture_output=True, text=True).stdout
    w, h = map(float, re.search(r"Page size:\s+([\d.]+) x ([\d.]+)", out).groups())
    return w / h
aspect = {n: _aspect(n) for n in ("fig1", "fig2", "fig3", "fig4")}
body_words = words(body)
cap_words = sum(words(c) for _, c in figs)
fig_names = re.findall(r"includegraphics\[[^\]]*\]\{figs/([a-z0-9_]+)\.pdf\}", src)
fig_equiv = 0
for (env, _), name in zip(figs, fig_names):
    a = aspect[name]
    fig_equiv += (300 / (0.5 * a) + 40) if env == "figure*" else (150 / a + 20)
eq_words = 16 * rows
bib = open("refs.bib").read()
smnote = re.search(r"@misc\{SM,\s*note = \{(.*?)\}\s*\}", bib, re.S).group(1)
sm_words = len(smnote.split())
# the SM endnote sits in the reference list, which APS excludes; it is reported but not added
total = body_words + cap_words + fig_equiv + eq_words
print(f"SM endnote {sm_words} (reference list, not counted)")
print(f"body {body_words}  captions {cap_words}  eq rows {rows} (= {eq_words})  figures {fig_equiv:.0f}")
print(f"TOTAL ~ {total:.0f} word equivalents (PRL limit 3750)")

abs_txt = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", src, re.S).group(1).strip()
abs_txt = re.sub(r"\\looseness=-?\d+\s*", "", abs_txt)   # typesetting only, not text
def _render(m):
    t = m.group(1)
    for k, v in ((r"\pi", "π"), (r"\le", "≤"), (r"\ge", "≥"), (r"\sqrt", "√")):
        t = t.replace(k, v)
    return re.sub(r"[{}\\ ]", "", t)
abs_plain = re.sub(r"\s+", " ", re.sub(r"\$([^$]*)\$", _render, abs_txt))
print(f"abstract: {len(abs_plain)} characters incl. spaces, math as rendered (limit 600)")
print(abs_plain)
