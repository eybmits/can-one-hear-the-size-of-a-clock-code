Can One Hear the Size of a Clock? -- source code and data
Markus Baumann and Claudia Linnhoff-Popien
Quantum Applications and Research Laboratory (QAR-Lab), Ludwig-Maximilians-Universität München
Contact: markus.baumann@campus.lmu.de

This archive accompanies the Letter and its Supplemental Material. It contains the LaTeX sources of both,
the figures, the scripts and data that compute and draw every figure, and the reproduction files cited in
Sec. VII of the Supplemental Material.

CONTENTS

main.tex            Letter (REVTeX 4.2, two-column); it ends with \input{sm}, so one build gives the Letter
                    followed by the Supplemental Material
sm.tex              Supplemental Material (pages S1-S12)
refs.bib, main.bbl  references; one list, printed after the Letter and cited by the same numbers in the SM
build.sh            full build: pdflatex, bibtex (twice), pdflatex (three times)
wordcount.py        APS length estimate of the Letter (body, captions, displayed equations, figures)
figs/               the four figures of the Letter as printed (vector PDF)
figsrc/             scripts that compute the data of Figs. 1-4 and draw them; stored data in figsrc/data/;
                    see figsrc/README.md (networks and parameters behind every panel: SM Sec. VII, table)
repro/              reproduction scripts and logs cited in the SM:
                      standalone_clocks.py   standalone clocks with graded entry rates, Eq. (S4), in
                                             60-digit arithmetic (SM Sec. I A)
                      sharing_resolution.py  distinguishability of the rhythms of nested cycles (SM Sec. III B)
                      finite_d.log           finite chains against standalone clocks, and the drive of
                                             shared paths (SM Secs. II C and IV F)
                      finite_d.py            prints finite_d.log
                      verify_certificates.py exact rational rates and inclusion discs of the two networks
                                             that break the per-rhythm drive bound (SM Sec. IV E)
                      models.py              network constructions used by sharing_resolution.py and
                                             finite_d.py
                      exactcert.py           exact-arithmetic helpers used by verify_certificates.py
data/               stored results cited in the SM: resolved_rhythms.json (optimized networks, SM Sec. III C),
                    resolved_windows.json (the 12-state search, SM Sec. III D), worked_example.json
                    (the hidden machine, SM Sec. VI C); the exact rates read by verify_certificates.py
                    (additive_law_witness_d27.json, additive_law_witness_stalled_d27.json,
                    c4_counterexample_d62.json, per_index_affinity_witness_d11.json,
                    near_equilibrium_witnesses.json)

REPRODUCE

  ./build.sh                                   # main.pdf: Letter (pages 1-6) and SM (S1-S12)
  python3 wordcount.py main.tex                # length estimate
  cd figsrc && python3 data.py                 # all figure data, about 10 min
  python3 fig1.py; python3 fig2.py; python3 fig3.py; python3 fig4.py
  cp out/fig*.pdf ../figs/ && cd ..
  python3 repro/standalone_clocks.py           # about 6 s
  python3 repro/sharing_resolution.py          # about 7 min; writes data/sharing_resolution.json
  python3 repro/finite_d.py                    # about 4 min; prints finite_d.log
  python3 repro/verify_certificates.py         # about 1 min

NAMES IN THE REPRODUCTION FILES

The reproduction scripts and finite_d.log keep the names of an earlier draft: theta chain = nested cycles;
independent rings, or rings glued at a hub = standalone clocks. In finite_d.log, "Eq. (2)" is the sharp form
of the capacity (Corollary S1), F~_Sigma is the total drive F_tot of Eq. (5), F_Sigma is the plain sum of
cycle affinities, and "Eq. (7)" is the earlier bound F_Sigma > pi N (Q - 1/Q).

REQUIREMENTS

pdfLaTeX and BibTeX with REVTeX 4.2 and the orcidlink package; Python 3 with numpy, scipy, mpmath, and
matplotlib 3.8-3.10 (3.11.0 drops minus signs from LaTeX labels in PDFs); for the figure scripts also LaTeX
with the type1cm and cm-super packages (text.usetex).
