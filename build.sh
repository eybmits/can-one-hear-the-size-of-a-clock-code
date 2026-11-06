#!/bin/bash
# full build: Letter + Supplemental Material (one reference list, in the Letter)
cd "$(dirname "$0")"
pdflatex -interaction=nonstopmode main.tex > /dev/null 2>&1
bibtex main > bib.log 2>&1
pdflatex -interaction=nonstopmode main.tex > /dev/null 2>&1
bibtex main > bib.log 2>&1          # second pass: citation order inside the SM endnote
pdflatex -interaction=nonstopmode main.tex > /dev/null 2>&1
pdflatex -interaction=nonstopmode main.tex > /dev/null 2>&1
pdflatex -interaction=nonstopmode main.tex > /dev/null 2>&1
grep -E "^!|Warning|Overfull|Underfull" main.log | grep -v "nameref\|No type size" | sort | uniq -c
pdfinfo main.pdf | grep Pages
