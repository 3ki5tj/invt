REM MS-DOS version of `make zip`

SET prj=invt
SET /P ver=Please enter the version number, like 0.2:
SET name="%prj%V%ver%"
SET namea="%prj%_annotatedV%ver%"

mkdir "old/V%ver%"
rsync -avz %prj%.tex "old/V%ver%/%namea%.tex"
rsync -avz Makefile "old/V%ver%/"
rsync -avz README.md "old/V%ver%/"
rsync -avz deannote.py "old/V%ver%/"
rsync -avz rmbiburl.py "old/V%ver%/"
python rmbiburl.py -i simul.bib -o "old/V%ver%/simul.bib"
mkdir "old/%ver%/fig"
cp "fig/*" "old/%ver%/fig/"

cd "old/V%ver%"
python deannote.py -v -c -i %namea%.tex -o %name%.tex

pdflatex %name%
bibtex   %name%
pdflatex %name%
pdflatex %name%

pdflatex %namea%
bibtex   %namea%
pdflatex %namea%
pdflatex %namea%

7z a "../../%prj%docV%ver%.zip" *.tex *.pdf simul.bib *.py "fig/*.gp" "fig/*.eps" "fig/*.pdf" Makefile README.md
cd "../.."
