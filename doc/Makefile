prj=invt

ver=8.9

tikzfigs=fig/massq.pdf fig/vardep.pdf
figgps=$(wildcard fig/*.gp)
pdffigs=$(patsubst %.gp,%.pdf, $(figgps)) $(tikzfigs)
epsfigs=$(patsubst %.gp,%.eps, $(figgps))

texpdf = $(prj).pdf $(prj)_notes.pdf

alltexpdf = $(patsubst %.tex,%.pdf, $(wildcard $(prj)*.tex))

regular: $(texpdf)

all: $(alltexpdf)

$(alltexpdf) : %.pdf : %.tex simul.bib $(pdffigs) Makefile
	make -C fig
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<

$(pdffigs): $(epsfigs) $(figgps) fig/figs.tex
	make -C fig

clean:
	$(RM) -rf *.aux *.bbl *.blg *.dvi *.log *Notes.bib *~
	rstrip.py -Rlv

zip: $(prj)docV$(ver).zip

# make a zip for the current version

$(prj)docV$(ver).zip::
	mkdir -p old/V$(ver)
	rsync -avz $(prj).tex old/V$(ver)/$(prj)_annotatedV$(ver).tex
	rsync -avz Makefile old/V$(ver)/
	rsync -avz README.md old/V$(ver)/
	rsync -avz deannote.py old/V$(ver)/
	rsync -avz rmbiburl.py old/V$(ver)/
	rsync -avz checkeq.py old/V$(ver)/
	python rmbiburl.py -i simul.bib -o old/V$(ver)/simul.bib
	cd old/V$(ver) && \
	  ln -sf ../../fig .
	cd old/V$(ver) && \
	  python deannote.py -v -c -i $(prj)_annotatedV$(ver).tex -o $(prj)V$(ver).tex
	make -C old/V$(ver)/ $(prj)V$(ver).pdf
	make -C old/V$(ver)/ $(prj)_annotatedV$(ver).pdf
	make diff
	cd old/V$(ver) && \
	  zip -r ../../$(prj)docV$(ver).zip *.tex *.pdf simul.bib *.py \
	       $(epsfigs) $(pdffigs) $(figgps) \
	       fig/figs.tex fig/Makefile Makefile README.md

#	rsync -avz $(prj)_notes.tex old/V$(ver)/$(prj)_notesV$(ver).tex


$(prj)doc.zip: $(prj).tex $(prj).pdf \
  	       $(prj)_notes.tex $(prj)_notes.pdf \
  	       simul.bib \
	       $(epsfigs) $(pdffigs) $(figgps) \
	       Makefile README.md
	zip -r $@ $^

refver=8.8

diff::
	mkdir -p old/V$(ver)
	#module load perl
	# consider the flag --math-markup=0 if there is trouble for formulae
	cd old/V$(ver) && latexdiff \
	  --exclude-textcmd="footnote" \
	  --exclude-textcmd="bibinfo" \
	  --append-safecmd="bibinfo" \
	  --config="PICTUREENV=(?:picture|DIFnomarkup|align|equation|tikzpicture|table|tabular|enumerate)[\w\d*@]*" \
	  ../V$(refver)/$(prj)_annotatedV$(refver).tex $(prj)_annotatedV$(ver).tex > $(prj)_diffV$(ver).tex
	make -C old/V$(ver) $(prj)_diffV$(ver).pdf

refver2=4.1_jd

diff2::
	mkdir -p old/V$(ver)
	cd old/V$(ver) && latexdiff \
	  --config="PICTUREENV=(?:picture|DIFnomarkup|displaymath|align|equation|tikzpicture|table|tabular|enumerate)[\w\d*@]*" \
	  ../V$(refver2)/invtV$(refver2).tex invt_annotatedV$(ver).tex > invt_diffV$(ver)_$(refver2).tex
	make -C old/V$(ver) invt_diffV$(ver)_$(refver2).pdf

Bossman: zip
	cp $(prj)docV$(ver).zip /Bossman/cz1/invt/doc/
	#rsync -avL $(prj)docV$(ver).zip cz1@129.109.88.204:/Bossman/cz1/invt/doc/
