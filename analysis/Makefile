R_OPT=--vanilla

SCRIPTS=scripts
FIG_SRC=scripts/figures
TAB_SRC=scripts/tables

MANUSCRIPT=manuscript
FIG_DIR=manuscript/figures
TAB_DIR=manuscript/tables

DATA=data
LOG=log

all: install dir_tree $(DATA)/wgcna.rda $(TAB_DIR)/* $(FIG_DIR)/* clean

# install/update the functions
install:
	Rscript -e "devtools::install('..')"

# make directory tree
dir_tree:
	test ! -d $(LOG) && mkdir $(LOG) || echo '$(LOG)/ is already there.'
	test ! -d $(DATA) && mkdir $(DATA) || echo '$(DATA)/ is already there.'
	test ! -d $(MANUSCRIPT) && mkdir $(MANUSCRIPT) || echo '$(MANUSCRIPT)/ is already there.'
	test ! -d $(TAB_DIR) && mkdir $(TAB_DIR) || echo '$(TAB_DIR)/ is already there.'
	test ! -d $(FIG_DIR) && mkdir $(FIG_DIR) || echo '$(FIG_DIR)/ is already there.'

# run the main part of the analysis
$(DATA)/wgcna.rda: $(SCRIPTS)/01.analysis.R
	R CMD BATCH $(R_OPT) $(SCRIPTS)/01.analysis.R $(LOG)/01.analysis.R.Rout

# generate tables
$(TAB_DIR)/%.tex: $(TAB_SRC)/%.R $(DATA)/wgcna.rda
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# generate figures
$(FIG_DIR)/%.png: $(FIG_SRC)/%.R $(DATA)/wgcna.rda
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# clean up
clean:
	test -f *.pdf && rm *.pdf || echo 'All clean.'
	test -f *.pdf && rm *.RData || echo 'All clean.'
