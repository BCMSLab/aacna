R_OPT=--vanilla

SCRIPTS=scripts
FIG_SRC=scripts/figures
TAB_SRC=scripts/tables

MANUSCRIPT=manuscript
FIG_DIR=manuscript/figures
TAB_DIR=manuscript/tables

DATA=data
LOG=log

all: install \
	dir_tree \
	$(DATA)/wgcna.rda \
	$(FIG_DIR)/differentiation_markers.png \
	$(FIG_DIR)/expression_averages.png \
	$(FIG_DIR)/gene_membership.png \
	$(FIG_DIR)/module_cor_overrep.png \
	$(FIG_DIR)/module_networks.png \
	$(FIG_DIR)/overrep_by_modules.png \
	$(FIG_DIR)/preservation_ranks.png \
	$(FIG_DIR)/preservation_zsummary.png \
	$(FIG_DIR)/qc_eda.png \
	$(FIG_DIR)/sft_power.png \
	$(FIG_DIR)/similarity_centrality.png \
	$(FIG_DIR)/workflow.png \
	$(FIG_DIR)/workflow.png \
	$(TAB_DIR)/datasets.tex \
	$(TAB_DIR)/go_genes.tex \
	$(TAB_DIR)/module_interactions.tex \
	$(TAB_DIR)/module_interactions_term.tex \
	$(TAB_DIR)/module_members.tex \
	$(TAB_DIR)/network_hubs.tex \
	$(TAB_DIR)/time_points.tex \
	clean

# install/update the functions
install:
	Rscript -e "devtools::install('..')"

#show:
#	for i in `ls $(FIG_SRC) | cut -d '.' -f1`; do echo $(FIG_DIR)/$$i.png; done
#	for i in `ls $(TAB_SRC) | cut -d '.' -f1`; do echo $(TAB_DIR)/$$i.tex; done

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

# generate figures
$(FIG_DIR)/%.png: $(FIG_SRC)/%.R $(DATA)/wgcna.rda
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# generate tables
$(TAB_DIR)/%.tex: $(TAB_SRC)/%.R $(DATA)/wgcna.rda
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# generate workflow
$(FIG_DIR)/workflow.png: $(FIG_SRC)/workflow.dot
	dot -Tpng $(FIG_SRC)/workflow.dot -o $(FIG_DIR)/workflow.png

# clean up
clean:
	rm -f *.pdf
	rm -f *.RData
