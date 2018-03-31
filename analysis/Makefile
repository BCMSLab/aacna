all: analysis tables figures

analysis:
	Rscript scripts/01.analysis.R 2> log.analysis
tables:
	Rscript scripts/02.tables.R 2> log.tables
figures:
	Rscript scripts/03.figures.R 2> log.figures

