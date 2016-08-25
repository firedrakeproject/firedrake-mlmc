all:

docs:
	@echo "    generating firedrake-mlmc documentation"
	@pdflatex -interaction batchmode firedrake_mlmc_doc
	@bibtex firedrake_mlmc_doc
	@pdflatex -interaction batchmode firedrake_mlmc_doc
	@pdflatex -interaction batchmode firedrake_mlmc_doc
