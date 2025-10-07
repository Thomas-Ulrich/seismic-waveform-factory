# Minimal Makefile for Sphinx documentation

SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: help clean html

help:
	@echo "    make html      to build HTML docs"
	@echo "    make clean     to remove build dir"

clean:
	rm -rf $(BUILDDIR)

autogen:
	python ../tools/generate_ini_doc_from_scripts.py ../scripts/config_schema.py .

html: autogen
	$(SPHINXBUILD) -b html $(SPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/html

# You can add more targets for other formats (latexpdf, epub, etc.)
