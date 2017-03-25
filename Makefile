.PHONY: all
all: reports/genomic_data_processing.html \
     reports/ma_ssr_data_cleaning.html \
     reports/ssr_data_cleaning.html \
     reports/ROC_Calculation.html \
     reports/ROC_Curve.html \
     reports/equilibrium_analysis.html \
     reports/genomic_data_analysis.html \
     reports/ma_ROC_Calculation.html \
     reports/ma_ROC_Curve.html \
     reports/ma_post_analysis.html \
     reports/negative_ia.html \
     reports/post_analysis.html
# Karl Broman saves the day: http://kbroman.org/minimal_make/#automatic-variables
# If you simply just use $@ for the output_file, you're going to have a bad time
# because knitr will automatically assume that the directory in the output file
# is supposed to be nested inside of the input directory. Using the $(@F) and
# $(@D) variables prevents this unfortunate behavior!
reports/%.html: analysis/%.Rmd \
	data/genomic_rda_files/ \
	data/rda_files/ \
	data/diversity_rda_files \
	data/locus_rda_files \
	data/locus_contribution_rda_files \
	data/jack_rda_files/ \
	data/ma_rda_files/ \
	data/ma_diversity_rda_files \
	data/ma_locus_rda_files \
	data/ma_locus_contribution_rda_files \
	data/ma_jack_rda_files/
	R --slave -e "rmarkdown::render(input = '$<', \
	              output_file = '$(@F)', \
	              output_dir = '$(@D)')"
