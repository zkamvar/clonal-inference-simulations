SSR_DATA    := data/processed_results.rda data/full_results.rda
MA_SSR_DATA := data/ma_processed_results.rda data/ma_full_results.rda
TARGETS     := $(SSR_DATA) \
               $(MA_SSR_DATA) \
               data/genomic_data.rda \
               data/ROC_data.rda \
               data/ma_ROC_data.rda \
               reports/ROC_Curve.html \
               reports/ma_ROC_Curve.html \
               reports/equilibrium_analysis.html \
               reports/genomic_data_analysis.html \
               reports/post_analysis.html \
               reports/ma_post_analysis.html \
               reports/negative_ia.html \
               reports/jackknife_analysis.html \
               reports/ma_jackknife_analysis.html

.PHONY: all

all : $(TARGETS) manuscript/clonal-inference.pdf

$(SSR_DATA) : reports/ssr_data_cleaning.html
reports/ssr_data_cleaning.html : data/rda_files/ \
                                 data/diversity_rda_files/ \
                                 data/locus_rda_files/ \
                                 data/locus_contribution_rda_files/

$(MA_SSR_DATA) : reports/ma_ssr_data_cleaning.html
reports/ma_ssr_data_cleaning.html : data/ma_rda_files/ \
                                    data/ma_diversity_rda_files/ \
                                    data/ma_locus_rda_files/ \
                                    data/ma_locus_contribution_rda_files/

data/genomic_data.rda : reports/genomic_data_processing.html
reports/genomic_data_processing.html : data/genomic_rda_files/

data/ROC_data.rda : reports/ROC_Calculation.html \
                    $(SSR_DATA) \
                    data/genomic_data.rda

results/ROC_Curve.html : data/ROC_data.rda

data/ma_ROC_data.rda : reports/ma_ROC_Calculation.html \
                       $(MA_SSR_DATA)

results/ma_ROC_Curve.html : data/ma_ROC_data.rda

results/jackknife_analysis.html : $(SSR_DATA) \
                                  data/jack_rda_files/ \
                                  data/jack_psex_rda_files/

results/ma_jackknife_analysis.html : $(MA_SSR_DATA) \
                                     data/ma_jack_rda_files/ \
                                     data/ma_jack_psex_rda_files/

# Karl Broman saves the day: http://kbroman.org/minimal_make/#automatic-variables
# If you simply just use $@ for the output_file, you're going to have a bad time
# because knitr will automatically assume that the directory in the output file
# is supposed to be nested inside of the input directory. Using the $(@F) and
# $(@D) variables prevents this unfortunate behavior!

reports/%.html : analysis/%.Rmd
	-R --slave -e "rmarkdown::render(input = '$<', \
	               output_file = '$(@F)', \
	               output_dir = '$(@D)')"

manuscript/%.pdf : manuscript/%.Rmd $(TARGETS)
	-R --slave -e "rmarkdown::render(input = '$<')"
