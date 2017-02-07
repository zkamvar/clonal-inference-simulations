# Karl Broman saves the day: http://kbroman.org/minimal_make/#automatic-variables
# If you simply just use $@ for the output_file, you're going to have a bad time
# because knitr will automatically assume that the directory in the output file
# is supposed to be nested inside of the input directory. Using the $(@F) and
# $(@D) variables prevents this unfortunate behavior!
reports/ROC_Curve.html: analysis/ROC_Curve.Rmd
	R --slave -e "rmarkdown::render(input = '$^', output_file = '$(@F)', output_dir = '$(@D)')"
