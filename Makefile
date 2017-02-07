reports/ROC_Curve.html: analysis/ROC_Curve.Rmd
	R --slave -e "rmarkdown::render(input = '$^', output_file = '$@', knit_root_dir = '.')"
