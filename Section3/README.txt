This package creates four individual figures.

Input files:
- paper_figure_data/bulk_modulus.csv
- paper_figure_data/volumetric_stress_strain.csv
- paper_figure_data/internal_variables.csv
- paper_figure_data/figure4b_300K.csv
- paper_figure_data/figure4b_400K.csv
- paper_figure_data/figure4b_500K.csv

Output files:
- paper_figures_output/bulk_modulus.pdf
- paper_figures_output/bulk_modulus.png
- paper_figures_output/stress_strain_plot.pdf
- paper_figures_output/stress_strain_plot.png
- paper_figures_output/internal_variables.pdf
- paper_figures_output/internal_variables.png
- paper_figures_output/relative_error_epoch.pdf
- paper_figures_output/relative_error_epoch.png

Figure 4(a) and Figure 4(b) both use "Relative error" with no percent sign.
LaTeX rendering is ON by default. For testing without LaTeX:
PAPER_USE_TEX=0 python plot_paper_figures_3_4.py
