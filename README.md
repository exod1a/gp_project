
# Code for the paper: "Reconstructing ecological community dynamics from limited observations"

The preprint of the article is published [here](https://arxiv.org/abs/2501.03820). In brief, it is a Bayesian model implemented in Stan that looks at decomposing the dynamics of a system into its constituent deterministic and stochastic parts (drift and diffusion respectively) using Gaussian process priors and the Euler-Maruyama approximation. The model can output a posterior on the number of stable states in the system, the locations of the stable and unstable points, the stationary density, stability landscape, and exit time if alternative stable states were found. 

# Main Figures
Before running any other file, run "scripts/main.R". It will load the necessary packages, functions, and compile the Stan model.

## Figure 1

Mostly made in Inkscape but Figure 1e can be made by running "scripts/figures/main_text/short_vs_long_ts.R"

<img src="./output/figures/main_text/figure 1/short_vs_long_ts_plot.pdf" width="350" height="200">
<img src="./output/figures/main_text/figure 1/densities.pdf" width="200" height="200">

<img src="./output/figures/main_text/figure 1/fig1.pdf" width="550" height="400">

## Figure 2

Made entirely in Inkscape

<img src="./output/figures/main_text/figure 2/fig2.pdf" width="350" height="350">

## Figure 3

Run "scripts/figures/main_text/fig3_code.R"

<img src="./output/figures/main_text/figure 3/fig3.pdf" width="400" height="400">

The final figure used in the publication was edited in Inkscape. The code to create the heatmaps (Fig. 3f,l) was run on the computing cluster Puhti and is therefore not included here.

## Figure 4

Run "scripts/figures/main_text/fig4_code.R"

<img src="./output/figures/main_text/figure 4/fig4.pdf" width="700" height="350">

The final figure used in the publication was edited in Inkscape

## Figure 5

Run "scripts/figures/main_text/fig5_code.R"

<img src="./output/figures/main_text/figure 5/fig5.pdf" width="250" height="350">

The final figure used in the publication was edited in Inkscape

# Extended Data Figures

## Figure 1

Run "scripts/figures/extended_data/ext_data_fig1code.R"

<img src="./output/figures/extended_data/figure 1/fig1.pdf" width="350" height="350">

The final figure used in the publication was edited in Inkscape

## Figure 2

Made entirely in Inkscape

<img src="./output/figures/extended_data/figure 2/fig2.pdf" width="700" height="350">

## Figure 3

Run "scripts/figures/extended_data/ext_data_fig3code.R"

<img src="./output/figures/extended_data/figure 3/fig3.pdf" width="500" height="350">

The final figure used in the publication was edited in Inkscape

## Figure 4

Run "scripts/figures/extended_data/ext_data_fig4code.R"

<img src="./output/figures/extended_data/figure 4/fig4.pdf" width="500" height="350">

The final figure used in the publication was edited in Inkscape

## Figure 5

Run "scripts/figures/extended_data/ext_data_fig5code.R"

<img src="./output/figures/extended_data/figure 5/fig5.pdf" width="250" height="350">

The final figure used in the publication was edited in Inkscape

# Supplementary Video

Run "scripts/figures/supplementary video/short_vs_long_ts_video.R"

<img src="./output/videos/ts_comparison_fast_high_res.gif" width="250" height="350">
