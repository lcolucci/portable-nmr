# Portable NMR


# Analysis Code 

Here is a summary of the analysis scripts (all are written in Matlab unless otherwise specified): 

## Bioimpedance
1. **process_bioimpedance** - this script takes raw BI data (.xlsx), discards bad data, and saves two tidy master BI tables (one containing data from each individual measurement and one containing averaged data from the multiple measurements of each time point). 
2. **analysis_bi** - this script takes one of the tidy master BI tables generated from the previous script (process_bioimpedance) and facilitates data exploration and visualization. This script makes it easy to generate (1) scatter plots of BI data vs demographics data, (2) boxplots comparing BI values from HC and HD participants, and (3) a summary table of BI results (including the pre-to-post changes). NOTE: The boxplots generated in this script can be seen in fig 7 of the STM paper.  

## MR Sensors
1. **process_mrSensor_rawData** - this script takes raw NMR sensor data (.csv), averages multiple measurements from the same time point together, and saves the processed sensor data as a .mat file
2. **process_t2decays_table** - this script takes the processed sensor data (.mat) from the previous script and aggregates all the data into a single table (mrSensor_T2Decays_table.mat). 
3. **analyze_forced3expFits** - this script loads the previously generated table (mrSensor_T2Decays_table.mat), applies a constrained 3-exponential fit to each T2 decay (described in the paper), and saves the result in another table (legSensor_forced3exp_40ms_250ms.mat). An abbreviated version of this table (mrSensor_R2.csv) is saved for statistics to be performed in R (statistics_for_paper.R). 
4. **visualize_MRsensor_results** - scatter plot of NMR sensor results (same as fig. 3A but for NMR sensor data)

## MRI
1. **process_mri** - this script takes raw NiFTI MRI images, performs multiple multi-exponential fits on each pixel, and saves the results (e.g. mri_pixelByPixel_[n]exp.mat) as well as the T2-decay for each pixel (mri_t2PixelDecays.mat). 
1. **process_noise** - this script calculates the noise (std dev of the background) of each MRI scan and saves the result (mri_noise_summary.mat).
1. **merge_results_with_masks_snr** - this script takes the results of the previous scripts (mri_pixelByPixel_[n]exp.mat) and adds columns for each ROI that have a 0/1 indicating whether or not the pixel belongs to that ROI. This enables anlyses to be done at an ROI level. This script can also calculate SNR for each pixel. 
2. **process_pixelbypixel_averages** - this script processes raw pixelwise results from the previous script into summary tables for both pre and post values as well as the pre-to-post change in values (mri_pixelByPixel_[n]exp_summary_all.mat, mri_pixelByPixel_[n]exp_summary_AMPMchange.mat) 
2. **process_pixelData_to_ROImeans** - this script takes the raw T2 decays for each pixel (mri_t2PixelDecays.mat) and averages them together to produce an average T2 decay for each ROI (ROImeans_table.mat)
3. **process_ROImeans** - this script takes the average T2 decay of each ROI (generated in previous script, ROImeans_table.mat) and fits them to multi-exponential models (ROIs_2expFits.mat) 
3. **process_ROImeans_to_GroupMeans** - This script loads ROImeans_table.mat and averages decays together to produce an average HC pre, HC post, HD pre, and HD post decay for each ROI.
5. **process_mri_roi_averages** - this script loads ROIs_2expFits.mat and generates a summary table that includes pre-to-post changes and accomanying statistics (mri_roi_2exp_summary_AMPMchange.mat). This resulting table is used to create figs. 1D-F. 
6. **process_mask_sizes** - this script calculates the number of pixels of each mask of each MRI scan (mask_sizes.mat) 
7. **process_subcu_thickness** - this script loads the excel sheet with subcu/skin thickness measurements and processes it into (1) a tidy master table and (2) a tidy summary table for future analyses. 

## Paper
1. **fig6AC_tableS8_mrSensor** - boxplots of NMR sensor data for figs 6A,C and table summary of NMR sensor data for table S8
2. **fig6BD_mrSensor_vs_thickness_vs_BI** - creates figs 6B - Change in NMR Sensor R2 vs Change in BI Re - and 6D - NMR Sensor RA_3 vs Subcu. Thickness
3. **figS13_raw_mrSensor** - plots representative raw data from an HC and HD participant
4. **statistics_for_paper** - R script to implement the various statistical tests that are displayed in the paper

