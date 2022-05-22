### Instructions ###
# 1. Copy this script to the directory where you want to store
#    the supporting data files.
# 2. In your ~/.bashrc or ~/.bash_profile file, point the environment
#    variable SOC_ROMAN_FILES to this directory. Then copy the following
#    environment variables into the same file:
#       export WEBBPSF_PATH="${SOC_ROMAN_FILES}/webbpsf-data"
#       export PYSYN_CDBS="${SOC_ROMAN_FILES}/grp/redcat/trds"
#       export pandeia_refdata="${SOC_ROMAN_FILES}/pandeia_data-1.7_roman"
# 3. Run the script: bash retrieve_data.sh

### WebbPSF Supporting Data Files ###
curl -OL https://stsci.box.com/shared/static/34o0keicz2iujyilg4uz617va46ks6u9.gz
tar zvxf 34o0keicz2iujyilg4uz617va46ks6u9.gz
rm 34o0keicz2iujyilg4uz617va46ks6u9.gz

### Synphot Models and Atlases ###
curl -OL https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_everything_multi_v11_sed.tar
curl -OL https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_star-galaxy-models_multi_v3_synphot2.tar
curl -OL https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_pheonix-models_multi_v2_synphot5.tar
curl -OL https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_jwst_multi_etc-models_multi_v1_synphot7.tar

tar xvzf hlsp_reference-atlases_hst_multi_everything_multi_v11_sed.tar
tar xvzf hlsp_reference-atlases_hst_multi_star-galaxy-models_multi_v3_synphot2.tar
tar xvzf hlsp_reference-atlases_hst_multi_pheonix-models_multi_v2_synphot5.tar
tar xvzf hlsp_reference-atlases_jwst_multi_etc-models_multi_v1_synphot7.tar

rm hlsp_reference-atlases_hst_multi_everything_multi_v11_sed.tar
rm hlsp_reference-atlases_hst_multi_star-galaxy-models_multi_v3_synphot2.tar
rm hlsp_reference-atlases_hst_multi_pheonix-models_multi_v2_synphot5.tar
rm hlsp_reference-atlases_jwst_multi_etc-models_multi_v1_synphot7.tar

curl -OL https://archive.stsci.edu/hlsps/reference-atlases/cdbs/current_calspec/alpha_lyr_stis_010.fits
mkdir grp/redcat/trds/calspec
mv alpha_lyr_stis_010.fits grp/redcat/trds/calspec

### Pandeia Data Files ###
curl -OL https://stsci.box.com/shared/static/ycbm34uxhzafgb7te74vyl2emnr1mdty.gz
tar xvzf ycbm34uxhzafgb7te74vyl2emnr1mdty.gz
rm ycbm34uxhzafgb7te74vyl2emnr1mdty.gz
