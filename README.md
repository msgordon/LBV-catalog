# 07/15/14
python photo_plot.py Drout_list_ZOMG.fits -pdf BB_M31_Drout.pdf --gal M31
python photo_plot.py Drout_list_ZOMG.fits -pdf BB_M33_Drout.pdf --gal M33

# 07/14/14
ran Drout_list.py on Drout_list_init.tsv
  ->Drout_list_FINAL.fits/tsv
ran final_check.py on Drout_list_FINAL.fits to see if any remaining matches in Mould/McQuinn IRAS photometry
   matches to 2MASS things

wrote mag2flux.py

# uses: http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4c.html#wappco   for aperture correction
run final_check.py
  ->  GENERATES Drout_list_ZOMG.tsv/fits


## July 7, 2014
python photo_plot.py photometry_ZOMG.tsv -pdf allSED.pdf
python photo_plot.py photometry_ZOMG.tsv -pdf M31A_SED.pdf --idlist M31A_IDs.list
python photo_plot.py photometry_ZOMG.tsv -pdf M31B_SED.pdf --idlist M31B_IDs.list



#in tables/
table_converter.py
## write 'include' as metadata keyword to search on
## write 'title' as metadata keyword as catalog title

python LBV_db.py tables/*.fits --master tables/MasseyXL.fit
python LBV_query.py catalog.json
python LBV_photo.py catalog.json -2MASS dMASStables -WISE dWISEtables
#python LBV_photo_corr.py

python LBV_SEDplot.py photometry.fits
#python LBV_html.py photometry.fits -o LBV_M31_photo.html

python LBV_specmatch.py
   -> generates .list files of IDs


# D 4304 is missing? in M31A_Blue

##M31
see readme in m312013
python LBV_SEDspec.py photometry_corr_FINAL.tsv M31A_FINAL.list -pdf BB_spec_M31A.pdf
python LBV_SEDspec.py photometry_corr_FINAL.tsv M31B_FINAL.list -pdf BB_spec_M31B.pdf

---> 
stitch_spec.py --list M31A_FINAL.list --outdir M312013/M31A/
stitch_spec.py --list M31B_FINAL.list --outdir M312013/M31B/

