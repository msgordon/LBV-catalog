#in tables/
table_converter.py
## write 'include' as metadata keyword to search on
## write 'title' as metadata keyword as catalog title

python LBV_db.py tables/*.fits --master tables/MasseyXL.fit
python LBV_query.py catalog.json
python LBV_photo.py catalog.json -2MASS dMASStables -WISE dWISEtables
python LBV_photo_corr.py

python LBV_SEDplot.py photometry.fits
#python LBV_html.py photometry.fits -o LBV_M31_photo.html

python LBV_specmatch.py
   -> generates .list files of IDs

python LBV_SEDspec.py photometry.tsv *.list
