#in tables/
table_converter.py
## write 'include' as metadata keyword to search on
## write 'title' as metadata keyword as catalog title

python LBV_db.py tables/*.fits --master tables/MasseyXL.fit
python LBV_query.py catalog.json
python LBV_photo.py catalog.json -2MASS dMASStables -WISE dWISEtables
