from transitleastsquares import catalog_info

tic_id = int('459350219') #466265409
ab, mass, mass_min, mass_max, radius, radius_min, radius_max = catalog_info(TIC_ID=tic_id)

print(ab, mass, mass_min, mass_max, radius, radius_min, radius_max)