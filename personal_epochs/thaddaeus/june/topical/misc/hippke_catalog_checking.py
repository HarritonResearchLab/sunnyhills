from transitleastsquares import catalog_info

tic_id = int('459350219') #466265409
ab, mass, mass_min, mass_max, radius, radius_min, radius_max = catalog_info(TIC_ID=tic_id)

radius_max = radius+radius_max 
radius_min = radius-radius_min 
mass_min = mass-mass_min 
mass_max = mass+mass_max 

print(ab, mass, mass_min, mass_max, radius, radius_min, radius_max)