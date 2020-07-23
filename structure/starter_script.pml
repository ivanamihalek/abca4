load abca4_swissmodel_w_nbd_substrate.pdb

bg_color white

select ATP, hetatm
show spheres, ATP
color magenta, ATP

select tmd1, resi 1-45 or resi 645-870
select tmd2, resi 1340-1395 or resi 1665-1905
select Rdomain, resi 1141-1271 or resi 2161-2260
select nbd1, resi 960-1140
select nbd2, resi 1940-2160
select ecd1, resi 50-335 or resi 365-640
select ecd2, resi 1405-1660

color blue, tmd1

deselect




# this rotates the view, you can ignore it otherwise
set_view (\
     0.947501481,   -0.089684971,    0.306905657,\
    -0.299787045,    0.084594771,    0.950247645,\
    -0.111185089,   -0.992369115,    0.053267162,\
    -0.000337111,    0.000172377, -585.428283691,\
   143.356201172,  153.000122070,  123.770782471,\
   428.610961914,  742.422241211,  -20.000000000 )
