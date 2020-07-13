structure_home = "/home/ivana/Dropbox/Anne/abca4/structure"

pos_in_STGD1 = [1038,1488,818,407,765,1805,971,2077,607,1430,1408,1562,244,1843,68,1433,1898,1733,2263,1486,550,974,
				2071,2139,716,54,220,1022,1729,2035,608,60,821,58,1890,249,1693,602,1526,863,1300,1072,1532,931,1087,
				1108,541,643,340,1129,2038,537,1406,572,1537,1063,1780,965,767,1703,18,1975,943,333,1820,764,445,959,
				1399,212,1689,686,247,957,1885,2096,1388,797,2229,1748,2106,1071,336,636,1112,1631,1763,1838,2150,230,
				75,1031,935,972,2160,854,72,549,2050,1036,156,300,1443,77,873,1380,1439,1440,851,1799,989,653,1794,
				1055,1250,471,1736,1884,2027,1705,824,897,2241,1652,978,1513,1525,1961,1019,1014,1977,1696,2107,206,
				328,1940,380,2030,65,645,1204,100,2149,192,1201,190,1776,635,2131,1640,1896,1122,1490,1886,525,849,1097,
				96,1429,2128,523,309,24,1681,1682,1683,1684,1685,13,14,15,1761,1762,1763]
other_phenotype_positions = {
				 "FFM" :    [1038,863,1488,1961,1108,541,1091,1562,1971,1940,2030,1970,943,11,1640,1508,339,2027,991,1253,2106],
				 "ARMD2" :  [762,2137,1970,1129,1517,818,1562,1898,2047,1724,471,1578,1428,1977,2177],
				 "CORD3" :  [1038,863,2060,1122,2150,1562,407,1490,2146,1640,212,1598,65,541],
				 "gnomad100"  : [423,2255,212,1868,1948,943,2177,915,915,1288,1288],
				 "gnomad1000" : [1201,1300,1961,1428,863,1209,1970,2050,552,152,1591,537,2107,901,1038,643,1433,643,
				                  1898,206,156, 849,1642,1562,897,106,1501]}
patients = [1038,863,989,2027,541,1526,54,1556,75,1087,2096,602,821,1776,1250,1961,418,1380,2035,1490,1868,1443,803]

other_phenotype_color = {"FFM" : "slate",  "ARMD2" : "pink", "CORD3" :  "sand",
						"gnomad100"   : "gray", "gnomad1000" : "palegreen"}

region_range = {"tmd1":[[1, 45], [645, 870]] , "tmd2":[[1340, 1395], [1665, 1905]],
				"Rdomain":[[1141, 1271], [2161, 2260]],
				"nbd1":[[960, 1140]]  , "nbd2":[[1940, 2160]],
				"ecd1":[[50,  335], [365,  640]], "ecd2":[[1405, 1660]]}


def find_region(res):
	rmin = 1000
	rmax = -1
	for name, ranges in region_range.items():
		for [rfrom, rto] in ranges:
			if rmin>rfrom: rmin = rfrom
			if rmax<rto: rmax = rto
			if rfrom<=res<=rto: return name
	if res<rmin: return "Nterm"
	if res>rmax: return "Cterm"
	return "linker"

region_range_pymol = {"tmd1":["resi 1-45","resi 645-870"], "tmd2":["resi 1340-1395","resi 1665-1905"],
				"Rdomain":["resi 1141-1271","resi 2161-2260"],
				"nbd1":["resi 960-1140"], "nbd2":["resi 1940-2160"],
				"ecd1":["resi 50-335","resi 365-640"], "ecd2":["resi 1405-1660"]}

region_color = {"tmd1":"blue", "tmd2":"green", "Rdomain":"grey",
				"nbd1":"yellow", "nbd2":"magenta", "ecd1":"cyan", "ecd2":"orange"}

home_view = "  -0.002808416,   -0.999822557,    0.018627511,\
     0.144136816,    0.018028343,    0.989393473,\
    -0.989553809,    0.005463540,    0.144060597,\
     0.000000000,    0.000000000, -671.049865723,\
   129.639907837,  134.362030029,  141.107208252,\
   529.060791016,  813.038940430,  -20.000000000 "

nbd_zoom = "   -0.022720490,   -0.999729455,    0.004949572,\
    -0.114353731,    0.007517095,    0.993412316,\
    -0.993180454,    0.022004750,   -0.114493623,\
    -0.000059098,    0.000020376, -235.058837891,\
   180.633300781,  138.383773804,  141.106155396,\
   109.819557190,  360.313537598,  -20.000000000 "

tm_zoom = "    -0.022720490,   -0.999729455,    0.004949572,\
    -0.114353731,    0.007517095,    0.993412316,\
    -0.993180454,    0.022004750,   -0.114493623,\
    -0.000082755,    0.000042200, -230.380935669,\
   148.598373413,  138.746856689,  142.714340210,\
   105.159667969,  355.653656006,  -20.000000000  "

ecd_zoom = "    -0.022720490,   -0.999729455,    0.004949572,\
    -0.114353731,    0.007517095,    0.993412316,\
    -0.993180454,    0.022004750,   -0.114493623,\
    -0.000066902,    0.000004679, -332.334716797,\
    75.022254944,  135.463134766,  143.767776489,\
   207.126342773,  457.620208740,  -20.000000000 "


def find_nbd_res(region_range):
	nbd_residues = []
	for domain in ["nbd1", "nbd2", "Rdomain"]:
		for region_str in region_range[domain]:
			[frm,to] = [int(i) for i in region_str.replace("resi ", "").split("-")]
			nbd_residues.extend(list(range(frm,to+1)))
	return set(nbd_residues)

